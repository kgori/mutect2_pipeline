nextflow.enable.dsl=2

params.reference               = "genome.fa"
params.normals                 = "normals_folder"
params.tumours                 = "tumours_folder"
params.outdir                  = "results"
params.fineSplitNumIntervals   = 8
params.coarseSplitNumIntervals = 4
params.intervals               = "$projectDir/NO_FILE"
params.realign                 = false

def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}

include { makeBamToSampleNameMap }                         from "./preprocessSamples.nf"
include { indexReference }                                 from "./preprocessReference.nf"
include { splitIntervals }                                 from "./preprocessReference.nf"
include { makeReferenceDict }                              from "./preprocessReference.nf"
include { runMutectOnNormal }                              from "./initialVariantCalling.nf"
include { runMutectOnTumour }                              from "./initialVariantCalling.nf"
include { runHaplotypeCallerOnNormal }                     from "./initialVariantCalling.nf"
include { runPlatypus }                                    from "./initialVariantCalling.nf"
include { runFreeBayes }                                   from "./initialVariantCalling.nf"
include { concatMutectVcfParts as concatMutectTumour }     from "./vcfConcatenator.nf"
include { concatMutectVcfParts as concatMutectNormal }     from "./vcfConcatenator.nf"
include { concatVcfParts as concatGvcfs }                  from "./vcfConcatenator.nf"
include { concatVcfParts as concatPlatypus }               from "./vcfConcatenator.nf"
include { concatVcfParts as concatFreeBayes }              from "./vcfConcatenator.nf"

process realign {
    input:
    tuple val(tumour_normal_tag), val(sample_name), path(reference_files), path(sample_bam)

    output:
    tuple val(tumour_normal_tag), val(sample_name), path("*.remapped.cram*")

    publishDir "${params.outdir}/RealignedBams", mode: 'symlink'

    script:
    """
    realign.sh \
        -o ${sample_name}.remapped.cram \
        --output-format CRAM \
        --threads ${task.cpus} \
        --write-index \
        ${sample_bam[0]} \
        ${reference_files[0]}
    """
}

workflow {
    /////////////////////////////////
    //         Preanalysis setup
    //

    // Load and process the reference
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true)
    fai_ch = indexReference(ref_ch)
    dict_ch = makeReferenceDict(ref_ch)
    ref_files = ref_ch.merge(fai_ch).merge(dict_ch)

    // Chop into intervals for scattering-gathering
    intervals_ch = Channel.fromPath(params.intervals)
    ivls = splitIntervals(
        ref_files,
        params.fineSplitNumIntervals,
        params.coarseSplitNumIntervals,
        intervals_ch)

    fineIvls_ch = ivls.fineIntervals.flatten()
        .map { interval ->
            def intervalNumberMatch = interval.getName() =~ /^(\d+)/
                def intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
            return tuple(intervalNumber, interval) }
    
    coarseIvls_ch = ivls.coarseIntervals.flatten()
        .map { interval ->
            def intervalNumberMatch = interval.getName() =~ /^(\d+)/
                def intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
            return tuple(intervalNumber, interval) }

    // Load normal bams
    bam_bai_inputs   = Channel.fromFilePairs("${params.normals}/*.bam{,.bai}")
    bam_csi_inputs   = Channel.fromFilePairs("${params.normals}/*.bam{,.csi}")
    bam_inputs = remove_duplicate_filepair_keys(bam_csi_inputs, bam_bai_inputs)
    cram_inputs = Channel.fromFilePairs("${params.normals}/*.cram{,.crai}")
    normals = bam_inputs.concat(cram_inputs)

    // Load tumour bams
    bam_bai_inputs   = Channel.fromFilePairs("${params.tumours}/*.bam{,.bai}")
    bam_csi_inputs   = Channel.fromFilePairs("${params.tumours}/*.bam{,.csi}")
    bam_inputs = remove_duplicate_filepair_keys(bam_csi_inputs, bam_bai_inputs)
    cram_inputs = Channel.fromFilePairs("${params.tumours}/*.cram{,.crai}")
    tumours = bam_inputs.concat(cram_inputs)

    // Extract sample names from BAMs
    normal_name_map_ch = normals
        .map { label, bam -> bam[0] }
    tumour_name_map_ch = tumours
        .map { label, bam -> bam[0] }
    makeBamToSampleNameMap(normal_name_map_ch.concat(tumour_name_map_ch).collect())

    // If doing realignment, do it here
    if (params.realign) {
        bwa_index_ch = Channel.fromFilePairs(params.reference + "{,.amb,.ann,.bwt.2bit.64,.pac,.0123}", size: 6, checkIfExists: true)

        normals_realign_ch = normals.combine(bwa_index_ch)
            .map { samplename, bam, ref_label, bwa_index_files ->
                tuple("normal", samplename, bwa_index_files, bam) }

        tumours_realign_ch = tumours.combine(bwa_index_ch)
            .map { samplename, bam, ref_label, bwa_index_files ->
                tuple("tumour", samplename, bwa_index_files, bam) }

        realign_ch = normals_realign_ch.concat(tumours_realign_ch)

        realigned_ch = realign(realign_ch)

        normals = realigned_ch.filter { type, label, bam -> type == "normal" }
            .map { type, label, bam -> tuple(label, bam) }

        tumours = realigned_ch.filter { type, label, bam -> type == "tumour" }
         .map { type, label, bam -> tuple(label, bam) }
    }


    // Construct a channel for every combination of sample and interval,
    // separately for normals and tumours
    normal_fine_interval_ch   = normals.combine(fineIvls_ch)
    tumour_fine_interval_ch   = tumours.combine(fineIvls_ch)
    normal_coarse_interval_ch = normals.combine(coarseIvls_ch)
    tumour_coarse_interval_ch = tumours.combine(coarseIvls_ch)

    ///////////////////////////////////////////////////////
    //         Stage 1: First round of calling
    //

    // Run Mutect2 on normals
    normals_for_calling_ch = ref_files.combine(normal_fine_interval_ch)
        .map { fa, fai, dict, sample, normal_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, normal_bam, interval) }
    normals_called_by_mutect_ch =
        runMutectOnNormal(normals_for_calling_ch) // Channel: sample, vcf, tbi, stats

    // Run Mutect2 on tumours
    tumours_for_calling_ch = ref_files.combine(tumour_fine_interval_ch)
        .map { fa, fai, dict, sample, tumour_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, tumour_bam, interval) }
    tumours_first_called_by_mutect_ch =
        runMutectOnTumour(tumours_for_calling_ch) // Channel: sample, vcf, tbi, stats

    // Run HaplotypeCaller on normals
    gvcfs_ch = runHaplotypeCallerOnNormal(normals_for_calling_ch) // Channel: sample, gvcf, tbi

    // Prepare channels for alternative variant callers
    normals_for_alt_calling_ch = ref_files.combine(normal_coarse_interval_ch)
        .map { fa, fai, dict, sample, normal_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, normal_bam, interval) }
    tumours_for_alt_calling_ch = ref_files.combine(tumour_coarse_interval_ch)
        .map { fa, fai, dict, sample, tumour_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, tumour_bam, interval) }
    
    alt_calling_ch = normals_for_alt_calling_ch.mix(tumours_for_alt_calling_ch)

    // Run platypus on normals and tumours
    platypus_calls_ch = runPlatypus(alt_calling_ch) // Channel: sample, vcf, tbi

    // Run freebayes on normals and tumours
    freebayes_calls_ch = runFreeBayes(alt_calling_ch) // Channel: sample, vcf, tbi

    // Concatenate all variant calls - all call result channels start with tuple(sample, vcf, tbi...)
    grouped_mutect_tumours_ch = tumours_first_called_by_mutect_ch
        .map { it -> tuple(it[0], it[1], it[2]) }
        .groupTuple(size: params.fineSplitNumIntervals)
        .map { sample, vcfs, tbis -> tuple(sample,
                                           "mutect2_tumour",
                                           vcfs.sort { it.name },
                                           tbis.sort { it.name }) }
    grouped_mutect_tumours_with_ref_ch = ref_files.combine(grouped_mutect_tumours_ch)
        .map { fa, fai, dict, sample, label, vcfs, tbis ->
            tuple(sample, label, [fa, fai, dict], vcfs, tbis) }
    concatMutectTumour(grouped_mutect_tumours_with_ref_ch)
    
    grouped_mutect_normals_ch = normals_called_by_mutect_ch
        .map { it -> tuple(it[0], it[1], it[2]) }
        .groupTuple(size: params.fineSplitNumIntervals)
        .map { sample, vcfs, tbis -> tuple(sample,
                                           "mutect2_normal",
                                           vcfs.sort { it.name },
                                           tbis.sort { it.name }) }
    grouped_mutect_normals_with_ref_ch = ref_files.combine(grouped_mutect_normals_ch)
        .map { fa, fai, dict, sample, label, vcfs, tbis ->
            tuple(sample, label, [fa, fai, dict], vcfs, tbis) }
    concatMutectNormal(grouped_mutect_normals_with_ref_ch)

    grouped_haplotypecaller_ch = gvcfs_ch
        .map { it -> tuple(it[0], it[1], it[2]) }
        .groupTuple(size: params.fineSplitNumIntervals)
        .map { sample, vcfs, tbis -> tuple(sample,
                                           "haplotypecaller_normal",
                                           vcfs.sort { it.name },
                                           tbis.sort { it.name }) }
    grouped_haplotypecaller_with_ref_ch = ref_files.combine(grouped_haplotypecaller_ch)
        .map { fa, fai, dict, sample, label, vcfs, tbis ->
            tuple(sample, label, [fa, fai, dict], vcfs, tbis) }
    concatGvcfs(grouped_haplotypecaller_with_ref_ch)

    grouped_platypus_ch = platypus_calls_ch
        .map { it -> tuple(it[0], it[1], it[2]) }
        .groupTuple(size: params.coarseSplitNumIntervals)
        .map { sample, vcfs, tbis -> tuple(sample,
                                           "platypus",
                                           vcfs.sort { it.name },
                                           tbis.sort { it.name }) }
    grouped_platypus_with_ref_ch = ref_files.combine(grouped_platypus_ch)
        .map { fa, fai, dict, sample, label, vcfs, tbis ->
            tuple(sample, label, [fa, fai, dict], vcfs, tbis) }
    concatPlatypus(grouped_platypus_with_ref_ch)

    grouped_freebayes_ch = freebayes_calls_ch
        .map { it -> tuple(it[0], it[1], it[2]) }
        .groupTuple(size: params.coarseSplitNumIntervals)
        .map { sample, vcfs, tbis -> tuple(sample,
                                           "freebayes",
                                           vcfs.sort { it.name },
                                           tbis.sort { it.name }) }
    grouped_freebayes_with_ref_ch = ref_files.combine(grouped_freebayes_ch)
        .map { fa, fai, dict, sample, label, vcfs, tbis ->
            tuple(sample, label, [fa, fai, dict], vcfs, tbis) }
    concatFreeBayes(grouped_freebayes_with_ref_ch)
}
