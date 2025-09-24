nextflow.enable.dsl=2

params.reference    = "genome.fa"
params.normals      = "normals_folder"
params.tumours      = "tumours_folder"
params.outdir       = "results"
params.numIntervals = 8

def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}

include { makeBamToSampleNameMap }                             from "./preprocessSamples.nf"
include { indexReference }                                     from "./preprocessReference.nf"
include { splitIntervals }                                     from "./preprocessReference.nf"
include { makeReferenceDict }                                  from "./preprocessReference.nf"
include { runMutectOnNormal }                                  from "./variantCalling.nf"
include { runMutectOnTumour }                                  from "./variantCalling.nf"
include { runHaplotypeCallerOnNormal }                         from "./variantCalling.nf"
include { genomicsDBImport }                                   from "./genomicsDB.nf"
include { genomicsDBImport_PON }                               from "./genomicsDB.nf"
include { genomicsDBImport_Somatic }                           from "./genomicsDB.nf"
include { genotypeGVCFs }                                      from "./germlineResource.nf"
include { mergeGenotypedVCFs }                                 from "./germlineResource.nf"
include { createPanelOfNormals }                               from "./panelOfNormals.nf"
include { mergePonVCFs }                                       from "./panelOfNormals.nf"
include { extractSomaticCandidates }                           from "./somaticCandidates.nf"
include { normalizeSomaticCandidates }                         from "./somaticCandidates.nf"
include { mergeSomaticCandidates }                             from "./somaticCandidates.nf"
include { finalizeSomaticCandidates }                          from "./somaticCandidates.nf"
include { callSomaticVariants }                                from "./somaticCalling.nf"
include { callSomaticVariants as callSomaticVariantsOnNormal } from "./somaticCalling.nf"
include { recallGermlineVariants }                             from "./somaticCalling.nf"


process getPileupSummaries {
    input:
    tuple val(sample),
        path(reference),
        path(bam),
        path(germline_resource)

    output:
    tuple val(sample), path("${sample}.pileups")

    script:
    """
    gatk GetPileupSummaries \
        --reference ${reference[0]} \
        --input ${bam[0]} \
        --variant ${germline_resource[0]} \
        --intervals ${germline_resource[0]} \
        --minimum-population-allele-frequency 0.01 \
        --maximum-population-allele-frequency 0.5 \
        --output ${sample}.pileups
    """
}

process calculateContamination {
    input:
    tuple val(sample), path(pileup)

    output:
    tuple val(sample), path("${sample}.contamination.table")

    script:
    """
    gatk CalculateContamination \
        --input ${pileup} \
        --output ${sample}.contamination.table
    """
}

process filterMutectCalls {
    input:
    tuple val(sample),
        path(reference),
        path(vcf),
        val(interval_id),
        path(intervals),
        path(contamination),
        path(orientation)

    output:
    tuple val(sample),
        path("${interval_id}.${sample}.filtered.vcf.gz"),
        path("${interval_id}.${sample}.filtered.vcf.gz.tbi")

    script:
    """
    gatk FilterMutectCalls \
        --reference ${reference[0]} \
        --variant ${vcf[0]} \
        --intervals ${intervals} \
        --contamination-table ${contamination} \
        --orientation-bias-artifact-priors ${orientation} \
        --unique-alt-read-count 3 \
        --output ${interval_id}.${sample}.filtered.vcf.gz
    """
}

process concatFilteredCalls {
    input:
    tuple val(sample),
        path(reference),
        path(vcfs),
        path(tbis)

    output:
    path(reference), emit: ref
    path("*.concatenated.vcf.gz"), emit: vcf
    path("*.concatenated.vcf.gz.tbi"), emit: tbi

    publishDir "${params.outdir}/Filtered/Samples", mode: 'symlink', pattern: '*.concatenated.vcf.gz*'

    script:
    """
    bcftools concat *.vcf.gz \
        | bcftools sort \
        | normalise_mutect2_vcf.py - \
        | bcftools norm -f ${reference[0]} \
            -Oz -o "${sample}.concatenated.vcf.gz"
    bcftools index -t "${sample}.concatenated.vcf.gz"
    """
}

process concatSecondHaplotypeCallerCalls {
    input:
    tuple val(sample),
        path(reference),
        path(vcfs),
        path(tbis)

    output:
    path(reference), emit: ref
    path("*.concatenated.vcf.gz"), emit: vcf
    path("*.concatenated.vcf.gz.tbi"), emit: tbi

    publishDir "${params.outdir}/SecondHaplotypeCallerCalls/Samples", mode: 'symlink', pattern: '*.concatenated.vcf.gz*'

    script:
    """
    bcftools concat *.vcf.gz \
        | bcftools sort \
        | bcftools norm -f ${reference[0]} -m -both \
            -Oz -o "${sample}.concatenated.vcf.gz"
    bcftools index -t "${sample}.concatenated.vcf.gz"
    """
}

process mergeFilteredCalls {
    input:
    path(reference)
    path(vcfs)
    path(tbis)

    output:
    path("Final.vcf.gz*")

    publishDir "${params.outdir}/Filtered", mode: 'symlink', pattern: 'Final.vcf.gz*'
    
    script:
    """
    bcftools merge \
        --threads ${task.cpus} \
        ${vcfs.sort { it.getName() }.join(' ')} \
        | bcftools norm -f ${reference[0]} -m + \
            -Oz -o Final.vcf.gz \
            -W=tbi
    """
}

process learnReadOrientationModel {
    input:
    tuple val(sample), path(f1r2)

    output:
    tuple val(sample), path("*.model.tar.gz")

    script:
    """
    gatk LearnReadOrientationModel \
        --input ${f1r2.join(' --input ')} \
        --output ${sample}.model.tar.gz
    """
}

workflow {
    // Load and process the reference
    ref_ch = Channel.fromPath(params.reference)
    fai_ch = indexReference(ref_ch)
    dict_ch = makeReferenceDict(ref_ch)
    ref_files = ref_ch.merge(fai_ch).merge(dict_ch)

    // Chop into intervals for scattering-gathering
    ivls = splitIntervals(ref_files, params.numIntervals).flatten()
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

    // Construct a channel for every combination of sample and interval,
    // separately for normals and tumours
    normal_interval_ch = normals.combine(ivls)
    tumour_interval_ch = tumours.combine(ivls)

    // Run Mutect2 on normals
    normals_for_calling_ch = ref_files.combine(normal_interval_ch)
        .map { fa, fai, dict, sample, normal_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, normal_bam, interval) }
    normals_called_by_mutect_ch = runMutectOnNormal(normals_for_calling_ch)
    
    // Collect the PON into a genomeDB
    ponImport_ch = normals_called_by_mutect_ch.groupTuple().combine(ivls, by: 0)
    pon_db = genomicsDBImport_PON(ponImport_ch)

    // Create the panel
    pon_ch = ref_files.combine(ivls.combine(pon_db, by: 0))
        .map { fa, fai, dict, interval_label, interval, genomeDB ->
            tuple(interval_label, [fa, fai, dict], interval, genomeDB) }
    panel_of_normals_vcfs = createPanelOfNormals(pon_ch)

    // Merge the Panel VCFs
    vcfs = panel_of_normals_vcfs.map { it -> it[0] }.collect()
    indices = panel_of_normals_vcfs.map { it -> it[1] }.collect()
    panel_of_normals = mergePonVCFs(vcfs, indices)

    // Run Mutect2 on tumours
    tumours_for_calling_ch = ref_files.combine(tumour_interval_ch)
        .map { fa, fai, dict, sample, tumour_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, tumour_bam, interval) }
    tumours_first_called_by_mutect_ch = runMutectOnTumour(tumours_for_calling_ch) 

    // Create genomicsDB from the first-round tumour mutect2 calls
    somatic_import_ch = tumours_first_called_by_mutect_ch.groupTuple().combine(ivls, by: 0)
    somatic_db = genomicsDBImport_Somatic(somatic_import_ch)

    // Use bcftools to merge somatic candidates
    // to_normalize_ch = ref_files.combine(tumours_first_called_by_mutect_ch).map {
    //     fa, fai, dict, interval_label, vcf, tbi, stats ->
    //     tuple(interval_label, [fa, fai, dict], vcf, tbi, stats)
    // }
    // somatic_norm_vcfs = bcftoolsNormalizeSomaticCandidates(to_normalize_ch)
    // somatic_merged_vcfs = bcftoolsMergeSomaticCandidatesByInterval(somatic_norm_vcfs.groupTuple())
    // vcfs = somatic_merged_vcfs.map { it -> it[0] }.collect()
    // indices = somatic_merged_vcfs.map { it -> it[1] }.collect()
    // bcftoolsConcatSomaticCandidates(vcfs, indices)

    // Extract somatic candidates from DB
    somatic_ch = ref_files.combine(ivls.combine(somatic_db, by: 0))
        .map { fa, fai, dict, interval_label, interval, genomeDB ->
            tuple(interval_label, [fa, fai, dict], interval, genomeDB) }
    somatic_candidates_by_interval = normalizeSomaticCandidates(extractSomaticCandidates(somatic_ch))
    candidate_vcfs = somatic_candidates_by_interval.map { it -> it[0] }.collect()
    candidate_indices = somatic_candidates_by_interval.map { it -> it[1] }.collect()
    somatic_candidates = mergeSomaticCandidates(candidate_vcfs, candidate_indices)                          

    // Run HaplotypeCaller on normals
    gvcfs_ch = runHaplotypeCallerOnNormal(normals_for_calling_ch)

    // Collect the gvcfs into a genomeDB
    dbImport_ch = gvcfs_ch.groupTuple().combine(ivls, by: 0)
    db = genomicsDBImport(dbImport_ch)

    // Genotype
    to_genotype_ch = ref_files.combine(ivls.combine(db, by: 0))
        .map { fa, fai, dict, interval_label, interval, genomeDB ->
            tuple(interval_label, [fa, fai, dict], interval, genomeDB) }
    genotyped_vcfs = genotypeGVCFs(to_genotype_ch)

    // Merge the genotyped VCFs
    vcfs = genotyped_vcfs.map { it -> it[0] }.collect()
    indices = genotyped_vcfs.map { it -> it[1] }.collect()
    germline_resource = mergeGenotypedVCFs(vcfs, indices)

    // Finalize the Somatic candidates by adding variants from the PON and GR
    collated_support_vcfs_ch = germline_resource
        .concat(panel_of_normals)
        .concat(somatic_candidates)
        .collate(3)
    to_finalize_ch = ref_files
        .combine(collated_support_vcfs_ch)
        .map { fa, fai, stats, gr, pon, sc ->
            tuple([fa, fai, stats], gr, pon, sc ) }
    
    candidates = finalizeSomaticCandidates(to_finalize_ch)

    // Call somatic variants (round 2 - final calling step before filtering)
    somatic_calling_ch = tumours_for_calling_ch
        .combine(candidates)
    unfiltered_calls_ch = callSomaticVariants(somatic_calling_ch)

    // Run on normals, too
    somatic_normals_ch = normals_for_calling_ch
        .combine(candidates)
    unfiltered_normals_ch = callSomaticVariantsOnNormal(somatic_normals_ch)

    // Rerun haplotype caller on normals at candidate sites
    rehaplotyped_normals = recallGermlineVariants(somatic_normals_ch)

    // Build a strand bias model
    orientation_model_ch = unfiltered_calls_ch.f1r2s
        .concat(unfiltered_normals_ch.f1r2s)
        .groupTuple(size: params.numIntervals)
        .map { key, files -> tuple(key, files.sort { it.name }) }
    orientation_model = learnReadOrientationModel(orientation_model_ch)
    
    // Build a contamination model
    to_pileup_ch = ref_files.combine(tumours.concat(normals)).combine(germline_resource)
        .map { fa, fai, dict, sample, bam, resource_vcf, resource_index ->
            tuple(sample, [fa, fai, dict], bam, [resource_vcf, resource_index]) }
    pileup = getPileupSummaries(to_pileup_ch)
    contamination_model = calculateContamination(pileup)

    // Filter calls
    to_filter_ch = ref_files.combine(unfiltered_calls_ch.vcfs
                                     .concat(unfiltered_normals_ch.vcfs))
            .map { fa, fai, dict, sample, vcf, tbi, stats, interval_id, intervals ->
                tuple(sample, [fa, fai, dict], [vcf, tbi, stats], interval_id, intervals) }
        .combine(contamination_model, by: 0)
        .combine(orientation_model, by: 0)

    filtered = filterMutectCalls(to_filter_ch)

    // Merge final variants
    grouped = filtered.groupTuple(size: params.numIntervals)
        .map { label, vcfs, indices -> tuple(label, vcfs.sort { it.name }, indices.sort { it.name }) }
    grouped_with_ref = ref_files.combine(grouped)
        .map { fa, fai, dict, sample, vcfs, tbis ->
            tuple(sample, [fa, fai, dict], vcfs, tbis) }
    concat = concatFilteredCalls(grouped_with_ref)
    mergeFilteredCalls(concat.ref, concat.vcf.collect(), concat.tbi.collect())

    // And the normals
    grouped_rehaplotyped = rehaplotyped_normals
        .map { it -> tuple(it[0], it[1], it[2]) } // sample, vcf, tbi
        .groupTuple(size: params.numIntervals)
        .map { label, vcfs, indices -> tuple(label, vcfs.sort { it.name }, indices.sort { it.name }) }
    grouped_rehaplotyped_with_ref = ref_files.combine(grouped_rehaplotyped)
        .map { fa, fai, dict, sample, vcfs, tbis ->
            tuple(sample, [fa, fai, dict], vcfs, tbis) }

    concat_rehaplotyped = concatSecondHaplotypeCallerCalls(grouped_rehaplotyped_with_ref)
}
