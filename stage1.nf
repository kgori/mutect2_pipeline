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
include { runMutectOnNormal }                                  from "./initialVariantCalling.nf"
include { runMutectOnTumour }                                  from "./initialVariantCalling.nf"
include { runHaplotypeCallerOnNormal }                         from "./initialVariantCalling.nf"


workflow {
    /////////////////////////////////
    //         Preanalysis setup
    //

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


    ///////////////////////////////////////////////////////
    //         Stage 1: First round of calling
    //

    // Run Mutect2 on normals
    normals_for_calling_ch = ref_files.combine(normal_interval_ch)
        .map { fa, fai, dict, sample, normal_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, normal_bam, interval) }
    normals_called_by_mutect_ch = runMutectOnNormal(normals_for_calling_ch)

    // Run Mutect2 on tumours
    tumours_for_calling_ch = ref_files.combine(tumour_interval_ch)
        .map { fa, fai, dict, sample, tumour_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, tumour_bam, interval) }
    tumours_first_called_by_mutect_ch = runMutectOnTumour(tumours_for_calling_ch)

    // Run HaplotypeCaller on normals
    gvcfs_ch = runHaplotypeCallerOnNormal(normals_for_calling_ch)
}
