nextflow.enable.dsl=2

params.reference    = "genome.fa"
params.normals      = "normals_folder"
params.tumours      = "tumours_folder"
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

process splitIntervals {
  input:
  path(reference)
  val(numIntervals)

  output:
  path("Intervals/*.interval_list")

  publishDir "results", mode: 'copy'

  script:
  """
  gatk SplitIntervals \
    --reference "${reference[0]}" \
    --scatter-count $numIntervals \
    --intervals 38 \
    --output Intervals
  """
}

process indexReference {
    input:
    path(reference)

    output:
    path("${reference.getName()}.fai")

    script:
    """
    samtools faidx ${reference}
    """
}

process makeReferenceDict {
    input:
    path(reference)

    output:
    path("*.dict")

    script:
    def refBase = reference.getName().replaceFirst(/\.[^.]+$/, '')
    """
    gatk CreateSequenceDictionary \
      -R ${reference} \
      -O ${refBase}.dict
    """
}

process runMutectOnNormal {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple path(reference), val(sample), path(normal_bam), path(interval)

    // output:
    // path("${normal_bam.getName()}.vcf.gz")

    script:
    """
    echo "${reference[0]} ${reference[1]} ${reference[2]}" > refs.txt
    echo "${sample} ${normal_bam[0]} ${normal_bam[1]} ${interval}" > inputs.txt
    gatk Mutect2 \
      --reference ${reference[0]} \
      --input ${normal_bam[0]} \
      --intervals ${interval} \
      --interval-padding 150 \
      --native-pair-hmm-threads ${task.cpus} \
      --max-mnp-distance 0 \
      --output ${sample}.normal_panel_calls.vcf.gz
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

    // Construct a channel for every combination of sample and interval,
    // separately for normals and tumours
    normal_interval_ch = normals.combine(ivls)
    tumour_interval_ch = tumours.combine(ivls)

    // Run Mutect2 on normals
    normals_for_pon_calling_ch = ref_files.combine(normal_interval_ch)
        .map { fa, fai, dict, sample, normal_bam, interval ->
            tuple([fa, fai, dict], sample, normal_bam, interval) }
    runMutectOnNormal(normals_for_pon_calling_ch)
}
