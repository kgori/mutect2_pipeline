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

    output:
    path("${sample}.*.normal.mutect2_panel_calls.vcf.gz")
    path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.tbi")
    path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.stats")

    script:
    intervalNumberMatch = interval.getName() =~ /^(\d+)/
    intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
    """
    gatk Mutect2 \
      --reference ${reference[0]} \
      --input ${normal_bam[0]} \
      --intervals ${interval} \
      --interval-padding 150 \
      --native-pair-hmm-threads ${task.cpus} \
      --max-mnp-distance 0 \
      --output ${sample}.${intervalNumber}.normal.mutect2_panel_calls.vcf.gz
    """
}

process runMutectOnTumour {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple path(reference), val(sample), path(tumour_bam), path(interval)

    output:
    path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz")
    path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.tbi")
    path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.stats")

    script:
    intervalNumberMatch = interval.getName() =~ /^(\d+)/
    intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
    """
    gatk Mutect2 \
      --reference ${reference[0]} \
      --input ${tumour_bam[0]} \
      --intervals ${interval} \
      --interval-padding 150 \
      --native-pair-hmm-threads ${task.cpus} \
      --max-mnp-distance 0 \
      --output ${sample}.${intervalNumber}.tumour.mutect2_candidate_discovery_calls.vcf.gz
    """
}

process runHaplotypeCallerOnNormal {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(normal_bam), path(interval)

    output:
    tuple val(interval_id),
        path("${sample}.${interval_id}.haplotypecaller.g.vcf.gz"),
        path("${sample}.${interval_id}.haplotypecaller.g.vcf.gz.tbi")

    script:
    """
    gatk HaplotypeCaller \
      --reference ${reference[0]} \
      --input ${normal_bam[0]} \
      --intervals ${interval} \
      --interval-padding 150 \
      --native-pair-hmm-threads ${task.cpus} \
      --emit-ref-confidence GVCF \
        --output ${sample}.${interval_id}.haplotypecaller.g.vcf.gz
    """
}

process genomicsDBImport {
    input:
    tuple val(interval_id), path(gvcfs), path(gvcf_index), path(interval)

    output:
    path("genomicsDB")

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path genomicsDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --intervals ${interval} \
      -V ${gvcfs.join(' -V ')}
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

    // Construct a channel for every combination of sample and interval,
    // separately for normals and tumours
    normal_interval_ch = normals.combine(ivls)
    tumour_interval_ch = tumours.combine(ivls)

    // Run Mutect2 on normals
    normals_for_calling_ch = ref_files.combine(normal_interval_ch)
        .map { fa, fai, dict, sample, normal_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, normal_bam, interval) }
    // runMutectOnNormal(normals_for_calling_ch)

    // Run Mutect2 on tumours
    tumours_for_candidate_discovery_ch = ref_files.combine(tumour_interval_ch)
        .map { fa, fai, dict, sample, tumour_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, tumour_bam, interval) }
    // runMutectOnTumour(tumours_for_candidate_discovery_ch) 

    // Run HaplotypeCaller on normals
    gvcfs_ch = runHaplotypeCallerOnNormal(normals_for_calling_ch)
    dbImport_ch = gvcfs_ch.groupTuple().combine(ivls, by: 0)
    // intervals_ch = ivls

    // Collect the gvcfs into a genomeDB
    db = genomicsDBImport(dbImport_ch)
}
