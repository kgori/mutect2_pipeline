process runMutectOnNormal {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(normal_bam), path(interval)

    output:
    tuple val(interval_id), path("${sample}.*.normal.mutect2_panel_calls.vcf.gz"),
        path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.tbi"),
        path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.stats")

    publishDir "${params.outdir}/InitialNormalMutectCalls", mode: 'copy'
    
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
    tuple val(interval_id), path(reference), val(sample), path(tumour_bam), path(interval)

    output:
    tuple val(interval_id),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz"),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.tbi"),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.stats")

    publishDir "${params.outdir}/InitialTumourMutectCalls", mode: 'copy'
    
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

    publishDir "${params.outdir}/InitialNormalHaplotypeCallerCalls", mode: 'copy'

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
