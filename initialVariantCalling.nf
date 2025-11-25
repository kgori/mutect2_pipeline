process runMutectOnNormal {
    cpus 8
    memory { 8.GB + 4.GB * (task.attempt - 1) }
    errorStrategy 'retry'
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(normal_bam), path(interval)

    output:
    tuple val(sample), path("${sample}.*.normal.mutect2_panel_calls.vcf.gz"),
        path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.tbi"),
        path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.stats")

    publishDir "${params.outdir}/InitialCalls/Mutect/Normal", mode: 'copy'
    
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
    cpus 8
    memory { 8.GB + 4.GB * (task.attempt - 1) }
    errorStrategy 'retry'
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(tumour_bam), path(interval)

    output:
    tuple val(sample),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz"),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.tbi"),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.stats")

    publishDir "${params.outdir}/InitialCalls/Mutect/Tumour", mode: 'copy'
    
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
    memory { 8.GB + 4.GB * (task.attempt - 1) }
    errorStrategy 'retry'
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(normal_bam), path(interval)

    output:
    tuple val(sample),
        path("${sample}.${interval_id}.haplotypecaller.g.vcf.gz"),
        path("${sample}.${interval_id}.haplotypecaller.g.vcf.gz.tbi")

    publishDir "${params.outdir}/InitialCalls/HaplotypeCaller", mode: 'copy'

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

process runPlatypus {
    cpus 4
    memory { 16.GB + 4.GB * (task.attempt - 1) }
    errorStrategy 'retry'
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'

    input:
    tuple val(interval_id), path(reference), val(sample), path(bam), path(interval)

    output:
    tuple val(sample),
        path("${sample}.${interval_id}.platypus.vcf.gz"),
        path("${sample}.${interval_id}.platypus.vcf.gz.tbi")

    publishDir "${params.outdir}/InitialCalls/Platypus", mode: 'copy'

    script:
    """
    awk 'BEGIN { OFS="\\t" } { if ( \$0 ~ /^@/ ) { next; } else { print \$1,\$2-1,\$3 } }' ${interval} > interval.bed
    platypus callVariants \
        --regions=interval.bed \
        --bamFiles=${bam[0]} \
        --refFile=${reference[0]} \
        --output=${sample}.${interval_id}.platypus.vcf \
        --nCPU=${task.cpus} \
        --logFile=${sample}.${interval_id}.platypus.log
    bgzip ${sample}.${interval_id}.platypus.vcf
    tabix ${sample}.${interval_id}.platypus.vcf.gz
    """
}

process runFreeBayes {
    cpus 1
    memory { 16.GB + 4.GB * (task.attempt - 1) }
    errorStrategy 'retry'
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'

    input:
    tuple val(interval_id), path(reference), val(sample), path(bam), path(interval)

    output:
    tuple val(sample),
        path("${sample}.${interval_id}.freebayes.vcf.gz"),
        path("${sample}.${interval_id}.freebayes.vcf.gz.tbi")

    publishDir "${params.outdir}/InitialCalls/FreeBayes", mode: 'copy'

    script:
    """
    awk 'BEGIN { OFS="\\t" } { if ( \$0 ~ /^@/ ) { next; } else { print \$1,\$2-1,\$3 } }' ${interval} > interval.bed
    freebayes \
        -f ${params.reference} \
        -t interval.bed \
        -b ${bam[0]} | bgzip -c  > ${sample}.${interval_id}.freebayes.vcf.gz
    bcftools index -t ${sample}.${interval_id}.freebayes.vcf.gz
    """
}
