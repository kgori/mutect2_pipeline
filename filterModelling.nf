process getPileupSummaries {
    memory { 4.GB + 4.GB * (task.attempt - 1) }
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'
    
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
    memory { 4.GB + 4.GB * (task.attempt - 1) }
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'

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
    memory { 4.GB + 4.GB * (task.attempt - 1) }
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'
    
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
