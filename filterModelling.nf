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
        path(stats),
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
        --stats ${stats} \
        --intervals ${intervals} \
        --contamination-table ${contamination} \
        --orientation-bias-artifact-priors ${orientation} \
        --unique-alt-read-count 3 \
        --output ${interval_id}.${sample}.filtered.vcf.gz
    """
}

process filterMutectCallsOnNormals {
    memory { 4.GB + 4.GB * (task.attempt - 1) }
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(sample),
        path(reference),
        path(vcf),
        path(stats),
        val(interval_id),
        path(intervals),
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
        --stats ${stats} \
        --intervals ${intervals} \
        --orientation-bias-artifact-priors ${orientation} \
        --unique-alt-read-count 3 \
        --output ${interval_id}.${sample}.filtered.vcf.gz
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
