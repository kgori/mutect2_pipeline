process runContaminationModelOnPaired {
    memory { 4.GB + 4.GB * (task.attempt - 1) }
    maxRetries 3
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(sample),
        path(reference),
        path(tumour_bam),
        path(normal_bam),
        path(germline_resource)

    output:
    tuple val(sample), path("${sample}.contamination.table")

    publishDir "${params.outdir}/ContaminationTables", mode: 'copy', pattern: '*table'

    script:
    """
    bcftools view -v snps -Oz -o snps.vcf.gz -W=tbi ${germline_resource[0]}
    
    gatk GetPileupSummaries \
        --reference ${reference[0]} \
        --input ${tumour_bam[0]} \
        --variant snps.vcf.gz \
        --intervals snps.vcf.gz \
        --minimum-population-allele-frequency 0.01 \
        --maximum-population-allele-frequency 0.5 \
        --output ${sample}.tumour.pileups
    
    gatk GetPileupSummaries \
        --reference ${reference[0]} \
        --input ${normal_bam[0]} \
        --variant snps.vcf.gz \
        --intervals snps.vcf.gz \
        --minimum-population-allele-frequency 0.01 \
        --maximum-population-allele-frequency 0.5 \
        --output ${sample}.normal.pileups

    gatk CalculateContamination \
        --input ${sample}.tumour.pileups \
        --matched-normal ${sample}.normal.pileups \
        --output ${sample}.contamination.table
    """
}

process runContaminationModelOnTumourOnly {
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
    tuple val(sample), path("${sample}.contamination.table")

    publishDir "${params.outdir}/ContaminationTables", mode: 'copy', pattern: '*table'

    script:
    """
    bcftools view -v snps -Oz -o snps.vcf.gz -W=tbi ${germline_resource[0]}
    
    gatk GetPileupSummaries \
        --reference ${reference[0]} \
        --input ${bam[0]} \
        --variant snps.vcf.gz \
        --intervals snps.vcf.gz \
        --minimum-population-allele-frequency 0.01 \
        --maximum-population-allele-frequency 0.5 \
        --output ${sample}.pileups

    gatk CalculateContamination \
        --input ${sample}.pileups \
        --output ${sample}.contamination.table
    """
}
