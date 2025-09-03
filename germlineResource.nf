process genotypeGVCFs {
    input:
    tuple val(interval_id), path(reference), path(interval), path(genomeDB)

    output:
    tuple path("${interval_id}.genotyped.vcf.gz"),
        path("${interval_id}.genotyped.vcf.gz.tbi")

    script:
    """
    gatk GenotypeGVCFs \
        --reference ${reference[0]} \
        --variant gendb://${genomeDB} \
        --intervals ${interval} \
        --output ${interval_id}.genotyped.vcf.gz
    """
}

process mergeGenotypedVCFs {
    executor 'local'
    input:
    path(genotyped_vcfs)
    path(vcf_indices)

    output:
    path("germline_resource.vcf.gz*")

    publishDir "${params.outdir}/GermlineResource", mode: 'copy'
    
    script:
    """
    bcftools concat -a -D ${genotyped_vcfs.join(' ')} \
        | bcftools sort -Oz -o germline_resource.vcf.gz
    bcftools index -t germline_resource.vcf.gz
    """
}
