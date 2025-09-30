process createPanelOfNormals {
    input:
    tuple val(interval_id), path(reference), path(interval), path(genomeDB)

    output:
    tuple path("${interval_id}.panel_of_normals.vcf.gz"),
        path("${interval_id}.panel_of_normals.vcf.gz.tbi")

    script:
    """
    gatk CreateSomaticPanelOfNormals \
        --reference ${reference[0]} \
        --variant gendb://${genomeDB} \
        --intervals ${interval} \
        --output ${interval_id}.panel_of_normals.vcf.gz
    """
}

process mergePonVCFs {
    executor 'local'
    input:
    path(pon_vcfs)
    path(vcf_indices)

    output:
    path("panel_of_normals.vcf.gz*")

    publishDir "${params.outdir}/PON", mode: 'copy'

    script:
    """
    bcftools concat -a -D ${pon_vcfs.join(' ')} \
        | bcftools sort -Oz -o panel_of_normals.vcf.gz
    bcftools index -t panel_of_normals.vcf.gz
    """
}
