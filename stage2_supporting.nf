process genotypeGvcfIntervals {
    input:
    tuple val(interval_id), path(ref), path(gvcfs), path(gvcf_index), path(interval)

    output:
    tuple val(interval_id),
        path("${interval_id}.genotyped.vcf.gz"),
        path("${interval_id}.genotyped.vcf.gz.tbi")

    publishDir "${params.outdir}/GermlineResource", mode: 'copy'

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path resourceDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --merge-input-intervals true \
      --intervals ${interval} \
      -V ${gvcfs.join(' -V ')}

    gatk GenotypeGVCFs \
        --reference ${ref[0]} \
        --variant gendb://resourceDB \
        --intervals ${interval} \
        --output ${interval_id}.genotyped.vcf.gz \
        && rm -rf resourceDB
    """
}


process makePonIntervals {
    input:
    tuple val(interval_id), path(ref), path(gvcfs), path(gvcf_index), path(interval)

    output:
    tuple val(interval_id),
        path("${interval_id}.panel_of_normals.vcf.gz"),
        path("${interval_id}.panel_of_normals.vcf.gz.tbi")

    publishDir "${params.outdir}/PanelOfNormals", mode: 'copy'

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path resourceDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --merge-input-intervals true \
      --intervals ${interval} \
      -V ${gvcfs.join(' -V ')}

    gatk CreateSomaticPanelOfNormals \
        --reference ${ref[0]} \
        --variant gendb://resourceDB \
        --intervals ${interval} \
        --output ${interval_id}.panel_of_normals.vcf.gz \
        && rm -rf resourceDB
    """
}

process makeSomaticCandidatesIntervals {
    input:
    tuple val(interval_id), path(ref), path(gvcfs), path(gvcf_index), path(interval)

    output:
    tuple val(interval_id),
        path("${interval_id}.somatic_candidates.vcf.gz"),
        path("${interval_id}.somatic_candidates.vcf.gz.tbi")

    publishDir "${params.outdir}/SomaticCandidates", mode: 'copy'

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path resourceDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --merge-input-intervals true \
      --intervals ${interval} \
      -V ${gvcfs.join(' -V ')}

    gatk SelectVariants \
        --reference ${ref[0]} \
        --variant gendb://resourceDB \
        --intervals ${interval} \
        --output ${interval_id}.intermediate.vcf.gz \
        && rm -rf resourceDB

    bcftools norm -m -both -f ${ref[0]} ${interval_id}.intermediate.vcf.gz \
        | bcftools view -e 'ALT="*"' \
            -Oz -o ${interval_id}.somatic_candidates.vcf.gz
    bcftools index -t ${interval_id}.somatic_candidates.vcf.gz
    """
}
