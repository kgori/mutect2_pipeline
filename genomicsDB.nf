process genomicsDBImport {
    input:
    tuple val(interval_id), path(gvcfs), path(gvcf_index), path(interval)

    output:
    tuple val(interval_id), path("resourceDB")

    publishDir "${params.outdir}/GenomicsDB/GermlineResource", mode: 'link'

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path resourceDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --intervals ${interval} \
      -V ${gvcfs.join(' -V ')}
    """
}

process genomicsDBImport_PON {
    input:
    tuple val(interval_id), path(vcfs), path(vcf_index), path(statsfiles), path(interval)

    output:
    tuple val(interval_id), path("ponDB")

    publishDir "${params.outdir}/GenomicsDB/PanelOfNormals", mode: 'link'

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path ponDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --intervals ${interval} \
      -V ${vcfs.join(' -V ')}
    """
}

process genomicsDBImport_Somatic {
    input:
    tuple val(interval_id), path(vcfs), path(vcf_index), path(statsfiles), path(interval)

    output:
    tuple val(interval_id), path("somaticDB")
    
    publishDir "${params.outdir}/GenomicsDB/SomaticCandidates", mode: 'link'

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path somaticDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --intervals ${interval} \
      -V ${vcfs.join(' -V ')}
    """
}
