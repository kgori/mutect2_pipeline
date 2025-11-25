process concatMutectVcfParts {
    input:
    tuple val(sample), val(label), path(reference), path(vcfs), path(tbis)

    output:
    path("${sample}.${label}.combined.vcf.gz"), emit: vcf
    path("${sample}.${label}.combined.vcf.gz.tbi"), emit: tbi

    publishDir "${params.outdir}/CombinedCalls/${label}", mode: 'copy'

    script:
    """
    bcftools concat -a -D *.vcf.gz \
        | bcftools sort \
        | split_mutect2_multiallelics.py - \
        | bcftools norm -f ${reference[0]} -d exact \
             -Oz -o "${sample}.${label}.combined.vcf.gz"
    bcftools index -t "${sample}.${label}.combined.vcf.gz"
    """
}

process concatVcfParts {
    input:
    tuple val(sample), val(label), path(reference), path(vcfs), path(tbis)

    output:
    path("${sample}.${label}.combined.vcf.gz"), emit: vcf
    path("${sample}.${label}.combined.vcf.gz.tbi"), emit: tbi

    publishDir "${params.outdir}/CombinedCalls/${label}", mode: 'copy'

    script:
    """
    bcftools concat -a -D *.vcf.gz \
        | bcftools sort \
        | bcftools norm -f ${reference[0]} -m -any -d exact \
        -Oz -o "${sample}.${label}.combined.vcf.gz"
    bcftools index -t "${sample}.${label}.combined.vcf.gz"
    """
}
