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

process concatGvcfParts {
    input:
    tuple val(sample), val(label), path(reference), path(vcfs), path(tbis)

    output:
    path("${sample}.${label}.combined.g.vcf.gz"), emit: vcf
    path("${sample}.${label}.combined.g.vcf.gz.tbi"), emit: tbi

    publishDir "${params.outdir}/CombinedCalls/${label}", mode: 'copy'

    script:
    """
    bcftools concat -a -D *.g.vcf.gz \
        -Oz -o "${sample}.${label}.intermediate.g.vcf.gz"
    bcftools index -t "${sample}.${label}.intermediate.g.vcf.gz"

    gatk ReblockGVCFs \
        --reference ${reference[0]} \
        --input "${sample}.${label}.intermediate.g.vcf.gz" \
        --output "${sample}.${label}.combined.g.vcf.gz"

    rm "${sample}.${label}.intermediate.g.vcf.gz"*
    """
}

process concatFilteredCalls {
    memory { 10.GB + 4.GB * (task.attempt - 1) }
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample),
        path(reference),
        path(pon),
        path(vcfs),
        path(tbis)

    output:
    path(reference), emit: ref
    path("*.concatenated.vcf.gz"), emit: vcf
    path("*.concatenated.vcf.gz.tbi"), emit: tbi

    publishDir "${params.outdir}/Filtered/Samples", mode: 'copy', pattern: '*.concatenated.vcf.gz*'

    script:
    """
    bcftools concat *.filtered.vcf.gz \
        | bcftools sort \
        | split_mutect2_multiallelics.py - \
        | bcftools norm -f ${reference[0]} \
        | fix_panel_of_normals_annotation.py  \
            --reference-fasta ${reference[0]} \
            --pon-vcf ${pon[0]} \
            -o "${sample}.concatenated.vcf.gz" -
    bcftools index -t "${sample}.concatenated.vcf.gz"
    """
}

process concatSecondHaplotypeCallerCalls {
    memory { 10.GB + 4.GB * (task.attempt - 1) }
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample),
        path(reference),
        path(vcfs),
        path(tbis)

    output:
    path("*.concatenated.vcf.gz"), emit: vcf
    path("*.concatenated.vcf.gz.tbi"), emit: tbi

    publishDir "${params.outdir}/SecondHaplotypeCallerCalls/Samples", mode: 'copy', pattern: '*.concatenated.vcf.gz*'

    script:
    """
    bcftools concat *.vcf.gz \
        | bcftools sort \
        | bcftools norm -f ${reference[0]} -m -both \
            -Oz -o "${sample}.concatenated.vcf.gz"
    bcftools index -t "${sample}.concatenated.vcf.gz"
    """
}
