process collectReadCounts {
    input:
    tuple path(reference), val(sample), path(bam), path(bins)

    output:
    path("${sample}.read_counts.tsv.gz")

    publishDir "${params.outdir}/CopyNumber/ReadCounts", mode: 'copy'
    
    script:
    """
    gatk CollectReadCounts \
        --intervals ${bins} \
        --reference ${reference[0]} \
        --input ${bam[0]} \
        --format TSV \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output ${sample}.read_counts.tsv

    gzip ${sample}.read_counts.tsv
    """
}
