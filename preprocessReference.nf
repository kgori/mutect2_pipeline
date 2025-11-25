process splitIntervals {
    input:
    path(reference)
    val(fineNum)
    val(coarseNum)
    path(intervals)

    output:
    path("FineIntervals/*.interval_list"), emit: fineIntervals
    path("CoarseIntervals/*.interval_list"), emit: coarseIntervals

    publishDir "${params.outdir}/Intervals", mode: 'copy'

    script:
    def intervalsArg = intervals.name == "NO_FILE" ? "" : "--intervals ${intervals}"
    """
    gatk SplitIntervals \
        --reference "${reference[0]}" \
        --scatter-count $fineNum \
        ${intervalsArg} \
        --output FineIntervals

    gatk SplitIntervals \
        --reference "${reference[0]}" \
        --scatter-count $coarseNum \
        ${intervalsArg} \
        --output CoarseIntervals
    """
}

process indexReference {
    input:
    path(reference)

    output:
    path("${reference.getName()}.fai")

    script:
    """
    samtools faidx ${reference}
    """
}

process makeReferenceDict {
    input:
    path(reference)

    output:
    path("*.dict")

    script:
    def refBase = reference.getName().replaceFirst(/\.[^.]+$/, '')
    """
    gatk CreateSequenceDictionary \
      -R ${reference} \
      -O ${refBase}.dict
    """
}
