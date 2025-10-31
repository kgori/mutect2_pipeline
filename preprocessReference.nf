process splitIntervals {
    input:
    path(reference)
    val(numIntervals)
    path(intervals)

    output:
    path("Intervals/*.interval_list")

    publishDir "${params.outdir}", mode: 'copy'

    script:
    def intervalsArg = intervals.name == "NO_FILE" ? "" : "--intervals ${intervals}"
    """
    gatk SplitIntervals \
        --reference "${reference[0]}" \
        --scatter-count $numIntervals \
        ${intervalsArg} \
        --output Intervals
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
