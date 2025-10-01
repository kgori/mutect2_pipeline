process splitIntervals {
    input:
    path(reference)
    val(numIntervals)

    output:
    path("Intervals/*.interval_list")

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    touch intervals.list
    for i in {4..7}; do echo "\${i}" >> intervals.list; done
      
    gatk SplitIntervals \
      --reference "${reference[0]}" \
      --scatter-count $numIntervals \
      --intervals intervals.list \
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
