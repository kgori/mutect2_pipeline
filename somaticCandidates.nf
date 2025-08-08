process extractSomaticCandidates {
    input:
    tuple val(interval_id), path(reference), path(interval), path(genomeDB)

    output:
    tuple val(interval_id), path(reference), path("${interval_id}.somatic_candidates_raw.vcf.gz"),
        path("${interval_id}.somatic_candidates_raw.vcf.gz.tbi")

    publishDir "results/SomaticCandidates/GATK", mode: 'copy'

    script:
    """
    gatk SelectVariants \
        --reference ${reference[0]} \
        --variant gendb://${genomeDB} \
        --intervals ${interval} \
        --output ${interval_id}.somatic_candidates_raw.vcf.gz
    """
}

process normalizeSomaticCandidates {
    input:
    tuple val(interval_id), path(reference), path(vcf), path(index)

    output:
    tuple path("${interval_id}.somatic_candidates.vcf.gz"),
        path("${interval_id}.somatic_candidates.vcf.gz.tbi")

    publishDir "results/SomaticCandidates/GATK", mode: 'copy'
    script:
    """
    bcftools norm -m -both -f ${reference[0]} ${vcf} \
        | bcftools view -e 'ALT="*"' \
            -Oz -o ${interval_id}.somatic_candidates.vcf.gz
    bcftools index -t ${interval_id}.somatic_candidates.vcf.gz
    """
}

process mergeSomaticCandidates {
    executor 'local'
    input:
    path(candidate_vcfs)
    path(vcf_indices)

    output:
    path("somatic_candidates.vcf.gz*")

    publishDir "results/SomaticCandidates/GATK", mode: 'copy'

    script:
    """
    bcftools concat -a -D ${candidate_vcfs.join(' ')} \
        | bcftools sort -Oz -o somatic_candidates.vcf.gz
    bcftools index -t somatic_candidates.vcf.gz
    """
}

process bcftoolsNormalizeSomaticCandidates {
    input:
    tuple val(interval_id), path(reference), path(vcf), path(index), path(stats)

    output:
    tuple val(interval_id),
        path("${vcf.getBaseName(2)}.norm.vcf.gz"),
        path("${vcf.getBaseName(2)}.norm.vcf.gz.tbi"),
        path(stats)

    publishDir "results/SomaticCandidates/Bcftools", mode: 'copy'

    script:
    """
    bcftools norm -m -both -f ${reference[0]} ${vcf} \
        | bcftools sort \
            -Oz -o ${vcf.getBaseName(2)}.norm.vcf.gz
    bcftools index -t ${vcf.getBaseName(2)}.norm.vcf.gz
    """
}

process bcftoolsMergeSomaticCandidatesByInterval {
    input:
    tuple val(interval_id), path(vcfs), path(index), path(stats)

    output:
    tuple path("${interval_id}.somatic_candidates.vcf.gz"),
        path("${interval_id}.somatic_candidates.vcf.gz.tbi")

    publishDir "results/SomaticCandidates/Bcftools", mode: 'copy'

    script:
    """
    bcftools merge -m none ${vcfs.sort { it.getName() }.join(' ')} \
        | bcftools sort -Oz -o ${interval_id}.somatic_candidates.vcf.gz
    bcftools index -t ${interval_id}.somatic_candidates.vcf.gz
    """
}

process bcftoolsConcatSomaticCandidates {
    input:
    path(vcfs)
    path(indexes)

    output:
    path("bcftoolsSomaticCandidates.vcf.gz*")

    publishDir "results/SomaticCandidates/BcfTools", mode: 'copy'

    script:
    """
    bcftools concat -a -D ${vcfs.sort { it.getName() }.join(' ')} \
        | bcftools sort -Oz \
        -o bcftoolsSomaticCandidates.vcf.gz
    bcftools index -t bcftoolsSomaticCandidates.vcf.gz
    """
}
    
