process extractSomaticCandidates {
    input:
    tuple val(interval_id), path(reference), path(interval), path(genomeDB)

    output:
    tuple val(interval_id), path(reference), path("${interval_id}.somatic_candidates_raw.vcf.gz"),
        path("${interval_id}.somatic_candidates_raw.vcf.gz.tbi")

    publishDir "${params.outdir}/SomaticCandidates/GATK", mode: 'copy'

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

    publishDir "${params.outdir}/SomaticCandidates/GATK", mode: 'copy'
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

    publishDir "${params.outdir}/SomaticCandidates/GATK", mode: 'copy'

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

    publishDir "${params.outdir}/SomaticCandidates/Bcftools", mode: 'copy'

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

    publishDir "${params.outdir}/SomaticCandidates/Bcftools", mode: 'copy'

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

    publishDir "${params.outdir}/SomaticCandidates/BcfTools", mode: 'copy'

    script:
    """
    bcftools concat -a -D ${vcfs.sort { it.getName() }.join(' ')} \
        | bcftools sort -Oz \
        -o bcftoolsSomaticCandidates.vcf.gz
    bcftools index -t bcftoolsSomaticCandidates.vcf.gz
    """
}
    
process finalizeSomaticCandidates {
    input:
    tuple path(reference), path(germline_resource), path(panel_of_normals), path(somatic_candidates)

    output:
    tuple path(germline_resource), path(panel_of_normals), path("candidates.vcf.gz*")

    publishDir "${params.outdir}/SomaticCandidates/Final", mode: 'copy', pattern: '*.vcf.gz*'

    script:
    """
    bcftools norm -m -both -f ${reference[0]} ${germline_resource[0]} \
        | bcftools view -Oz -o gr.vcf.gz -S ^<(bcftools query -l ${germline_resource[0]})
    bcftools index -t gr.vcf.gz
    bcftools norm -m -both -f ${reference[0]} ${panel_of_normals[0]} \
        | bcftools view -Oz -o pon.vcf.gz -S ^<(bcftools query -l ${panel_of_normals[0]})
    bcftools index -t pon.vcf.gz
    bcftools norm -m -both -f ${reference[0]} ${somatic_candidates[0]} \
        | bcftools view -Oz -o sc.vcf.gz -S ^<(bcftools query -l ${somatic_candidates[0]})
    bcftools index -t sc.vcf.gz

    bcftools concat -a -d exact gr.vcf.gz pon.vcf.gz sc.vcf.gz \
        | bcftools view -e 'TYPE="indel" && strlen(REF) - strlen(ALT) > 150' \
        | bcftools view -e 'ALT="*"' \
        | bcftools sort \
        | bcftools norm -f ${reference[0]} \
        | bcftools norm -d exact \
        | bcftools annotate -x INFO,QUAL -Oz -o candidates.vcf.gz
    bcftools index -t candidates.vcf.gz

    rm gr.vcf.gz* pon.vcf.gz* sc.vcf.gz*
    """
}
