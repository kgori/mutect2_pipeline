process callSomaticVariants {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'

    input:
    tuple val(interval_id),
        path(reference),
        val(sample),
        path(bam),
        path(intervals),
        path(germline_resource),
        path(panel_of_normals),
        path(candidates)

    output:
    tuple val(sample),
        path("${interval_id}.${sample}.unfiltered.vcf.gz"),
        path("${interval_id}.${sample}.unfiltered.vcf.gz.tbi"),
        path("${interval_id}.${sample}.unfiltered.vcf.gz.stats"),
        val(interval_id),
        path(intervals),
        emit: vcfs
    tuple val(sample), path("*.f1r2.tar.gz*"), emit: f1r2s

    publishDir "${params.outdir}/SecondTumourCalls", mode: 'symlink', pattern: '*gz*'
    script:
    """
    # Restrict candidate set to the current interval
    gatk SelectVariants \
        --reference ${reference[0]} \
        --variant ${candidates[0]} \
        --intervals ${intervals} \
        --output candidates.subset.vcf.gz

    # Call somatic variants with a tight focus on the candidate positions
    gatk Mutect2 \
        --reference ${reference[0]} \
        --input ${bam[0]} \
        --germline-resource ${germline_resource[0]} \
        --panel-of-normals ${panel_of_normals[0]} \
        --alleles candidates.subset.vcf.gz \
        --force-call-filtered-alleles \
        --f1r2-tar-gz ${interval_id}.${sample}.f1r2.tar.gz \
        --output ${interval_id}.${sample}.unfiltered.vcf.gz \
        --intervals candidates.subset.vcf.gz \
        --native-pair-hmm-threads ${task.cpus} \
        --max-mnp-distance 0 \
        --assembly-region-padding 300
    """
}

process recallGermlineVariants {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'

    input:
    tuple val(interval_id),
        path(reference),
        val(sample),
        path(bam),
        path(intervals),
        path(germline_resource),
        path(panel_of_normals),
        path(candidates)

    output:
    tuple val(sample),
        path("${interval_id}.${sample}.haplotypeCaller.vcf.gz"),
        path("${interval_id}.${sample}.haplotypeCaller.vcf.gz.tbi"),
        val(interval_id),
        path(intervals),
        emit: vcfs

    publishDir "${params.outdir}/SecondHaplotypeCallerCalls", mode: 'symlink', pattern: '*gz*'
    script:
    """
    # Restrict candidate set to the current interval
    gatk SelectVariants \
        --reference ${reference[0]} \
        --variant ${candidates[0]} \
        --intervals ${intervals} \
        --output candidates.subset.vcf.gz

    # Call somatic variants with a tight focus on the candidate positions
    gatk HaplotypeCaller \
        --reference ${reference[0]} \
        --input ${bam[0]} \
        --alleles candidates.subset.vcf.gz \
        --population-callset ${germline_resource[0]} \
        --force-call-filtered-alleles \
        --output ${interval_id}.${sample}.haplotypeCaller.vcf.gz \
        --intervals candidates.subset.vcf.gz \
        --native-pair-hmm-threads ${task.cpus} \
        --max-mnp-distance 0 \
        --assembly-region-padding 300
    """
}
