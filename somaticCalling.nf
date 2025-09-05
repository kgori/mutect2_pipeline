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

    publishDir "${params.outdir}/SecondTumourCalls", mode: 'symlink'
    script:
    """
    gatk Mutect2 \
        --reference ${reference[0]} \
        --input ${bam[0]} \
        --germline-resource ${germline_resource[0]} \
        --panel-of-normals ${panel_of_normals[0]} \
        --alleles ${candidates[0]} \
        --f1r2-tar-gz ${interval_id}.${sample}.f1r2.tar.gz \
        --output ${interval_id}.${sample}.unfiltered.vcf.gz \
        --intervals ${intervals} \
        --interval-padding 150 \
        --native-pair-hmm-threads ${task.cpus} \
        --assembly-region-padding 300
    """
}
