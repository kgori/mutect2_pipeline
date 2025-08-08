nextflow.enable.dsl=2

params.reference    = "genome.fa"
params.normals      = "normals_folder"
params.tumours      = "tumours_folder"
params.numIntervals = 8


def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}

include { makeBamToSampleNameMap }                   from "./preprocessSamples.nf"
include { indexReference }                           from "./preprocessReference.nf"
include { splitIntervals }                           from "./preprocessReference.nf"
include { makeReferenceDict }                        from "./preprocessReference.nf"
include { runMutectOnNormal }                        from "./variantCalling.nf"
include { runMutectOnTumour }                        from "./variantCalling.nf"
include { runHaplotypeCallerOnNormal }               from "./variantCalling.nf"
include { genomicsDBImport }                         from "./genomicsDB.nf"
include { genomicsDBImport_PON }                     from "./genomicsDB.nf"
include { genomicsDBImport_Somatic }                 from "./genomicsDB.nf"
include { genotypeGVCFs }                            from "./germlineResource.nf"
include { mergeGenotypedVCFs }                       from "./germlineResource.nf"
include { createPanelOfNormals }                     from "./panelOfNormals.nf"
include { mergePonVCFs }                             from "./panelOfNormals.nf"
include { extractSomaticCandidates }                 from "./somaticCandidates"
include { normalizeSomaticCandidates }               from "./somaticCandidates"
include { mergeSomaticCandidates }                   from "./somaticCandidates"
include { bcftoolsNormalizeSomaticCandidates }       from "./somaticCandidates"
include { bcftoolsMergeSomaticCandidatesByInterval } from "./somaticCandidates"
include { bcftoolsConcatSomaticCandidates }          from "./somaticCandidates"

process finalizeSomaticCandidates {
    input:
    tuple path(reference), path(germline_resource), path(panel_of_normals), path(somatic_candidates)

    output:
    path("candidates.vcf.gz*")

    publishDir "results/SomaticCandidates/Final", mode: 'copy'
    
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

    bcftools concat -a -d all gr.vcf.gz pon.vcf.gz sc.vcf.gz \
        | bcftools view -e 'TYPE="indel" && strlen(REF) - strlen(ALT) > 150' \
        | bcftools view -e 'ALT="*"' \
        | bcftools sort \
        | bcftools norm -d all \
        | bcftools annotate -x INFO,QUAL -Oz -o candidates.vcf.gz
    bcftools index -t candidates.vcf.gz

    rm gr.vcf.gz* pon.vcf.gz* sc.vcf.gz*
    """
}

process callSomaticVariants {
    input:
    tuple val(intervals_id),
        path(bam),
        path(reference),
        path(germline_resource),
        path(panel_of_normals),
        path(candidates),
        path(intervals)

    output:
    path("*.unfiltered.vcf.gz*")

    script:
    """
    gatk Mutect2 \
        --reference ${reference[0]} \
        --input ${bam} \
        --germline-resource ${germline_resource[0]} \
        --panel-of-normals ${panel_of_normals[0]} \
        --alleles ${candidates[0]} \
        --f1r2-tar-gz f1r2.tar.gz \
        --output mutect2.unfiltered.vcf.gz \
        --intervals ${intervals} \
        --interval-padding 150 \
        --assembly-region-padding 300
    """
}

workflow {
    // Load and process the reference
    ref_ch = Channel.fromPath(params.reference)
    fai_ch = indexReference(ref_ch)
    dict_ch = makeReferenceDict(ref_ch)
    ref_files = ref_ch.merge(fai_ch).merge(dict_ch)

    // Chop into intervals for scattering-gathering
    ivls = splitIntervals(ref_files, params.numIntervals).flatten()
        .map { interval ->
            def intervalNumberMatch = interval.getName() =~ /^(\d+)/
                def intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
            return tuple(intervalNumber, interval) }

    // Load normal bams
    bam_bai_inputs   = Channel.fromFilePairs("${params.normals}/*.bam{,.bai}")
    bam_csi_inputs   = Channel.fromFilePairs("${params.normals}/*.bam{,.csi}")
    bam_inputs = remove_duplicate_filepair_keys(bam_csi_inputs, bam_bai_inputs)
    cram_inputs = Channel.fromFilePairs("${params.normals}/*.cram{,.crai}")
    normals = bam_inputs.concat(cram_inputs)

    // Load tumour bams
    bam_bai_inputs   = Channel.fromFilePairs("${params.tumours}/*.bam{,.bai}")
    bam_csi_inputs   = Channel.fromFilePairs("${params.tumours}/*.bam{,.csi}")
    bam_inputs = remove_duplicate_filepair_keys(bam_csi_inputs, bam_bai_inputs)
    cram_inputs = Channel.fromFilePairs("${params.tumours}/*.cram{,.crai}")
    tumours = bam_inputs.concat(cram_inputs)

    // Extract sample names from BAMs
    normal_name_map_ch = normals
        .map { label, bam -> bam[0] }
    tumour_name_map_ch = tumours
        .map { label, bam -> bam[0] }
    makeBamToSampleNameMap(normal_name_map_ch.concat(tumour_name_map_ch).collect())

    // Construct a channel for every combination of sample and interval,
    // separately for normals and tumours
    normal_interval_ch = normals.combine(ivls)
    tumour_interval_ch = tumours.combine(ivls)

    // Run Mutect2 on normals
    normals_for_calling_ch = ref_files.combine(normal_interval_ch)
        .map { fa, fai, dict, sample, normal_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, normal_bam, interval) }
    normals_called_by_mutect_ch = runMutectOnNormal(normals_for_calling_ch)
    
    // Collect the PON into a genomeDB
    ponImport_ch = normals_called_by_mutect_ch.groupTuple().combine(ivls, by: 0)
    pon_db = genomicsDBImport_PON(ponImport_ch)

    // Create the panel
    pon_ch = ref_files.combine(ivls.combine(pon_db, by: 0))
        .map { fa, fai, dict, interval_label, interval, genomeDB ->
            tuple(interval_label, [fa, fai, dict], interval, genomeDB) }
    panel_of_normals_vcfs = createPanelOfNormals(pon_ch)

    // Merge the Panel VCFs
    vcfs = panel_of_normals_vcfs.map { it -> it[0] }.collect()
    indices = panel_of_normals_vcfs.map { it -> it[1] }.collect()
    panel_of_normals = mergePonVCFs(vcfs, indices)

    // Run Mutect2 on tumours
    tumours_for_candidate_discovery_ch = ref_files.combine(tumour_interval_ch)
        .map { fa, fai, dict, sample, tumour_bam, interval_label, interval ->
            tuple(interval_label, [fa, fai, dict], sample, tumour_bam, interval) }
    tumours_first_called_by_mutect_ch = runMutectOnTumour(tumours_for_candidate_discovery_ch) 

    // Create genomicsDB from the first-round tumour mutect2 calls
    somatic_import_ch = tumours_first_called_by_mutect_ch.groupTuple().combine(ivls, by: 0)
    somatic_db = genomicsDBImport_Somatic(somatic_import_ch)

    // Use bcftools to merge somatic candidates
    to_normalize_ch = ref_files.combine(tumours_first_called_by_mutect_ch).map {
        fa, fai, dict, interval_label, vcf, tbi, stats ->
        tuple(interval_label, [fa, fai, dict], vcf, tbi, stats)
    }
    somatic_norm_vcfs = bcftoolsNormalizeSomaticCandidates(to_normalize_ch)
    somatic_merged_vcfs = bcftoolsMergeSomaticCandidatesByInterval(somatic_norm_vcfs.groupTuple())
    vcfs = somatic_merged_vcfs.map { it -> it[0] }.collect()
    indices = somatic_merged_vcfs.map { it -> it[1] }.collect()
    bcftoolsConcatSomaticCandidates(vcfs, indices)

    // Extract somatic candidates from DB
    somatic_ch = ref_files.combine(ivls.combine(somatic_db, by: 0))
        .map { fa, fai, dict, interval_label, interval, genomeDB ->
            tuple(interval_label, [fa, fai, dict], interval, genomeDB) }
    somatic_candidates_by_interval = normalizeSomaticCandidates(extractSomaticCandidates(somatic_ch))
    candidate_vcfs = somatic_candidates_by_interval.map { it -> it[0] }.collect()
    candidate_indices = somatic_candidates_by_interval.map { it -> it[1] }.collect()
    somatic_candidates = mergeSomaticCandidates(candidate_vcfs, candidate_indices)                          

    // Run HaplotypeCaller on normals
    gvcfs_ch = runHaplotypeCallerOnNormal(normals_for_calling_ch)

    // Collect the gvcfs into a genomeDB
    dbImport_ch = gvcfs_ch.groupTuple().combine(ivls, by: 0)
    db = genomicsDBImport(dbImport_ch)

    // Genotype
    to_genotype_ch = ref_files.combine(ivls.combine(db, by: 0))
        .map { fa, fai, dict, interval_label, interval, genomeDB ->
            tuple(interval_label, [fa, fai, dict], interval, genomeDB) }
    genotyped_vcfs = genotypeGVCFs(to_genotype_ch)

    // Merge the genotyped VCFs
    vcfs = genotyped_vcfs.map { it -> it[0] }.collect()
    indices = genotyped_vcfs.map { it -> it[1] }.collect()
    germline_resource = mergeGenotypedVCFs(vcfs, indices)

    // Finalize the Somatic candidates by adding variants from the PON and GR
    to_finalize_ch = ref_files
        .combine(germline_resource
                 .concat(panel_of_normals)
                 .concat(somatic_candidates)
                 .collate(3))
        .map { fa, fai, stats, gr, pon, sc ->
            tuple([fa, fai, stats], gr, pon, sc ) }
    
    finalizeSomaticCandidates(to_finalize_ch)
}
