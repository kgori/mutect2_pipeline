nextflow.enable.dsl=2

params.reference        = "genome.fa"
params.intervals        = 50
params.panel_normals    = "panel_of_normals_sampes_folder"
params.resource_normals = "germline_resource_samples_folder"
params.tumours          = "somatic_tumour_samples_folder"
params.outdir           = "results"

def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}

def make_db_ch(filepairs_ch, ivls, ref_files) {
    vcf_files = filepairs_ch
        .map { _id, files -> tuple("vcf", files.find { it.name.endsWith(".vcf.gz") }) }
        .groupTuple()
    index_files = filepairs_ch
        .map { _id, files -> tuple("tbi", files.find { it.name.endsWith(".vcf.gz.tbi") }) }
        .groupTuple()
    return ref_files.combine(vcf_files.combine(index_files).combine(ivls))
        .map { fa, fai, dict, _vcf, vcfs, _tbi, tbis, interval_id, interval ->
            tuple(interval_id, [fa, fai, dict], vcfs, tbis, interval) }
}

include { indexReference }                           from "./preprocessReference.nf"
include { simpleSplitIntervals }                     from "./preprocessReference.nf"
include { makeReferenceDict }                        from "./preprocessReference.nf"
include { genotypeGvcfIntervals }                    from "./stage2_supporting.nf"
include { makePonIntervals }                         from "./stage2_supporting.nf"
include { makeSomaticCandidatesIntervals }           from "./stage2_supporting.nf"
// include { genomicsDBImport_PON }                     from "./genomicsDB.nf"
// include { genomicsDBImport_Somatic }                 from "./genomicsDB.nf"
// include { genotypeGVCFs }                            from "./germlineResource.nf"
// include { mergeGenotypedVCFs }                       from "./germlineResource.nf"
// include { createPanelOfNormals }                     from "./panelOfNormals.nf"
// include { mergePonVCFs }                             from "./panelOfNormals.nf"
// include { extractSomaticCandidates }                 from "./somaticCandidates.nf"
// include { normalizeSomaticCandidates }               from "./somaticCandidates.nf"
// include { mergeSomaticCandidates }                   from "./somaticCandidates.nf"
// include { finalizeSomaticCandidates }                from "./somaticCandidates.nf"
// include { callSomaticVariants }                      from "./finalVariantCalling.nf"


workflow {
    /////////////////////////////////
    //         Preanalysis setup
    //

    // Load and process the reference
    ref_ch = Channel.fromPath(params.reference)
    fai_ch = indexReference(ref_ch)
    dict_ch = makeReferenceDict(ref_ch)
    ref_files = ref_ch.merge(fai_ch).merge(dict_ch)

    // Construct intervals
    split_intervals = simpleSplitIntervals(ref_files, params.intervals)
    ivls = split_intervals.flatten()
        .map { interval ->
            def intervalNumberMatch = interval.getName() =~ /^(\d+)/
                def intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
            return tuple(intervalNumber, interval) }

    ///////////////////////////////////////////////////////
    //         Stage 2: Create panels and candidates
    // Load the germline resource normals
    gvcf_filepairs = Channel.fromFilePairs("${params.resource_normals}/*.{vcf.gz,vcf.gz.tbi}", checkIfExists: true)
    pon_filepairs = Channel.fromFilePairs("${params.panel_normals}/*.{vcf.gz,vcf.gz.tbi}", checkIfExists: true)
    som_filepairs = Channel.fromFilePairs("${params.tumours}/*.{vcf.gz,vcf.gz.tbi}", checkIfExists: true)

    // Construct appropriate nextflow channels for stage 2
    gvcf_ch = make_db_ch(gvcf_filepairs, ivls, ref_files)
    pon_ch = make_db_ch(pon_filepairs, ivls, ref_files)
    som_ch = make_db_ch(som_filepairs, ivls, ref_files)

    // Run on each interval
    genotyped_gvcfs = genotypeGvcfIntervals(gvcf_ch)
    // pon_intervals = makePonIntervals(pon_ch)
    // somatic_candidates_intervals = makeSomaticCandidatesIntervals(som_ch)

    // Merge interval results
    
    // Merge the Panel VCFs
    // vcfs = panel_of_normals_vcfs.map { it -> it[0] }.collect()
    // indices = panel_of_normals_vcfs.map { it -> it[1] }.collect()
    // panel_of_normals = mergePonVCFs(vcfs, indices)

    // candidate_vcfs = somatic_candidates_by_interval.map { it -> it[0] }.collect()
    // candidate_indices = somatic_candidates_by_interval.map { it -> it[1] }.collect()
    // somatic_candidates = mergeSomaticCandidates(candidate_vcfs, candidate_indices)

    // vcfs = genotyped_vcfs.map { it -> it[0] }.collect()
    // indices = genotyped_vcfs.map { it -> it[1] }.collect()
    // germline_resource = mergeGenotypedVCFs(vcfs, indices)

    // Finalize the Somatic candidates by adding variants from the PON and GR
    // collated_support_vcfs_ch = germline_resource
    //     .concat(panel_of_normals)
    //     .concat(somatic_candidates)
    //     .collate(3)
    // to_finalize_ch = ref_files
    //     .combine(collated_support_vcfs_ch)
    //     .map { fa, fai, stats, gr, pon, sc ->
    //         tuple([fa, fai, stats], gr, pon, sc ) }

    // candidates = finalizeSomaticCandidates(to_finalize_ch)
}
