nextflow.enable.dsl=2

params.reference          = "genome.fa"
params.normals            = "normals_folder"
params.tumours            = "tumours_folder"
params.samplesheet        = "samplesheet.csv"
params.germline_resource  = "germline_resource.vcf.gz"
params.panel_of_normals   = "panel_of_normals.vcf.gz"
params.somatic_candidates = "somatic_candidates.vcf.gz"
params.intervals          = 50
params.outdir             = "results"

include { indexReference }                                     from "./preprocessReference.nf"
include { simpleSplitIntervals }                               from "./preprocessReference.nf"
include { makeReferenceDict }                                  from "./preprocessReference.nf"
include { callMatchedSomaticVariants }                         from "./finalVariantCalling.nf"
include { callSomaticVariants }                                from "./finalVariantCalling.nf"
include { callSomaticVariants as callSomaticVariantsOnNormal } from "./finalVariantCalling.nf"
include { recallGermlineVariants }                             from "./finalVariantCalling.nf"
include { runContaminationModelOnPaired }                      from "./contaminationModelling.nf"
include { runContaminationModelOnTumourOnly }                  from "./contaminationModelling.nf"
include { filterMutectCalls }                                  from "./filterModelling.nf"
include { filterMutectCallsOnNormals }                         from "./filterModelling.nf"
include { learnReadOrientationModel }                          from "./filterModelling.nf"
include { concatFilteredCalls }                                from "./vcfConcatenator.nf"
include { concatSecondHaplotypeCallerCalls }                   from "./vcfConcatenator.nf"


def get_id_from_filename = { f ->
    f.getBaseName().replaceFirst(/\.(bam|cram)$/, "")
}

workflow {
    /////////////////////////////////
    //         Preanalysis setup
    //

    // Load and process the reference
    ref_ch = Channel.fromPath(params.reference)
    fai_ch = indexReference(ref_ch)
    dict_ch = makeReferenceDict(ref_ch)
    ref_files = ref_ch.merge(fai_ch).merge(dict_ch)

    // Chop into intervals for scattering-gathering
    split_intervals = simpleSplitIntervals(ref_files, params.intervals)
    ivls = split_intervals.flatten()
        .map { interval ->
            def intervalNumberMatch = interval.getName() =~ /^(\d+)/
                def intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
            return tuple(intervalNumber, interval) }
    
    // Load data  /// TODO: Column1 is redundant and can be removed
    // Parse samplesheet and resolve file paths
    // Reminder of the input format: the samplesheet has three columns and no header.
    // - Column 1 = sample name (any unique ID for the samples in the row)
    // - Column 2 = tumour BAM/CRAM file name (relative to params.tumours), or blank if no tumour
    // - Column 3 = normal BAM/CRAM file name (relative to params.normals), or blank if no normal
    // The following channel constructor extracts data according to the samplesheet. It branches into three channels:
    // - paired: rows with both tumour and normal [pair_id, tumour_id, [tumour_cram, tumour_crai], normal_id, [normal_cram, normal_crai]]
    // - tumour_only: rows with only a tumour [tumour_id, [tumour_cram, tumour_crai]]
    // - normal_only: rows with only a normal [normal_id, [normal_cram, normal_crai]]
    inputs = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: false)
        .filter { row -> row.any { it?.trim() } }           // skip empty rows
        .map { row ->
            def sample   = row[0]?.trim()
            def t_name   = row[1]?.trim()
            def n_name   = row[2]?.trim()

            def tumour = null
            def normal = null
            def tumour_id = null
            def normal_id = null

            if (t_name) {
                def cram = file(params.tumours).resolve(t_name)
                def crai
                if (cram.name.endsWith("cram")) crai = file("${cram}.crai")
                else if (cram.name.endsWith("bam")) crai = file("${cram}.bai")
                else error "Alignment files must end with .bam or .cram: ${cram}"
                    
                if (!cram.exists()) error "Tumour CRAM not found for ${sample}: ${cram}"
                if (!crai.exists()) error "Tumour CRAI not found for ${sample}: ${crai}"
                tumour = [cram, crai]
                tumour_id = get_id_from_filename(cram)
            }

            if (n_name) {
                def cram = file(params.normals).resolve(n_name)
                def crai
                if (cram.name.endsWith("cram")) crai = file("${cram}.crai")
                else if (cram.name.endsWith("bam")) crai = file("${cram}.bai")
                else error "Alignment files must end with .bam or .cram: ${cram}"

                if (!cram.exists()) error "Normal CRAM not found for ${sample}: ${cram}"
                if (!crai.exists()) error "Normal CRAI not found for ${sample}: ${crai}"
                normal = [cram, crai]
                normal_id = get_id_from_filename(cram)
            }

            tuple(sample, tumour_id, tumour, normal_id, normal)
        }
        .branch {
            paired:      it[2] != null && it[4] != null
            tumour_only: it[2] != null
            normal_only: it[4] != null
        }

    // A little bit of fine-tuning to the channels. Also construct channels for all tumours and all normals, which will be used when
    // building filtering models
    inputs.paired
        .map { sample, tumour_id, tumour, normal_id, normal ->
            def pair_id = "${tumour_id}_vs_${normal_id}"
            tuple(pair_id, tumour_id, tumour, normal_id, normal) }
        .set { paired }
    inputs.tumour_only
        .map { sample, tumour_id, tumour, _, _2 -> tuple(tumour_id, tumour) }
        .set { tumour_only }
    inputs.normal_only
        .map { sample, _, _2, normal_id, normal -> tuple(normal_id, normal) }
        .set { normal_only }
    all_normals = paired
        .map { sample, _, _2, normal_id, normal -> tuple(normal_id, normal) }
        .concat(normal_only)
        .unique { it[0] }
    all_tumours = paired
        .map { pair_id, tumour_id, tumour, _, _2 -> tuple(tumour_id, tumour) }
        .concat(tumour_only)
        .unique { it[0] } 

    // Combine inputs with intervals, for scattering+gathering
    paired_interval_ch = paired.combine(ivls)
    tumour_interval_ch = tumour_only.combine(ivls)
    normal_interval_ch = normal_only.combine(ivls)
    all_normals_interval_ch = all_normals.combine(ivls)

    // Load the panel, germline resource, and somatic candidates
    panel_of_normals = Channel
        .fromFilePairs("${params.panel_of_normals}{,.tbi}", size: 2, checkIfExists: true)
        .map { label, files -> files }
    germline_resource = Channel
        .fromFilePairs("${params.germline_resource}{,.tbi}", size: 2, checkIfExists: true)
        .map { label, files -> files }
    somatic_candidates = Channel
        .fromFilePairs("${params.somatic_candidates}{,.tbi}", size: 2, checkIfExists: true)
        .map { label, files -> files }


    ///////////////////////////////////////////////////////
    //         Stage 3: Final calling and filtering
    //

    // Call somatic variants (round 2 - final calling step before filtering)
    // Run mutect2 in paired mode
    paired_calling_ch = paired_interval_ch
        .combine(ref_files)
        .combine(germline_resource)
        .combine(panel_of_normals)
        .combine(somatic_candidates)
        .map { pair_id, tumour_id, tumour_bam, normal_id, normal_bam, interval_id,
              intervals, fa, fai, dict, germline_resource, germline_resource_index,
              panel_of_normals, panel_of_normals_index, candidates, candidates_index ->
            tuple(interval_id, [fa, fai, dict], tumour_id,
                  tumour_bam.toList(), normal_bam.toList(),
                  intervals,
                  [germline_resource, germline_resource_index],
                  [panel_of_normals, panel_of_normals_index],
                  [candidates, candidates_index]) }
    
    paired_calls_ch = callMatchedSomaticVariants(paired_calling_ch)

    // Run mutect2 in tumour only mode
    tumour_only_calling_ch = tumour_interval_ch
        .combine(ref_files)
        .combine(germline_resource)
        .combine(panel_of_normals)
        .combine(somatic_candidates)
        .map { sample, tumour_bam, interval_id,
              intervals, fa, fai, dict, germline_resource, germline_resource_index,
              panel_of_normals, panel_of_normals_index, candidates, candidates_index ->
            tuple(interval_id, [fa, fai, dict], sample,
                  tumour_bam.toList(),
                  intervals,
                  [germline_resource, germline_resource_index],
                  [panel_of_normals, panel_of_normals_index],
                  [candidates, candidates_index]) }

    tumour_only_calls_ch = callSomaticVariants(tumour_only_calling_ch)

    // Run mutect2 on normals, too
    all_normals_calling_ch = all_normals_interval_ch
        .combine(ref_files)
        .combine(germline_resource)
        .combine(panel_of_normals)
        .combine(somatic_candidates)
        .map { sample, normal_bam, interval_id,
              intervals, fa, fai, dict, germline_resource, germline_resource_index,
              panel_of_normals, panel_of_normals_index, candidates, candidates_index ->
            tuple(interval_id, [fa, fai, dict], sample,
                  normal_bam.toList(),
                  intervals,
                  [germline_resource, germline_resource_index],
                  [panel_of_normals, panel_of_normals_index],
                  [candidates, candidates_index]) }

    all_normals_calls_ch = callSomaticVariantsOnNormal(all_normals_calling_ch)

    // Rerun haplotype caller on normals at candidate sites
    rehaplotyped_normals = recallGermlineVariants(all_normals_calling_ch)

    // Build a strand bias model
    orientation_model_ch = paired_calls_ch.f1r2s
        .mix(tumour_only_calls_ch.f1r2s)
        .mix(all_normals_calls_ch.f1r2s)
        .groupTuple(size: params.intervals)
        .map { key, files -> tuple(key, files.sort { it.name }) }
    orientation_model = learnReadOrientationModel(orientation_model_ch)

    // Build a contamination model
    paired_contamination_ch = paired
        .combine(ref_files)
        .combine(germline_resource)
        .map { pair_id, tumour_id, tumour_bam, normal_id, normal_bam, 
               fa, fai, dict, germline_resource, germline_resource_index ->
            tuple(tumour_id, [fa, fai, dict],
                  tumour_bam.toList(),
                  normal_bam.toList(),
                  [germline_resource, germline_resource_index]) }

    run_contamination_model_on_paired_ch = runContaminationModelOnPaired(paired_contamination_ch)
    
    unpaired_contamination_ch = tumour_only
        .combine(ref_files)
        .combine(germline_resource)
        .map { tumour_id, tumour, fa, fai, dict, germline_resource, germline_resource_index ->
            tuple(tumour_id, [fa, fai, dict],
                  tumour.toList(),
                  [germline_resource, germline_resource_index]) }

    run_contamination_model_on_tumour_only_ch = runContaminationModelOnTumourOnly(unpaired_contamination_ch)
    
    // Filter calls
    tumour_calls = paired_calls_ch.vcfs
        .mix(tumour_only_calls_ch.vcfs)

    tumour_contamination = run_contamination_model_on_paired_ch
        .mix(run_contamination_model_on_tumour_only_ch)

    tumour_filter_ch = ref_files
        .combine(tumour_calls)
        .map { fa, fai, dict, sample, vcf, tbi, stats, interval_id, intervals ->
            tuple(sample, [fa, fai, dict], [vcf, tbi], stats, interval_id, intervals) }
        .combine(tumour_contamination, by: 0) 
        .combine(orientation_model, by: 0)

    filtered = filterMutectCalls(tumour_filter_ch)

    normal_calls = all_normals_calls_ch.vcfs
    normal_filter_ch = ref_files
        .combine(normal_calls)
        .map { fa, fai, dict, sample, vcf, tbi, stats, interval_id, intervals ->
            tuple(sample, [fa, fai, dict], [vcf, tbi], stats, interval_id, intervals) }
        .combine(orientation_model, by: 0)

    filtered_normals = filterMutectCallsOnNormals(normal_filter_ch)

    // Merge final variants
    grouped = filtered.mix(filtered_normals).groupTuple(size: params.intervals)
        .map { label, vcfs, indices -> tuple(label, vcfs.sort { it.name }, indices.sort { it.name }) }
    grouped_with_ref = ref_files.combine(grouped)
        .map { fa, fai, dict, sample, vcfs, tbis ->
            tuple(sample, [fa, fai, dict], vcfs, tbis) }
    grouped_with_pon = panel_of_normals.combine(grouped_with_ref)
        .map { pon_vcf, pon_tbi, sample, ref, vcfs, tbis ->
            tuple(sample, ref, [pon_vcf, pon_tbi], vcfs, tbis) }
    concat = concatFilteredCalls(grouped_with_pon)

    // And the normals
    grouped_rehaplotyped = rehaplotyped_normals
        .map { it -> tuple(it[0], it[1], it[2]) } // sample, vcf, tbi
        .groupTuple(size: params.intervals)
        .map { label, vcfs, indices -> tuple(label, vcfs.sort { it.name }, indices.sort { it.name }) }
    grouped_rehaplotyped_with_ref = ref_files.combine(grouped_rehaplotyped)
        .map { fa, fai, dict, sample, vcfs, tbis ->
            tuple(sample, [fa, fai, dict], vcfs, tbis) }

    concat_rehaplotyped = concatSecondHaplotypeCallerCalls(grouped_rehaplotyped_with_ref)
}
