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

process makeBamToSampleNameMap {
    input:
    path(bam_files)

    output:
    path("bam_to_sample_name_map.txt")

    script:
    """
    for f in ${bam_files.join(' ')}; do
        samplename=\$(samtools view -H "\$f" | awk '
    {
      for (i=1; i<=NF; i++)
        if (\$i ~ /^SM:/)
        {
            split(\$i, a, ":");
            print a[2];
            exit
        }
    }')
        filename=\$(basename "\$f")
        filename=\${filename%.*}
        printf "%s\t%s\n" \$filename \$samplename
    done >> bam_to_sample_name_map.txt
    """
}

process splitIntervals {
  input:
  path(reference)
  val(numIntervals)

  output:
  path("Intervals/*.interval_list")

  publishDir "results", mode: 'copy'

  script:
  """
  gatk SplitIntervals \
    --reference "${reference[0]}" \
    --scatter-count $numIntervals \
    --intervals 38 \
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

process runMutectOnNormal {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(normal_bam), path(interval)

    output:
    tuple val(interval_id), path("${sample}.*.normal.mutect2_panel_calls.vcf.gz"),
        path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.tbi"),
        path("${sample}.*.normal.mutect2_panel_calls.vcf.gz.stats")

    script:
    intervalNumberMatch = interval.getName() =~ /^(\d+)/
    intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
    """
    gatk Mutect2 \
      --reference ${reference[0]} \
      --input ${normal_bam[0]} \
      --intervals ${interval} \
      --interval-padding 150 \
      --native-pair-hmm-threads ${task.cpus} \
      --max-mnp-distance 0 \
      --output ${sample}.${intervalNumber}.normal.mutect2_panel_calls.vcf.gz
    """
}

process runMutectOnTumour {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(tumour_bam), path(interval)

    output:
    tuple val(interval_id),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz"),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.tbi"),
        path("${sample}.*.tumour.mutect2_candidate_discovery_calls.vcf.gz.stats")

    publishDir "results/InitialTumourCalls", mode: 'copy'
    
    script:
    intervalNumberMatch = interval.getName() =~ /^(\d+)/
    intervalNumber = intervalNumberMatch ? intervalNumberMatch[0][1] : 99999
    """
    gatk Mutect2 \
      --reference ${reference[0]} \
      --input ${tumour_bam[0]} \
      --intervals ${interval} \
      --interval-padding 150 \
      --native-pair-hmm-threads ${task.cpus} \
      --max-mnp-distance 0 \
      --output ${sample}.${intervalNumber}.tumour.mutect2_candidate_discovery_calls.vcf.gz
    """
}

process runHaplotypeCallerOnNormal {
    cpus 4
    memory '8 GB'
    time '12h'
    queue 'normal'
    executor 'lsf'
    
    input:
    tuple val(interval_id), path(reference), val(sample), path(normal_bam), path(interval)

    output:
    tuple val(interval_id),
        path("${sample}.${interval_id}.haplotypecaller.g.vcf.gz"),
        path("${sample}.${interval_id}.haplotypecaller.g.vcf.gz.tbi")

    script:
    """
    gatk HaplotypeCaller \
      --reference ${reference[0]} \
      --input ${normal_bam[0]} \
      --intervals ${interval} \
      --interval-padding 150 \
      --native-pair-hmm-threads ${task.cpus} \
      --emit-ref-confidence GVCF \
        --output ${sample}.${interval_id}.haplotypecaller.g.vcf.gz
    """
}

process genomicsDBImport {
    input:
    tuple val(interval_id), path(gvcfs), path(gvcf_index), path(interval)

    output:
    tuple val(interval_id), path("resourceDB")

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path resourceDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --intervals ${interval} \
      -V ${gvcfs.join(' -V ')}
    """
}

process genomicsDBImport_PON {
    input:
    tuple val(interval_id), path(vcfs), path(vcf_index), path(statsfiles), path(interval)

    output:
    tuple val(interval_id), path("ponDB")

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path ponDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --intervals ${interval} \
      -V ${vcfs.join(' -V ')}
    """
}

process genomicsDBImport_Somatic {
    input:
    tuple val(interval_id), path(vcfs), path(vcf_index), path(statsfiles), path(interval)

    output:
    tuple val(interval_id), path("somaticDB")

    script:
    """
    gatk GenomicsDBImport \
      --genomicsdb-workspace-path somaticDB \
      --batch-size 50 \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --intervals ${interval} \
      -V ${vcfs.join(' -V ')}
    """
}

process createPanelOfNormals {
    input:
    tuple val(interval_id), path(reference), path(interval), path(genomeDB)

    output:
    tuple path("${interval_id}.panel_of_normals.vcf.gz"),
        path("${interval_id}.panel_of_normals.vcf.gz.tbi")

    script:
    """
    gatk CreateSomaticPanelOfNormals \
        --reference ${reference[0]} \
        --variant gendb://${genomeDB} \
        --intervals ${interval} \
        --output ${interval_id}.panel_of_normals.vcf.gz
    """
}

process extractSomaticCandidates {
    input:
    tuple val(interval_id), path(reference), path(interval), path(genomeDB)

    output:
    tuple val(interval_id), path(reference), path("${interval_id}.somatic_candidates_raw.vcf.gz"),
        path("${interval_id}.somatic_candidates_raw.vcf.gz.tbi")

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

    script:
    """
    bcftools view -e 'ALT="*"' ${vcf} \
        | bcftools norm -m -both -f ${reference[0]} \
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
    bcftools index somatic_candidates.vcf.gz
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

    script:
    """
    bcftools merge ${vcfs.sort { it.getName() }.join(' ')} \
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
    bcftools index bcftoolsSomaticCandidates.vcf.gz
    """
}
    

process genotypeGVCFs {
    input:
    tuple val(interval_id), path(reference), path(interval), path(genomeDB)

    output:
    tuple path("${interval_id}.genotyped.vcf.gz"),
        path("${interval_id}.genotyped.vcf.gz.tbi")

    script:
    """
    gatk GenotypeGVCFs \
        --reference ${reference[0]} \
        --variant gendb://${genomeDB} \
        --intervals ${interval} \
        --output ${interval_id}.genotyped.vcf.gz
    """
}

process mergeGenotypedVCFs {
    executor 'local'
    input:
    path(genotyped_vcfs)
    path(vcf_indices)

    output:
    path("genome_resource.vcf.gz*")

    publishDir "results/GermlineResource", mode: 'copy'
    
    script:
    """
    bcftools concat -a -D ${genotyped_vcfs.join(' ')} \
        | bcftools sort -Oz -o genome_resource.vcf.gz
    bcftools index genome_resource.vcf.gz
    """
}

process mergePonVCFs {
    executor 'local'
    input:
    path(pon_vcfs)
    path(vcf_indices)

    output:
    path("panel_of_normals.vcf.gz*")

    publishDir "results/PON", mode: 'copy'

    script:
    """
    bcftools concat -a -D ${pon_vcfs.join(' ')} \
        | bcftools sort -Oz -o panel_of_normals.vcf.gz
    bcftools index panel_of_normals.vcf.gz
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
    mergePonVCFs(vcfs, indices)

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
    mergeSomaticCandidates(candidate_vcfs, candidate_indices)                          

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
    mergeGenotypedVCFs(vcfs, indices)
}
