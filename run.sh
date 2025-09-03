#!/bin/bash

nextflow run main.nf \
         --reference /lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2/refs/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
         --normals /lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2/bams/normals \
         --tumours /lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2/bams/tumours \
         --outdir /lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2/results_16_samples_20250903 \
         -resume
