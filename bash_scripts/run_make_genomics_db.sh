#!/bin/bash
PROJ_DIR="/lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2"

VCF_FOLDER=${1}
REFERENCE=${2}
INTERVALS=${3}
OUTPUT=${4}

if [[ -z "${VCF_FOLDER}" || -z "${REFERENCE}" || -z "${INTERVALS}" || -z "${OUTPUT}" ]]; then
  echo "Usage: $0 <vcf_folder:Path> <reference:Fasta file> <intervals:GATK Interval list> <output:GenomicsDB path>"
  exit 1
fi

if [ ! -d "${VCF_FOLDER}" ]; then
  echo "Error: VCF folder not found: ${VCF_FOLDER}"
  exit 1
fi

# Check that $VCF_FOLDER contains at least one file matching the wildcard *.vcf.gz
shopt -s nullglob
files=("${VCF_FOLDER}"/*.vcf.gz)
if ! (( ${#files[@]} )); then
  echo "Error: No VCF files found in the specified folder: ${VCF_FOLDER}"
  exit 1
fi
shopt -u nullglob

# Resolve full path to each VCF file. Then combine into a string of "-V file1 -V file2 ..."
inputs=""
for file in "${files[@]}"; do
  full_path=$(realpath "${file}")
  inputs+=" -V ${full_path}"
done

if [ ! -f "${REFERENCE}" ]; then
  echo "Error: Reference file not found: ${REFERENCE}"
  exit 1
fi

if [ ! -f "${REFERENCE}.fai" ]; then
  echo "Error: Reference index file not found: ${REFERENCE}.fai"
  exit 1
fi

if [ ! -f "${INTERVALS}" ]; then
  echo "Error: Intervals file not found: ${INTERVALS}"
  exit 1
fi

if [ -d "${OUTPUT}" ]; then
  echo "Error: Output path already exists and is a directory: ${OUTPUT}"
  exit 1
fi

REFERENCE=$(realpath "${REFERENCE}")
INTERVALS=$(realpath "${INTERVALS}")

singularity exec \
  -B '/nfs:/nfs' \
  -B '/lustre:/lustre' \
  --no-home \
  -e ${PROJ_DIR}/container/gatk_latest.sif \
  gatk GenomicsDBImport \
      -R ${REFERENCE} \
      ${inputs} \
      -L ${INTERVALS} \
      --genomicsdb-shared-posixfs-optimizations true \
      --bypass-feature-reader true \
      --genomicsdb-workspace-path ${OUTPUT} 
