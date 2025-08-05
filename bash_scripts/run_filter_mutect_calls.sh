#!/bin/bash
PROJ_DIR="/lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2"

REFERENCE=${1}
CALLS_VCF=${2}
CONTAMINATION_TABLE=${3}
ORIENTATION_BIAS=${4}
OUTPUT_PATH=${5}

if [[ -z "${REFERENCE}" || -z "${CALLS_VCF}" || -z "${CONTAMINATION_TABLE}" || -z "${ORIENTATION_BIAS}" ]]; then
  echo "Usage: $0 <reference:Fasta file> <calls:VCF file> <contamination_table:TSV file> <orientation_bias:TSV file> <output_path>"
  exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
  echo "Error: Reference file not found: ${REFERENCE}"
  exit 1
fi

if [ ! -f "${CALLS_VCF}" ]; then
  echo "Error: Calls VCF file not found: ${CALLS_VCF}"
  exit 1
fi

if [ ! -f "${CONTAMINATION_TABLE}" ]; then
  echo "Error: Contamination table file not found: ${CONTAMINATION_TABLE}"
  exit 1
fi

if [ ! -f "${ORIENTATION_BIAS}" ]; then
  echo "Error: Orientation bias file not found: ${ORIENTATION_BIAS}"
  exit 1
fi

REFERENCE=$(realpath "${REFERENCE}")
CALLS_VCF=$(realpath "${CALLS_VCF}")
CONTAMINATION_TABLE=$(realpath "${CONTAMINATION_TABLE}")
ORIENTATION_BIAS=$(realpath "${ORIENTATION_BIAS}")

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_PATH}"
if [ ! -d "${OUTPUT_PATH}" ]; then
  echo "Error: Unable to create output directory: ${OUTPUT_PATH}"
  exit 1
fi

OUT_NAME=$(basename "${CALLS_VCF}" .vcf.gz)
OUT_NAME=${OUT_NAME%.vcf}

singularity exec \
  -B '/nfs:/nfs' \
  -B '/lustre:/lustre' \
  --no-home \
  --env JAVA_OPTION='-Xmx16g' \
  -e ${PROJ_DIR}/container/gatk_latest.sif \
  gatk FilterMutectCalls \
    -R "${REFERENCE}" \
    -V "${CALLS_VCF}" \
    --contamination-table "${CONTAMINATION_TABLE}" \
    --ob-priors "${ORIENTATION_BIAS}" \
    -O "${OUTPUT_PATH}/${OUT_NAME}.filtered.vcf.gz" \
    --create-output-variant-index true
