#!/bin/bash
PROJ_DIR="/lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2"

BAM=${1}
REFERENCE=${2}
OUTPUT_PATH=${3}
INTERVALS=${4}

if [[ -z "${BAM}" || -z "${REFERENCE}" || -z "${OUTPUT_PATH}" ]]; then
  echo "Usage: $0 <alignment:BAM file> <reference:Fasta file> <output:Path> [intervals file]"
  exit 1
fi

# resolve the full path to the BAM file
BAM=$(realpath "${BAM}")

if [ ! -f "${BAM}" ]; then
  echo "Error: BAM file not found: ${BAM}"
  exit 1
fi

if [ ! -z "${INTERVALS}" ]; then
  IVLS_FILE=$(realpath "${INTERVALS}")
  if [ -f "${IVLS_FILE}" ]; then
    INTERVALS="-L ${IVLS_FILE}"
    echo "Using intervals file: ${IVLS_FILE}"
  else
    echo "Warning: Intervals file not found: ${IVLS_FILE}. Proceeding without intervals."
    INTERVALS=""
  fi
else
  echo "No intervals file provided. Proceeding without intervals."
  INTERVALS=""
fi

# Create the output directory if it doesn't exist
if [ -f "${OUTPUT_PATH}" ]; then
  echo "Error: Output path is a file, not a directory: ${OUTPUT_PATH}"
  exit 1
fi

mkdir -p "${OUTPUT_PATH}"
if [ ! -d "${OUTPUT_PATH}" ]; then
  echo "Error: Unable to create output directory: ${OUTPUT_PATH}"
  exit 1
fi

# Create an output filename from the BAM file name
BAM_NAME=$(basename "${BAM}")
# Remove any .bam or .cram extension
BAM_NAME=${BAM_NAME%.bam}
BAM_NAME=${BAM_NAME%.cram}

singularity exec \
  -B '/nfs:/nfs' \
  -B '/lustre:/lustre' \
  --no-home \
  -e ${PROJ_DIR}/container/gatk_latest.sif \
  gatk Mutect2 \
      --reference ${REFERENCE} \
      --input ${BAM} \
      --max-mnp-distance 0 \
      --output ${OUTPUT_PATH}/${BAM_NAME}.mutect2.vcf.gz \
      $INTERVALS \
      --interval-padding 150

