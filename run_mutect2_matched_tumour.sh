#!/bin/bash
PROJ_DIR="/lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2"

TUMOUR_BAM=${1}
NORMAL_BAM=${2}
REFERENCE=${3}
GERMLINE_RESOURCE=${4}
PANEL_OF_NORMALS=${5}
SOMATIC_CANDIDATES=${6}
INTERVALS=${7}
OUTPUT_PATH=${8:-$PROJ_DIR}

if [[ -z "${TUMOUR_BAM}" || -z "${NORMAL_BAM}" || -z "${REFERENCE}" || -z "${GERMLINE_RESOURCE}" || -z "${PANEL_OF_NORMALS}" || -z "${SOMATIC_CANDIDATES}" ]]; then
  echo "Usage: $0 <tumour:BAM file> <normal:BAM file> <reference:Fasta file> <germline_resource:VCF file> <panel_of_normals:VCF file> <somatic_candidates:VCF file> [intervals file]"
  exit 1
fi

if [ ! -f "${TUMOUR_BAM}" ]; then
  echo "Error: BAM file not found: ${TUMOUR_BAM}"
  exit 1
fi

if [ ! -f "${NORMAL_BAM}" ]; then
  echo "Error: BAM file not found: ${NORMAL_BAM}"
  exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
  echo "Error: Reference file not found: ${REFERENCE}"
  exit 1
fi

if [ ! -f "${GERMLINE_RESOURCE}" ]; then
  echo "Error: Germline resource file not found: ${GERMLINE_RESOURCE}"
  exit 1
fi

if [ ! -f "${PANEL_OF_NORMALS}" ]; then
  echo "Error: Panel of normals file not found: ${PANEL_OF_NORMALS}"
  exit 1
fi

if [ ! -f "${SOMATIC_CANDIDATES}" ]; then
  echo "Error: Somatic candidates file not found: ${SOMATIC_CANDIDATES}"
  exit 1
fi

# Create the output path if it doesn't exist
mkdir -p "${OUTPUT_PATH}"
if [ ! -d "${OUTPUT_PATH}" ]; then
  echo "Error: Unable to create output directory: ${OUTPUT_PATH}"
  exit 1
fi

# resolve the full path to the input files
TUMOUR_BAM=$(realpath "${TUMOUR_BAM}")
REFERENCE=$(realpath "${REFERENCE}")
NORMAL_BAM=$(realpath "${NORMAL_BAM}")
GERMLINE_RESOURCE=$(realpath "${GERMLINE_RESOURCE}")
PANEL_OF_NORMALS=$(realpath "${PANEL_OF_NORMALS}")
SOMATIC_CANDIDATES=$(realpath "${SOMATIC_CANDIDATES}")

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

# Create an output filename from the BAM file name
BAM_NAME=$(basename "${TUMOUR_BAM}")
# Remove any .bam or .cram extension
BAM_NAME=${BAM_NAME%.bam}
BAM_NAME=${BAM_NAME%.cram}

# Extract the normal sample name from the normal BAM file
NORMAL_NAME=$(samtools view -H ${NORMAL_BAM} | grep '^@RG' | sed -n 's/.*SM:\([^ \t]*\).*/\1/p' | sort -u)

singularity exec \
  -B '/nfs:/nfs' \
  -B '/lustre:/lustre' \
  -Ce ${PROJ_DIR}/container/gatk_latest.sif \
  gatk Mutect2 \
      --reference ${REFERENCE} \
      --input ${TUMOUR_BAM} \
      --input ${NORMAL_BAM} \
      --normal-sample ${NORMAL_NAME} \
      --germline-resource ${GERMLINE_RESOURCE} \
      --panel-of-normals ${PANEL_OF_NORMALS} \
      --genotype-germline-sites \
      --genotype-pon-sites \
      --alleles ${SOMATIC_CANDIDATES} \
      --f1r2-tar-gz ${OUTPUT_PATH}/${BAM_NAME}.mutect2.f1r2.tar.gz \
      --output ${OUTPUT_PATH}/${BAM_NAME}.mutect2.unfiltered.vcf.gz \
      $INTERVALS \
      --interval-padding 150
