#!/bin/bash
PROJ_DIR="/lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2"

TUMOUR_BAM=${1}
REFERENCE=${2}
GERMLINE_RESOURCE=${3}
PANEL_OF_NORMALS=${4}
GENOTYPING_CANDIDATES=${5}
INTERVALS=${6}
OUTPUT_PATH=${7}

if [[ -z "${REFERENCE}" || -z "${TUMOUR_BAM}" || -z "${GERMLINE_RESOURCE}" || -z "${PANEL_OF_NORMALS}" || -z "${GENOTYPING_CANDIDATES}" ]]; then
  echo "Usage: $0 <tumour:BAM file> <reference:Fasta file> <germline_resource:VCF file> <panel_of_normals:VCF file> <somatic_candidates:VCF file> [intervals file]"
  exit 1
fi

if [ ! -f "${TUMOUR_BAM}" ]; then
  echo "Error: BAM file not found: ${TUMOUR_BAM}"
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

if [ ! -f "${GENOTYPING_CANDIDATES}" ]; then
  echo "Error: Somatic candidates file not found: ${GENOTYPING_CANDIDATES}"
  exit 1
fi

# resolve the full path to the input files
TUMOUR_BAM=$(realpath "${TUMOUR_BAM}")
REFERENCE=$(realpath "${REFERENCE}")
GERMLINE_RESOURCE=$(realpath "${GERMLINE_RESOURCE}")
PANEL_OF_NORMALS=$(realpath "${PANEL_OF_NORMALS}")
GENOTYPING_CANDIDATES=$(realpath "${GENOTYPING_CANDIDATES}")

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_PATH}"
if [ ! -d "${OUTPUT_PATH}" ]; then
  echo "Error: Unable to create output directory: ${OUTPUT_PATH}"
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

# Create an output filename from the BAM file name
BAM_NAME=$(basename "${TUMOUR_BAM}")
# Remove any .bam or .cram extension
BAM_NAME=${BAM_NAME%.bam}
BAM_NAME=${BAM_NAME%.cram}

singularity exec \
  -B '/nfs:/nfs' \
  -B '/lustre:/lustre' \
  --no-home \
  --env JAVA_OPTIONS='-Xmx16g' \
  -e ${PROJ_DIR}/container/gatk_latest.sif \
  gatk Mutect2 \
      --reference ${REFERENCE} \
      --input ${TUMOUR_BAM} \
      --germline-resource ${GERMLINE_RESOURCE} \
      --panel-of-normals ${PANEL_OF_NORMALS} \
      --alleles ${GENOTYPING_CANDIDATES} \
      --f1r2-tar-gz ${OUTPUT_PATH}/${BAM_NAME}.mutect2.f1r2.tar.gz \
      --output ${OUTPUT_PATH}/${BAM_NAME}.mutect2.unfiltered.vcf.gz \
      $INTERVALS \
      --interval-padding 150
