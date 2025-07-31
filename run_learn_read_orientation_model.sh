#!/bin/bash
PROJ_DIR="/lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2"

F1R2=${1}

if [[ -z "${F1R2}" ]]; then
  echo "Usage: $0 <F1R2 file>"
  exit 1
fi

if [ ! -f "${F1R2}" ]; then
  echo "Error: F1R2 file not found: ${F1R2}"
  exit 1
fi

F1R2=$(realpath "${F1R2}")

# Create an output filename from the BAM file name
OUT_NAME=$(basename "${F1R2}")
# Remove any .bam or .cram extension
OUT_NAME=${OUT_NAME%.bam}
OUT_NAME=${OUT_NAME%.cram}

singularity exec \
  -B '/nfs:/nfs' \
  -B '/lustre:/lustre' \
  --no-home \
  --env JAVA_OPTIONS='-Xmx16g' \
  -e ${PROJ_DIR}/container/gatk_latest.sif \
  gatk LearnReadOrientationModel \
    -I "${F1R2}" \
    -O ${OUT_NAME}.model.tar.gz \
