#!/bin/bash
PROJ_DIR="/lustre/scratch125/casm/staging/team267_murchison/kg8/20250704_Mutect2"

BAM=${1}
INTERVALS=${2}
echo $INTERVALS

if [ -z "${BAM}" ]; then
  echo "Usage: $0 <BAM file> [intervals file]"
  exit 1
fi

# resolve the full path to the BAM file
BAM=$(realpath "${BAM}")

if [ ! -f "${BAM}" ]; then
  echo "Error: BAM file not found: ${BAM}"
  exit 1
fi

if [ ! -z "$INTERVALS" ]; then
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
  -Ce ${PROJ_DIR}/container/gatk_latest.sif \
  gatk HaplotypeCaller \
      -R $PROJ_DIR/refs/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
      -I ${BAM} \
      -O $PROJ_DIR/${BAM_NAME}.g.vcf.gz \
      -ERC GVCF \
      $INTERVALS \
      --interval-padding 150 

