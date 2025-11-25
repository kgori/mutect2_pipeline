#!/bin/bash
Help() {
  echo -e "\033[1;4;37mRealign Reads Script\033[0m"
  echo ""
  echo -e " Usage: \033[1m$(basename $0)\033[0m [\033[36moptions\033[0m] \033[32m<input.bam> <reference.fa>\033[0m"
  echo
  echo " Realign reads in the input BAM file to the reference genome using BWA-MEM2 or BWA."
  echo
  echo -e " \033[36mOptions\033[0m:"
  echo "  -o, --output        <file>               Output CRAM file name. If not provided, defaults to <input>_realigned.cram."
  echo "  -O, --output-format <BAM|CRAM|CRAM3.1>   Output format (choose BAM or CRAM; default: CRAM, aka CRAMv3.0)."
  echo "  -t, --threads       <num>                Number of threads to use (default: 16)."
  echo "      --dry-run                            Print the commands without executing them."
  echo "      --write-index                        Write an index for the output file."
  echo "  -h, --help                               Display this help message."
}

# Defaults
OUTPUT_FORMAT="CRAM"
THREADS=16
DRY_RUN=0
WRITE_INDEX=0
OUT=""

# Getopt
ARGS=$(getopt -o o:O:t:h --long output:,output-format:,threads:,dry-run,write-index,help -n "$0" -- "$@")
if [ $? -ne 0 ]; then
    Help
    exit 1
fi
eval set -- "$ARGS"

while true; do
  case "$1" in
    -o|--output)
      OUT="$2"
      shift 2
      ;;
    -O|--output-format)
      OUTPUT_FORMAT="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --write-index)
      WRITE_INDEX=1
      shift
      ;;
    -h|--help)
      Help
      exit 0
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "Error: Invalid option '$1'" >&2
      Help
      exit 1
      ;;
  esac
done

# Positional arguments
BAM="${1:-}"
REF="${2:-}"

if [ -z "$BAM" ] || [ -z "$REF" ]; then
  echo -e "\033[1;31mError: Missing required arguments.\033[0m" >&2
  Help
  exit 1
fi

if [ ! -e "$BAM" ]; then
  echo -e "\033[1;31mError: Input BAM file '$BAM' not found!\033[0m" >&2
  exit 1
fi

if [ ! -e "$REF" ]; then
  echo -e "\033[1;31mError: Reference genome file '$REF' not found!\033[0m" >&2
  exit 1
fi

# Check which programs are available, bwa-mem2 or bwa?
HAVE_BWA2=$(command -v bwa-mem2)
HAVE_BWA=$(command -v bwa)

if [ -z "$HAVE_BWA2" ] && [ -z "$HAVE_BWA" ]; then
  echo -e "\033[1;31mError: Neither bwa-mem2 nor bwa is installed or in PATH\033[0m" >&2
  exit 1
fi

# Check which index files are available
echo "Checking for "${REF}.bwt.2bit.64" index files" >&2
HAVE_BWA2_INDEX=$([ -f "${REF}.bwt.2bit.64" ] && echo 1 || echo 0)
HAVE_BWA_INDEX=$([ -f "${REF}.bwt" ] && echo 1 || echo 0)

# Select bwa program based on what's available. Prefer bwa-mem2.
if [ -n "$HAVE_BWA2" ] && [ "$HAVE_BWA2_INDEX" -eq 1 ]; then
  BWA="bwa-mem2"
elif [ -n "$HAVE_BWA" ] && [ "$HAVE_BWA_INDEX" -eq 1 ]; then
  BWA="bwa"
else
  echo -e "\033[1;31mError: Index files are missing for available aligners.\033[0m" >&2
  if [ -n "$HAVE_BWA2" ]; then
    echo "Please run bwa-mem2 index on the reference genome." >&2
  else
    echo "Please run bwa index on the reference genome." >&2
  fi
  exit 1
fi

re='^[0-9]+$'
if ! [[ $THREADS =~ $re ]] ; then
  echo -e "\033[1;31mError: Number of threads is not a number\033[0m" >&2
  exit 1
fi

if [ $THREADS -lt 1 ]; then
  echo -e "\033[1;31mError: Number of threads must be at least 1\033[0m" >&2
  exit 1
fi

if [ -z "$OUT" ]; then
  if [ "$OUTPUT_FORMAT" == "BAM" ]; then
    OUT="${BAM%.*am}_realigned.bam"
  else
    OUT="${BAM%.*am}_realigned.cram"
  fi
fi

if [ "$OUTPUT_FORMAT" != "CRAM" ] && [ "$OUTPUT_FORMAT" != "BAM" ] && [ "$OUTPUT_FORMAT" != "CRAM3.1" ]; then
  echo -e "\033[1;31mError: Unsupported output format '$OUTPUT_FORMAT'. Use BAM, CRAM or CRAM3.1.\033[0m" >&2
  exit 1
fi

FMT="-O CRAM,version=3.0 --reference ${REF}"
if [ "$OUTPUT_FORMAT" == "BAM" ]; then
  FMT="-O BAM"
elif [ "$OUTPUT_FORMAT" == "CRAM3.1" ]; then
  FMT="-O CRAM,version=3.1 --reference ${REF}"
fi

ALIGN_THREADS=$THREADS
SORT_THREADS=$((THREADS / 2))
if [ $SORT_THREADS -lt 1 ]; then
  SORT_THREADS=1
fi

# All options are valid, so make a temporary directory
TMPDIR="$(mktemp -d -t realign_tmp_XXXXXX)"

cleanup() {
  rm -rf "$TMPDIR"
}

trap 'exit_code=$?; cleanup; exit $exit_code' EXIT

echo -e "\033[4;32m$(basename $0)\033[0m" >&2
echo -e "Using aligner:    \033[1;34m$BWA\033[0m" >&2
echo -e "Input BAM:        \033[1;34m$BAM\033[0m" >&2
echo -e "Reference genome: \033[1;34m$REF\033[0m" >&2
echo -e "Output file:      \033[1;34m$OUT\033[0m" >&2
echo -e "Output format:    \033[1;34m$OUTPUT_FORMAT\033[0m" >&2
echo -e "Total threads:    \033[1;34m$THREADS \033[34m(Align: $ALIGN_THREADS, Sort: $SORT_THREADS)\033[0m" >&2

CMD="samtools collate -Oun128 -T ${TMPDIR}/collate ${BAM} \
    | samtools fastq -OT RG,BC - \
    | ${BWA} mem -p -Y -K 100000000 -t ${ALIGN_THREADS} -CH <(samtools view -H ${BAM} | grep ^@RG) ${REF} - \
    | samtools sort -n -m 2G -@ ${SORT_THREADS} -T ${TMPDIR}/namesort -o - \
    | samtools fixmate -m - - \
    | samtools sort -m 2G -@ ${SORT_THREADS} -T ${TMPDIR}/possort -o - \
    | samtools markdup -@ ${SORT_THREADS} -T ${TMPDIR}/markdup ${FMT} - ${OUT}"

if [ $DRY_RUN -eq 1 ]; then
  echo ""
  echo -e "\033[1;31m**Dry run mode. The following command would be executed:**\033[0m" >&2
  echo "$CMD" >&2
  if [ "$WRITE_INDEX" -eq 1 ]; then
      echo "samtools index ${OUT}" >&2
  fi
  exit 0
else
  echo -e "\033[1;33mExecuting realignment command\n${CMD}\033[0m" >&2

  eval "$CMD"

  if [ $? -ne 0 ]; then
    echo -e "\033[1;31mError: Realignment command failed\033[0m" >&2
    exit 1
  fi

  # Write the index and check it was created
  if [ "$WRITE_INDEX" -eq 1 ]; then
    samtools index "${OUT}"
    if [ $? -ne 0 ]; then
      echo -e "\033[1;31mError: Failed to write index for output file '${OUT}'\033[0m" >&2
      exit 1
    fi

    index_exists() {
      local f=$1
      [[ -f ${f}.bai || -f ${f}.csi || -f ${f}.crai ]]
    }

    if ! index_exists "${OUT}"; then
      echo -e "\033[1;31mError: Index file for output '${OUT}' was not created!\033[0m" >&2
      exit 1
    fi
  fi

  # Check we produced some output
  [ -s "$OUT" ] || { echo -e "\033[1;31mError: Output file '$OUT' is empty!\033[0m" >&2; exit 1; }

  # Check output file integrity
  samtools quickcheck "$OUT" || { echo -e "\033[1;31mError: Output file '$OUT' is corrupted or incomplete!\033[0m" >&2; exit 1; }

  # Check the output contains any reads
  OUT_READS=$(samtools view -c -F 0x900 "$OUT")
  (( OUT_READS > 0 )) || {
    echo -e "\033[1;31mError: No reads found in output file '$OUT'!\033[0m" >&2
    exit 1
  }

  # Check we aligned the same number of reads as input
  IN_READS=$(samtools view -c -F 0x900 "$BAM")
  (( OUT_READS == IN_READS )) || {
    echo -e "\033[1;31mError: Mismatch in read counts! Input BAM has $IN_READS reads, but output has $OUT_READS reads.\033[0m" >&2
    exit 1
  }

  # All looks good
  echo -e "\033[1;32mRealignment completed successfully. ${OUT_READS} were realigned. Output written to '$OUT'.\033[0m" >&2
  exit 0
fi
