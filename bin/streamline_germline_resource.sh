#!/bin/bash

VCF=$1
REF=$2

if [ -z "$VCF" ] || [ -z "$REF" ]; then
  echo "Usage: $0 variants.vcf.gz ref.fa"
  exit 1
fi

if [ ! -f "$VCF" ]; then
  echo "$VCF file not found"
  exit 1
fi

if [ ! -f "$REF" ]; then
  echo "$REF file not found"
  exit 1
fi

TMP_BODY="$(mktemp)"
TMP_HEADER="$(mktemp)"
echo "${TMP_BODY} ${TMP_HEADER}"

cleanup() {
  rm -f "$TMP_BODY"
  rm -f "$TMP_HEADER"
}

trap 'exit_code=$?; cleanup; exit $exit_code' EXIT

bcftools view -h "$VCF" \
  | sed 's/\tFORMAT.*//' \
  | grep -v '^##FORMAT=' \
  > "$TMP_HEADER"

bcftools norm --multiallelics -both \
  --fasta-ref "$REF" \
  "$VCF" \
  | bcftools view -e 'ALT="*"' \
  | bcftools view -m2 -M2 -G \
    -i 'INFO/AF>=0.01 && INFO/AF!="."' \
  | bcftools query -f \
    '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\tAF=%INFO/AF\n' \
  > "$TMP_BODY"

cat "$TMP_HEADER" "$TMP_BODY" | bcftools view -W=tbi -Oz -o "${VCF%.vcf.gz}.biallelic.vcf.gz"
