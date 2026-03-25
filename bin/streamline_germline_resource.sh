#!/bin/bash

VCF=$1

if [ -z "$VCF" ]; then
  echo "No VCF file provided"
  exit 1
fi

if [ ! -f "$VCF" ]; then
  echo "$VCF file not found"
  exit 1
fi

bcftools view -h "$VCF" \
  | sed 's/\tFORMAT.*//' \
  | grep -v '^##FORMAT=' \
  > tmp_header

bcftools view -m2 -M2 -G -v snps \
  -i 'INFO/AF>0.01 && INFO/AF!="."' "$VCF" \
  | bcftools query -f \
  '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\tAF=%INFO/AF\n' \
  > tmp_body

cat tmp_header tmp_body | bcftools view -W=tbi -Oz -o "${VCF%.vcf.gz}.biallelic_snps.vcf.gz" \
  && rm tmp_{header,body}
