#!/usr/bin/env python3
import pysam
from pysam.libcbcf import VariantHeader


def fix_header(vcf_header: VariantHeader) -> VariantHeader:
    samplenames = sorted(vcf_header.samples)
    new_header = pysam.VariantHeader()

    added = False
    for rec in vcf_header.records:
        if rec.key != "tumor_sample":
            new_header.add_record(rec)
        else:
            if not added:
                for sample in samplenames:
                    new_header.add_meta("tumor_sample", value=sample)
                added = True

    for sample in samplenames:
        new_header.add_sample(sample)

    return new_header


def write_fixed_vcf(input_vcf: int, output_vcf: str) -> None:
    vcf_in = pysam.VariantFile(input_vcf, "r")
    new_header = fix_header(vcf_in.header)
    vcf_out = pysam.VariantFile(output_vcf, "w", header=new_header)

    for rec in vcf_in:
        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Fix the tumor_sample header in a somatic candidates VCF."
    )
    parser.add_argument(
        "input_vcf",
        type=str,
        help="Input VCF file with potentially incorrect tumor_sample header.",
    )
    parser.add_argument(
        "output_vcf",
        type=str,
        help="Output VCF file with fixed tumor_sample header.",
    )

    args = parser.parse_args()
    write_fixed_vcf(args.input_vcf, args.output_vcf)
