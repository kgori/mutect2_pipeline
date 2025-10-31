#!/usr/bin/env python3
import pysam
from pysam.libcfaidx import FastaFile
from dataclasses import dataclass
import sys

import signal

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


@dataclass(eq=True, frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str


def left_align_and_normalise(variant: Variant, reference: FastaFile) -> Variant:
    """
    Normalisation algorithm used by vt:
    see https://academic.oup.com/bioinformatics/article/31/13/2202/196142
    """
    ref = variant.ref
    alt = variant.alt
    pos = variant.pos

    changed = True
    while changed:  # loop until no more changes
        changed = False
        if ref[-1] == alt[-1]:  # truncate last base if identical
            ref = ref[:-1]
            alt = alt[:-1]
            changed = True
        if len(alt) == 0 or len(ref) == 0:  # shift pos left by 1
            pos -= 1
            ref_base = reference.fetch(variant.chrom, pos - 1, pos)
            ref = ref_base + ref
            alt = ref_base + alt
            changed = True

    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    return Variant(variant.chrom, pos, ref, alt)


def get_pon_variants(pon_vcf: str, reference_fasta: str) -> frozenset[Variant]:
    print("Loading Panel of Normals variants", file=sys.stderr)
    variants = set()
    vcf = pysam.VariantFile(pon_vcf, "r")
    ref = pysam.FastaFile(reference_fasta)
    for rec in vcf:
        for alt in rec.alts:
            variant = Variant(chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=alt)
            variant = left_align_and_normalise(variant, ref)
            if alt != "*":
                variants.add(variant)
    return frozenset(variants)


def check_vcf(vcf, outfile, pon_set: set[Variant]):
    print("Checking VCF for Panel of Normals annotations", file=sys.stderr)
    write_mode = 'wb'
    if vcf == "-":
        print("Reading VCF from standard input", file=sys.stderr)
    if outfile == "-":
        print("Writing output VCF to standard output", file=sys.stderr)
        write_mode = 'w'
    edits = 0
    records = 0
    vcf_in = pysam.VariantFile(vcf, "r")
    vcf_out = pysam.VariantFile(outfile, write_mode, header=vcf_in.header)
    for rec in vcf_in:
        records += 1
        if records % 100000 == 0:
            print(f"Progress - {rec.chrom}:{rec.pos}", file=sys.stderr)
        assert len(rec.alts) == 1, f"Record {rec.chrom}:{rec.pos} has multiple alts"
        pon_in_filter = "panel_of_normals" in rec.filter
        pon_in_info = "PON" in rec.info
        if pon_in_filter and pon_in_info:
            variant = Variant(
                chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=rec.alts[0]
            )
            if variant not in pon_set:
                del rec.filter["panel_of_normals"]
                if len(rec.filter) == 0:
                    rec.filter.add("PASS")
                del rec.info["PON"]
                edits += 1

        elif pon_in_filter and not pon_in_info:
            raise ValueError(
                f"Record {rec.chrom}:{rec.pos} has panel_of_normals in FILTER but not in INFO PON"
            )
        elif pon_in_info and not pon_in_filter:
            raise ValueError(
                f"Record {rec.chrom}:{rec.pos} has PON in INFO but not panel_of_normals in FILTER"
            )
        vcf_out.write(rec)
    print(f"Total edits made: {edits}", file=sys.stderr)
    print(f"Output written to {outfile}", file=sys.stderr)
    return edits


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Finalise Panel of Normals annotations in a filtered Mutect2 VCF"
    )
    parser.add_argument("--pon-vcf", required=True, help="Panel of Normals VCF file")
    parser.add_argument(
        "--reference-fasta",
        required=True,
        help="Reference FASTA file (requires fai index)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite output VCF if it exists"
    )
    parser.add_argument("input_vcf", help="Input VCF file to check", default="-")
    parser.add_argument("-o", "--output_vcf", help="Output VCF file", default="-")

    return parser.parse_args()


def check_args(args):
    import os

    if not os.path.exists(args.reference_fasta):
        raise FileNotFoundError(
            f"Reference FASTA file {args.reference_fasta} not found"
        )
    if args.input_vcf != "-" and not os.path.exists(args.input_vcf):
        raise FileNotFoundError(f"Input VCF file {args.input_vcf} not found")
    if not os.path.exists(args.pon_vcf):
        raise FileNotFoundError(f"Panel of Normals VCF file {args.pon_vcf} not found")
    if os.path.exists(args.output_vcf) and not args.overwrite:
        raise FileExistsError(
            f"Output VCF file {args.output_vcf} already exists and --overwrite option is not set"
        )

    if args.output_vcf != "-":
        try:
            open(args.output_vcf, "wb").close()
            os.remove(args.output_vcf)
        except Exception as e:
            raise PermissionError(
                f"Cannot write to output VCF file {args.output_vcf}: {e}"
            )


if __name__ == "__main__":
    args = parse_args()

    try:
        check_args(args)
    except FileNotFoundError as e:
        print(f"FATAL ERROR: FILE NOT FOUND: {e}", file=sys.stderr)
        sys.exit(1)
    except FileExistsError as e:
        print(f"FATAL ERROR: NO OVERWRITE: {e}", file=sys.stderr)
        sys.exit(1)
    except PermissionError as e:
        print(f"FATAL ERROR: NOT WRITEABLE: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"FATAL ERROR: UNKNOWN: {e}", file=sys.stderr)
        sys.exit(1)

    pon_variants = get_pon_variants(args.pon_vcf, args.reference_fasta)
    check_vcf(args.input_vcf, args.output_vcf, pon_variants)
    print("Finished processing VCF.", file=sys.stderr)
