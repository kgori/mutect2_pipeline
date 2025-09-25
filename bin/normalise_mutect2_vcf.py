#!/usr/bin/env python3
import pysam
from pysam.libcbcf import VariantRecordFilter
from dataclasses import dataclass

import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def number_g_indexes(allele_index: int):
    """ Return the indexes required to extract the 'allele_index'th
    allele from a Number=G field.
    """
    i = 0
    j = (allele_index * (allele_index + 1)) // 2
    k = j + allele_index
    return i, j, k 

def split_number_g_field(array: tuple[int | float, ...], allele_index: int, total_alleles: int):
    """ Number=G fields have entries for all possible genotypes.
    The order is specified in the VCF spec (>=4.3), under the
    "Genotype Ordering" heading. They go 0/0, 0/1, 1/1, 0/2, 1/2, 2/2,
    0/3, 1/3, 2/3, 3/3, etc.
    allele_index is the index of the allele being extracted (0 is the
    reference allele).
    Total alleles is the total number of alleles (reference + alternates).
    """
    assert len(array) == (total_alleles * (total_alleles + 1)) // 2
    assert 0 < allele_index < total_alleles
    i, j, k = number_g_indexes(allele_index)
    return array[i], array[j], array[k]
    
def split_number_r_field(array: tuple[int | float, ...], allele_index: int):
    """ Number=R fields have entries for all alleles, including the
    reference allele.
    allele_index is the index of the allele being extracted (0 is the
    reference allele).
    """
    assert 0 < allele_index < len(array) 
    return array[0], array[allele_index]

def split_number_a_field(array: tuple[int | float, ...], allele_index: int):
    """ Number=A fields have entries for each alternate allele.
    allele_index is the index of the allele being extracted (0 is the
    reference allele).
    """
    assert 0 < allele_index <= len(array) 
    return array[allele_index - 1]

@dataclass
class ASFilterStatus:
    INFO: str
    FILTERKEYS: tuple[str, ...]

def split_as_filter_status(array: tuple[str, ...], allele_index: int, site_flags: VariantRecordFilter):
    """ AS_FilterStatus is a GATK-specific field that is not
    standards-compliant. It contains a comma-separated list of
    filter flags for each allele, and each allele is separated
    by '|'. This gets decomposed "incorrectly" by pysam/bcftools,
    so it needs to be reconstructed and resplit.
    The special flag "SITE" is used to indicate that the allele
    should get the same filter status as the site itself.
    Returns: A tuple of (allele_info_flags, site_flags), where
    allele_info_flags is a comma-separated string to use in the INFO
    field, and site_flags is the tuple of flags found in FILTER,
    as can be obtained from VariantRecordFilter.
    """
    alleles = ','.join(array).split('|')
    assert 0 < allele_index <= len(alleles)
    allele = alleles[allele_index - 1]
    if allele == "SITE":
        return ASFilterStatus(INFO="SITE", FILTERKEYS=site_flags.keys())
    else:
        as_filter_status = ASFilterStatus(INFO=allele, FILTERKEYS=allele.split(','))
        return as_filter_status
    
def split_as_sb_table(array: tuple[str, ...], allele_index: int):
    """ AS_SB_TABLE is a GATK-specific field that is not
    standards-compliant. It contains a table of strand-bias
    metrics for each allele, with each row separated by '|'
    and each column separated by ','.
    The first row is for the reference allele, and subsequent
    rows are for each alternate allele.
    Returns: A tuple of (ref_metrics, alt_metrics), where each
    is a comma-separated string of metrics.
    """
    alleles = ','.join(array).split('|')
    assert 0 < allele_index < len(alleles)
    return '|'.join((alleles[0], alleles[allele_index]))

def split_multiallelic_record(rec, bcf_out, info_fields, format_fields):
    if len(rec.alleles) <= 2:
        return [rec]
    
    records = []
    for ai in range(1, len(rec.alleles)):
        new_rec = bcf_out.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=(rec.alleles[0], rec.alleles[ai]),
            id=rec.id,
            qual=rec.qual
        )
        new_filter_keys = set()
        for k, v in rec.info.items():
            if k not in info_fields:
                continue
            if k == "AS_FilterStatus":
                as_filter_status = split_as_filter_status(v, ai, rec.filter)
                new_rec.info[k] = as_filter_status.INFO
                new_filter_keys.update(as_filter_status.FILTERKEYS)
            elif k == "AS_SB_TABLE":
                new_rec.info[k] = split_as_sb_table(v, ai)
            elif info_fields[k] == "R":
                new_rec.info[k] = split_number_r_field(v, ai)
            elif info_fields[k] == "A":
                new_rec.info[k] = split_number_a_field(v, ai)
            elif info_fields[k] == "G":
                new_rec.info[k] = split_number_g_field(v, ai, len(rec.alleles))
            else:
                new_rec.info[k] = v
            if k == "PON" and v:
                new_filter_keys.add("panel_of_normals")

        for fk in sorted(new_filter_keys):
            new_rec.filter.add(fk)

        for sample_name, sample_data in rec.samples.items():
            new_sample_data = new_rec.samples[sample_name]
            alleles = sample_data.alleles
            new_sample_data.alleles = tuple(sample_data.alleles[ai]
                                            if i == ai
                                            else sample_data.alleles[0]
                                            for i in range(len(sample_data.alleles)))
            for k, v in sample_data.items():
                if k not in format_fields:
                    continue
                if k == 'GT':
                    continue
                if format_fields[k] == "R":
                    new_sample_data[k] = split_number_r_field(v, ai)
                elif format_fields[k] == "A":
                    new_sample_data[k] = split_number_a_field(v, ai)
                elif format_fields[k] == "G":
                    new_sample_data[k] = split_number_g_field(v, ai, len(rec.alleles))
                else:
                    new_sample_data[k] = v
            
        records.append(new_rec)
    return records

def get_fields_from_header(bcf_in):
    """ Get the INFO and FORMAT fields from the header and return
    what kind of "Number" format they are (R, A, G, <scalar>, etc)
    """
    info_fields = {}
    format_fields = {}
    for key, metadata in bcf_in.header.info.items():
        info_fields[key] = metadata.number
    for key, metadata in bcf_in.header.formats.items():
        format_fields[key] = metadata.number
    return info_fields, format_fields

def split_multiallelics(input_vcf: str):
    bcf_in = pysam.VariantFile(input_vcf, "r")
    header = bcf_in.header.copy()
    header.add_meta("split_multiallelic", value="normalise_mutect2_vcf.py")
    bcf_out = pysam.VariantFile("-", "w", header=header)
    info_fields, format_fields = get_fields_from_header(bcf_in)
    for rec in bcf_in:
        split_recs = split_multiallelic_record(rec, bcf_out, info_fields, format_fields)
        for split_rec in split_recs:
            bcf_out.write(split_rec)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Split multiallelic records in a Mutect2 VCF")
    parser.add_argument("input_vcf", help="Input VCF file, or '-' for stdin")
    args = parser.parse_args()
    split_multiallelics(args.input_vcf)

if __name__ == "__main__":
    main()
