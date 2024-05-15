#!/usr/bin/env python3
"""
Written by Pedro Feijao.
Canada's Michael Smith Genome Sciences Centre.
"""

__author__ = "Pedro feijao"
__version__ = "0.1.0"
import argparse
import sys
import vcfpy
from typing import Dict, List, Optional
import pandas as pd
import re
from enum import Enum


class VariantStatus(str, Enum):
    """Variant status after parsing, with respect to parent of origin."""

    NotFound = "Not Found"
    UnmatchedVariant = "Unmatched variant"
    UnphasedVariant = "Unphased variant"
    MaternalVariant = "Maternal variant"
    PaternalVariant = "Paternal variant"
    ComplexVariant = "Complex variant"

def single_parent_alt_status(vcf_alt, annotation_alt, ps_annotation):
    if vcf_alt == annotation_alt:
        if ps_annotation == "Mat|Pat":
            return VariantStatus.PaternalVariant
        elif ps_annotation == "Pat|Mat":
            return VariantStatus.MaternalVariant
    else:
        return VariantStatus.UnmatchedVariant
    
def both_parents_alt_status(alternate_alleles, annotation_alt, ps_annotation):
    assert len(alternate_alleles) == 2
    if alternate_alleles[0] == annotation_alt:
        if ps_annotation == "Ref_Pat|Mat":
            return VariantStatus.PaternalVariant
        else:
            return VariantStatus.MaternalVariant
    elif alternate_alleles[1] == annotation_alt:
        if ps_annotation == "Ref_Mat|Pat":
            return VariantStatus.PaternalVariant
        else:
            return VariantStatus.MaternalVariant
    else:
        return VariantStatus.UnmatchedVariant


def parse_parents(alternate_alleles, annotation_alt, call):
    if call['GT'] in ["0|1", "1|0"]:
        return single_parent_alt_status(alternate_alleles[0], annotation_alt, call['PS'])
    elif call['GT'] == "1|2":
        return both_parents_alt_status(alternate_alleles, annotation_alt, call['PS'])
    else:
        return VariantStatus.UnmatchedVariant


def parse_variant_status(vcf_record, alt, variant_type):
    
    # Parse varian status, depending on case:
    call = vcf_record.calls[0].data
    # Check if it is a insertion:
    if vcf_record.ALT[0].type != variant_type:
        return VariantStatus.UnmatchedVariant

    # Phasing info:
    if "PS" not in call.keys():
        return VariantStatus.UnphasedVariant

    # Phased, try to find parent:
    if variant_type == "SNV":
        alternate_alleles = [alt.value for alt in vcf_record.ALT]
    elif variant_type == "INS":
        alternate_alleles = [alt.value[1] for alt in vcf_record.ALT]
    elif variant_type == "DEL":
        alternate_alleles = [vcf_record.REF[1:]]
    else:
        raise NotImplementedError(f"Variant type {variant_type} is not supported.")
    return parse_parents(alternate_alleles, alt, call)


def fetch_record(reader, chrom, start_pos, end_pos, single=True):
    records = [r for r in reader.fetch(chrom, start_pos, end_pos)]
    if len(records) == 0:
        return None
    if single:
        assert len(records) == 1
        return records[0]
    return records


def main(raw_arguments=None):
    """Do the thing.

    Args:
        raw_arguments (List[str], optional): Argument list, defaults to sys.argv is not given.
    """
    if raw_arguments is None:
        raw_arguments = sys.argv[1:]

    args = _parse_arguments(raw_arguments)

    # Read all annotations:
    annotations = pd.concat(
        [
            pd.read_csv(annotation_file, sep="\t")
            for annotation_file in args.annotated_variants
        ]
    ).sort_values(['coordinates(gDNA/cDNA/protein)'])

    # List to store status of each annotation:
    variant_status_list = []

    # VCF parser/writer:
    reader = vcfpy.Reader.from_path(args.vcf)
    writer = vcfpy.Writer.from_path(args.output_vcf, reader.header)
    # Iterate over all annotations: 
    for _, row in annotations.iterrows():
        # Get full gDNA annotation:
        chromosome_gdna_position = row["coordinates(gDNA/cDNA/protein)"].split("/")[0]
        chromosome, gdna_position = chromosome_gdna_position.split(":")
        # get left_aligned string if present, works better for indels:
        left_align_re = re.match(".+left_align_gDNA=(.*?);", row['info'])
        if left_align_re:
            gdna_position = left_align_re.group(1)

        # SPECIAL CASE: delins
        # Most likely >1 record ; hard to parse. For now just output the VCF records and 
        # mark output as "Complex Variant", for manual inspection.
        delins_re = re.match("[gG]\.(\d+)(?:_(\d+))?([GCTAgcta]+)?delins([GCTAgcta]+)", gdna_position)
        if delins_re:
            start_pos, end_pos, deletion, insertion = delins_re.groups()
            start_pos = int(start_pos) - 2
            end_pos = int(end_pos) if end_pos else start_pos + 1

            # INDELS might have more than one record:
            records = fetch_record(reader, chromosome, start_pos, end_pos, single=False)
            if records:
                variant_status = VariantStatus.ComplexVariant
                for vcf_record in records:
                    writer.write_record(vcf_record)
            else:
                variant_status = VariantStatus.NotFound
            variant_status_list.append(variant_status)
            continue

        # Match SNV. Example: "g.20625337T>G"
        snv_re = re.match("[gG]\.(\d+)([GCTAgcta])?>([GCTAgcta])", gdna_position)
        if snv_re:  # SNV found:
            pos, ref, alt = snv_re.groups()
            start_pos = int(pos) - 1
            end_pos = start_pos + 1
            variant_type = 'SNV'

        # Match insertion OR deletion
        ins_del_re = re.match("[gG]\.(\d+)(?:_(\d+))?(ins|del)([GCTAgcta]+)?", gdna_position)
        if ins_del_re:
            start_pos, end_pos, variant_type, alt = ins_del_re.groups()
            start_pos = int(start_pos) - 1
            end_pos = int(end_pos) if end_pos else start_pos + 1
            variant_type = variant_type.upper()  # INS or DEL

        # Fetch record from VCF, by gdna position
        vcf_record = fetch_record(reader, chromosome, start_pos, end_pos)
        
        # If a VCF record was obtained, try to parse the variant to match the annotation;
        # Always write the variant to output VCF.           
        if vcf_record:
            variant_status = parse_variant_status(vcf_record, alt, variant_type)
            writer.write_record(vcf_record)
        else:    
            variant_status = VariantStatus.NotFound
        variant_status_list.append(variant_status)
        

    # Add all status to annotations and write:
    annotations["PofoStatus"] = [s.value for s in variant_status_list]
    annotations.drop(["region", "info"], axis="columns", inplace=True)
    annotations.to_csv(args.output_status, sep="\t", index=False)


def _parse_arguments(raw_arguments: List[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        raw_arguments (List[str]): list of arguments

    Returns:
        argparse.Namespace: namespace object with parsed arguments.
    """
    parser = argparse.ArgumentParser()

    # Input(s) for building the samplesheet:
    parser.add_argument(
        "-v",
        "--vcf",
        help="PoFo VCF file.",
    )
    parser.add_argument(
        "-a",
        "--annotated_variants",
        nargs="+",
        help="Variants of interest, annotated by TransVar.",
    )

    parser.add_argument(
        "-o", "--output_vcf", required=True, help="VCF file with the filtered variants."
    )
    parser.add_argument(
        "-s",
        "--output_status",
        required=True,
        help="Status file with results for each variant.",
    )

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )
    parsed_arguments = parser.parse_args(raw_arguments)
    return parsed_arguments


if __name__ == "__main__":
    main()
