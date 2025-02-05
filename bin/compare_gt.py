#!/usr/bin/env python
import argparse
import sys
import typing

import pandas as pd


def _parse_arguments(raw_arguments: typing.Sequence[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        raw_arguments (Sequence[str]): list of arguments

    Returns:
        argparse.Namespace: namespace object with parsed arguments.
    """
    parser = argparse.ArgumentParser()

    # Input(s) for building the samplesheet:
    parser.add_argument(
        "-n",
        "--nanopore_vcf",
        required=True,
        help="Nanopore VCF file.",
    )
    parser.add_argument(
        "-s",
        "--strandseq_vcf",
        required=True,
        help="StrandSeq VCF from mpileup.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file with SNP concordance metrics.",
    )
    parsed_arguments = parser.parse_args(raw_arguments)
    return parsed_arguments


def main(raw_arguments=None):
    """Measure the SNP concordance between a Nanopore VCF file and a Strand-Seq VCF file.

    Args:
        raw_arguments (List[str], optional): Argument list, defaults to sys.argv is not given.
    """
    if raw_arguments is None:
        raw_arguments = sys.argv[1:]

    args = _parse_arguments(raw_arguments)

    nano_vcf = pd.read_csv(
        args.nanopore_vcf,
        sep="\t",
        comment="#",
        names=[
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "SAMPLE",
        ],
    ).query("REF.str.len() == 1 & ALT.str.len() == 1 & FILTER == 'PASS'")[
        ["CHROM", "POS", "REF", "ALT", "FILTER", "SAMPLE"]
    ]
    nano_vcf["SAMPLE"] = nano_vcf["SAMPLE"].apply(
        lambda x: list(set(x.split(":")[0].replace("/", "|").split("|")))
    )
    nano_vcf["GT"] = nano_vcf.apply(
        lambda df: [df["REF"] if x == "0" else df["ALT"] for x in df["SAMPLE"]], axis=1
    )

    sseq_vcf = pd.read_csv(args.strandseq_vcf, sep="\t", names=["CHROM", "POS", "NUC"])

    sseq_vcf["NUC"] = sseq_vcf["NUC"].apply(
        lambda x: [s for s in x.split(",") if s != "<*>"]
    )

    merged = pd.merge(nano_vcf, sseq_vcf, on=["CHROM", "POS"], how="inner")
    merged["concor"] = merged.apply(
        lambda df: all(item in df["GT"] for item in df["NUC"]), axis=1
    )

    df = pd.DataFrame(
        [
            (
                merged["concor"].shape[0],
                merged["concor"].value_counts().get(False),
                merged["concor"].value_counts(True).get(False, 0.0).round(5),
            )
        ],
        columns=["Total", "Discordant", "P"],
    )
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
