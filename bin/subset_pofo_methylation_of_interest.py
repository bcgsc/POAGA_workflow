#!/usr/bin/env python3
"""
Written by Pedro Feijao.
Canada's Michael Smith Genome Sciences Centre.
"""

__author__ = "Pedro feijao"
import argparse
import sys
from typing import Dict, List, Optional
import pandas as pd
import re
from enum import Enum

def main(raw_arguments=None):
    """Given a BED file and the PofO TSV methylation files, output subset versions with the BED chromosomes.

    Args:
        raw_arguments (List[str], optional): Argument list, defaults to sys.argv is not given.
    """
    if raw_arguments is None:
        raw_arguments = sys.argv[1:]

    args = _parse_arguments(raw_arguments)

    # Open bed file and get chromosomes:
    bed_file = pd.read_csv(args.bed_file, sep="\t", header=None)
    chromosomes = bed_file[0].unique()

    # open and subset Scores file:
    scores = pd.read_csv(args.scores_file, sep="\t")
    scores[scores.Chromosome.isin(chromosomes)].to_csv(args.scores_file_output, sep="\t", index=False)

    # open and subset DMRs file:
    dmr_status = pd.read_csv(args.dmr_status, sep="\t")
    dmr_status[dmr_status.Chromosome.isin(chromosomes)].to_csv(args.dmr_status_output, sep="\t", index=False)



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
        "-b",
        "--bed_file",
        help="BED file with subset chromosomes.",
    )
    parser.add_argument(
        "-s",
        "--scores_file",
        help="Methylation scores file.",
    )
    parser.add_argument(
        "-d",
        "--dmr_status",
        help="DMR status file.",
    )
    parser.add_argument(
        "-so",
        "--scores_file_output",
        help="Output for the Methylation scores subset file.",
    )
    parser.add_argument(
        "-do",
        "--dmr_status_output",
        help="Output for the DMR status subset file.",
    )
    parsed_arguments = parser.parse_args(raw_arguments)
    return parsed_arguments


if __name__ == "__main__":
    main()
