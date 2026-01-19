#!/usr/bin/env python3
"""
Written by Pedro Feijao.
Canada's Michael Smith Genome Sciences Centre.
"""

__author__ = "Pedro feijao"
__version__ = "0.1.0"
import argparse
import sys
from typing import Dict, List, Optional

import pandas as pd


def main(raw_arguments=None):
    """

    Args:
        raw_arguments (List[str], optional): Argument list, defaults to sys.argv is not given.
    """
    if raw_arguments is None:
        raw_arguments = sys.argv[1:]

    args = _parse_arguments(raw_arguments)

    # read chr states:
    states = pd.read_csv(args.chr_states_input, sep="\t").fillna("?")

    # generate summary counts dataframe
    summary = (
        states.set_index("library")
        .melt(var_name="Chromosome", value_name="State")
        .groupby(["Chromosome"])
        .agg(N=("State", "value_counts"))
    )

    # write
    summary.to_csv(args.output, sep="\t")


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
        "chr_states_input",
        help="Stranseq Chr states file.",
    )
    parser.add_argument(
        "output",
        help="Output file name.",
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
