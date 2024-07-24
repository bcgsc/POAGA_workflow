#!/usr/bin/env bash
set -u
set -e
trap 'echo "Aborting due to errexit on line $LINENO. Exit code: $?" >&2' ERR
set -E
set -o pipefail

# This program's basename.
_ME="$(basename "${0}")"

# Print the program help information.
_print_help() {
  cat <<HEREDOC
Prepare VCF and BAM files for Haplotagging with whatshap.
Usage:
  ${_ME} ANNOTATION1 [ ANNOTATION2 ...]
  ${_ME} -h | --help

Options:
  -h --help  Show this screen.
HEREDOC
}


_parse_regions_of_interest() {
    FLANKING=$1
    shift
    IMPRINT=$1
    shift
    ANNOTATIONS=$@
    cut -f3 $ANNOTATIONS | grep -v gene | sort | uniq > genes.txt
    transvar canno --refseq -l genes.txt --longest > gene_annotations.txt
    cut -f3,5 gene_annotations.txt | grep -v coordinates | cut -f1 -d"/" | sed 's/:/\t/' | sed 's/g\.//' | sed 's/_/\t/' | awk '{print $2"\t"$3"\t"$4"\t"$1}' > genes.bed
    # get chromosomes
    cut -f1 genes.bed | sort | uniq | sed 's/$/\t/' > chr.txt
    # get imprint regions from these chromosomes:
    grep -f chr.txt $IMPRINT | cut -f1-4 > imprint_regions.bed
    # combine:
    awk -v FLANK="$FLANKING" '{print $1"\t"$2-FLANK"\t"$3+FLANK"\t"$4}' genes.bed imprint_regions.bed
}

_main() {
  if [[ "${1:-}" =~ ^-h|--help$ ]]  # -h or --help
  then
    _print_help
  else
    _parse_regions_of_interest "$@"
  fi
}

# Call `_main` after everything has been defined.
_main "$@"

