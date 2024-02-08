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
  ${_ME} VCFFILE BAMFILE THREADS
  ${_ME} -h | --help

Options:
  -h --help  Show this screen.
HEREDOC
}


_prepare_files() {
    VCFFILE=$1
    BAMFILE=$2
    THREADS=$3

    # Version log
    samtools --version

    # Create a header for the tagged bam file
    samtools view -H $BAMFILE | grep "^@SQ" > sam_header.txt

    # Create a vcf for whatshap
    (grep "^#" $VCFFILE && awk -F'\t' '$1 !~ /^#/ {print $1,$2,$3,$4,$5,$6,$7,".","GT",$10}' OFS='\t' $VCFFILE | cut -d':' -f1) | bgzip > PofO_Assignment_ForWhatsHap.vcf.gz

    tabix -f -p vcf PofO_Assignment_ForWhatsHap.vcf.gz

    # Create bam file for whatshap 
    samtools view -@ $THREADS $BAMFILE | sed 's/\\tHP:i:2//g; s/\\tHP:i:1//g' | cat sam_header.txt - | samtools sort -@ $THREADS -o whatshap.bam
    samtools index -@ $THREADS whatshap.bam

    # remove temp files
    rm sam_header.txt
}

_main() {
  if [[ "${1:-}" =~ ^-h|--help$ || $# != 3 ]]  # Asked for help or passed wrong number of parameters
  then
    _print_help
  else
    _prepare_files "$@"
  fi
}

# Call `_main` after everything has been defined.
_main "$@"

