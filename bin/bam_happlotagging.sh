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
Haplotagging the nanopore bam file using the phased vcf file and whatshap.
Usage:
  ${_ME} VCFFILE BAMFILE REF_GENOME THREADS
  ${_ME} -h | --help

Options:
  -h --help  Show this screen.
HEREDOC
}


_haplotag_bam_file() {
    VCFFILE=$1
    BAMFILE=$2
    REF_GENOME=$3
    THREADS=$4

    # Version log
    whatshap --version
    samtools --version

    # Create a header for the tagged bam file
    samtools view -H $BAMFILE | grep "^@SQ" > sam_header.txt

    # Create a vcf for whatshap
    (grep "^#" $VCFFILE && awk -F'\t' '$1 !~ /^#/ {print $1,$2,$3,$4,$5,$6,$7,".","GT",$10}' OFS='\t' $VCFFILE | cut -d':' -f1) | bgzip > PofO_Assignment_ForWhatsHap.vcf.gz

    tabix -f -p vcf PofO_Assignment_ForWhatsHap.vcf.gz

    # Create bam file for whatshap 
    samtools view -@ $THREADS $BAMFILE | sed 's/\\tHP:i:2//g; s/\\tHP:i:1//g' | cat sam_header.txt - | samtools sort -@ $THREADS -o temp.bam
    samtools index -@ $THREADS temp.bam

    # haplotag with whatshap
    whatshap haplotag --ignore-read-groups --tag-supplementary --output-threads $THREADS --reference $REF_GENOME -o PofO_HaploTagged.bam PofO_Assignment_ForWhatsHap.vcf.gz temp.bam 

    # Sort and index bam file
    samtools sort -@ $THREADS PofO_HaploTagged.bam -o PofO_HaploTagged.sorted.bam
    samtools index -@ $THREADS PofO_HaploTagged.sorted.bam

    # remove temp files
    rm sam_header.txt temp.bam temp.bam.bai PofO_Assignment_ForWhatsHap.vcf.gz PofO_Assignment_ForWhatsHap.vcf.gz.tbi PofO_HaploTagged.bam
}

_main() {
  if [[ "${1:-}" =~ ^-h|--help$ || $# != 4 ]]  # Asked for help or passed wrong number of parameters
  then
    _print_help
  else
    _haplotag_bam_file "$@"
  fi
}

# Call `_main` after everything has been defined.
_main "$@"

