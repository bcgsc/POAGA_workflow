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
Runs ASHLEYS (https://github.com/friendsofstrandseq/ashleys-qc) to select high quality strandseq libraries.
Usage:
  ${_ME} BAM_FOLDER PASS_PROBABILITY
  ${_ME} -h | --help

Options:
  -h --help  Show this screen.
HEREDOC
}

_run_ashleys() {
  BAM_FOLDER=$1
  PASS_PROBABILITY=$2

  ashleys.py -j 48 \
  features \
  -f $BAM_FOLDER \
  -w 5000000 2000000 1000000 800000 600000 400000 200000 \
  -o ./features.tsv

  ashleys.py predict -p ./features.tsv -o ./quality.txt -m /opt/ashleys-qc/models/svc_default.pkl

  mkdir ashleys_pass
  mkdir ashleys_fail
  
  awk -v pval=$PASS_PROBABILITY 'NR>1 && $3>= pval { print $1 }' ./quality.txt | xargs -i cp $BAM_FOLDER/{} $BAM_FOLDER/{}.bai ashleys_pass
  awk -v pval=$PASS_PROBABILITY 'NR>1 && $3< pval { print $1 }' ./quality.txt | xargs -i cp $BAM_FOLDER/{} $BAM_FOLDER/{}.bai ashleys_fail

}
_main() {
  if [[ "${1:-}" =~ ^-h|--help$ || $# != 2 ]]  # Asked for help or passed wrong number of parameters
  then
    _print_help
  else
    _run_ashleys "$@"
  fi
}

# Call `_main` after everything has been defined.
_main "$@"

