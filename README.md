# POAGA Parent of Origina scripts

using strand-seq and nanopore long reads this project leverages known imprinted regions
to identify which variants came from which parent

### Assumptions made, future changes to make

1. Read names have the format \*_*R{1,2}*.fastq.gz
2. Fragment count minimum for usable strand-seq cells is 5000
3. Currently, all strand-seq fastqs in a folder are used.    Should be parsed from a file.
4. 