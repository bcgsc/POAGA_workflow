# POAGA Parent of Origina scripts

Using strand-seq and nanopore long reads this project leverages known imprinted regions
to identify which variants in a child were inherited from each parent

## Example command

nextflow run poaga_workflow.nf --input_csv /home/rcorbett/RCORBETT/poaga_workflow/example_samples.csv --output_folder aaa_255_test -profile hg38_no_alt  -process.executor=local -resume

### Assumptions made, future changes to make

1. Read names have the format \*_*R{1,2}*.fastq.gz
2. Fragment count minimum for usable strand-seq cells is 5000
3. Currently, all strand-seq fastqs in a folder are used.    Should be parsed from a file.
4. Needs patmat phasing put into the bam
5. 