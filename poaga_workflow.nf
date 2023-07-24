#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include {find_fastq_files} from './poaga_processes.nf'

/*
* pipeline input parameters
*/
		
params.input_csv = "/projects/rcorbettprj2/COVID19/ARTIC_V3"
    
log.info """\
QC whole genome data in BAM files
===================================
output_folder	    : ${params.output_folder}
samples_file        : ${params.samples_file}
"""
.stripIndent()


/*
* Main workflow
*/

//parses CSV file and runs the run_artic
workflow {
    //parse the json or csv input (TBD)
    input = Channel
        .fromPath(params.samples_file)
        .splitCsv(header:true)
        .map{ row-> tuple(row.flowcell, row.library, row.fast5, row.fastq, row.sequencing_summary) }

    //align, sort, dupmark strand-seq data

    //ashley's qc

    //strandseq_phase.R

    //patmat.py

    //put patmat phased resuts into the bamå
}