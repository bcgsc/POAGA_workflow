#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {strandseq_concat_fastq} from './poaga_processes.nf'
include {strandseq_cutadapt} from './poaga_processes.nf'
include {strandseq_bowtie} from './poaga_processes.nf'
include {strandseq_tidy_reads} from './poaga_processes.nf'
include {strandseq_mark_duplicates} from './poaga_processes.nf'
include {strand_pass_fail_cells} from './poaga_processes.nf'
include {ashleys_qc} from './poaga_processes.nf'
include {strandseq_phase} from './poaga_processes.nf'
include {patmat} from './poaga_processes.nf'

/*
* pipeline input parameters
*/
		
params.input_csv = "/projects/rcorbettprj2/COVID19/ARTIC_V3"
params.output_folder
    
log.info """\
Run the POAGA pipeline using strand-seq and nanopore data
===================================
output_folder	    : ${params.output_folder}
samples_file        : ${params.input_csv}
"""
.stripIndent()

/*
* Main workflow
*/

//parses CSV file and runs the run_artic
workflow {
    //parse the json or csv input (TBD)
    input = Channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map{ row-> tuple(row.patient,row.strand_seq_fastq,row.clair3_vcf,row.bam) }

    strandseq(input)

    //ashley's qc
    ashleys_qc(input.join(strandseq.out.pass_strandseq_bam_folder))

    //strandseq_phase.R
    strandseq_phase(input.join(ashleys_qc.out.ashleys_pass_folder))

    //patmat
    patmat(input.join(strandseq_phase.out)

    //put patmat phased resuts into the bamå
}



workflow strandseq {
    take:
        input

    main:
        //concat_fastq
        strandseq_concat_fastq(input)
       
        //will process each fastq pair as its own job
        //strandseq_concat_fastq.out.transpose().view()
        strandseq_cutadapt(strandseq_concat_fastq.out.transpose())

        //bowtie2 - outputs sam
        strandseq_bowtie(strandseq_cutadapt.out, params.bowtie_index)

        //tidybam
        strandseq_tidy_reads(strandseq_bowtie.out)
        
        //remove bams with very few reads
        strandseq_mark_duplicates(strandseq_tidy_reads.out)

        //split bams into pass/fail sets
        strand_pass_fail_cells(strandseq_mark_duplicates.out.groupTuple())

    emit:
        passed_strand_seq_list = strand_pass_fail_cells.out.passed_strand_seq_list
        failed_strand_seq_list = strand_pass_fail_cells.out.failed_strand_seq_list
        pass_strandseq_bam_folder = strand_pass_fail_cells.out.pass_strandseq_bams
        fail_strandseq_bam_folder = strand_pass_fail_cells.out.fail_strandseq_bams

}