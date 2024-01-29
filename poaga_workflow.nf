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
include {haplotag} from './poaga_processes.nf'


// HELP message
if ( params.help ) {
    help = """
             |poaga_processes.nf: Parent of Origin (POAGA) pipeline using strand-seq and nanopore data.
             |
             |Usage: nextflow run poaga_processes.nf --input_csv <CSV> --output_folder <DIR>
             |
             |Required arguments:
             |  --input_csv             Location of the input CSV file. Each row corresponds to a sample to be processed.
             |                          The CSV file must have a header with 5 columns:
             |     - patient            Sample label.
             |     - strand_seq_fastq   Folder with the Strand-seq FastQ files for the sample.
             |     - clair3_vcf         Clair3 VCF output for the sample's Nanopore data.
             |     - bam                Nanopore BAM file.
             |     - bai                Index file for the Nanopore BAM file.
             |
             |  --output_folder         Output folder to store POAGA results.
             |
             |Optional arguments:
             |  --min_reads_threshold   Minimum number of high quality reads to pass a StrandSeq library file.
             |                [default: ${params.min_reads_threshold}]
             |  --chromosomes           Comma-separated list of chromosomes to process.
             |                [default: ${params.chromosomes}]
             |  --min_ashleys_qc_prob   Minimum probability score for a Strand-seq library to pass the ASHLEYS QC step.
             |                [default: ${params.min_ashleys_qc_prob}]
             |  -w                      The NextFlow work directory. 
             |                [default: ${workDir}]
             |  -profile                Nextflow config profile(s) with pre-defined parameters (comma-separated profiles to use more than one). 
             |                Available profiles include 'hg38_no_alt' (hg38 reference genome) and 'test' (run minimal dataset test)
             |
             |Useful Nextflow arguments:
             |  -resume                Does not execute steps that have already been executed with the same parameters. Useful for resuming
             |                         failed runs or re-runs with some parameter changes.
             |  -with-report           Creates an HTML report with many nextflow execution metrics.
             |  -with-trace            Creates a tab-separated file with time, cpu and memory stats of each process.
             |  -with-timeline         Creates an HTML execution timeline of the pipeline
             |""".stripMargin()

    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

if (!params.input_csv) {
    error("Error: Please specify the input CSV file with the --input_csv parameter.")
}
if (!params.output_folder) {
    error("Error: Please specify the output location with the --output_folder parameter.")
}


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
        .map{ row-> tuple(row.patient,row.strand_seq_fastq,row.clair3_vcf,row.bam,row.bai) }

    //implemented below. Does the basic aligning, tidying and read count qc for strand-seq
    strandseq(input)

    //ashley's qc - scores the cells, returns folders with the cells that pass and fail
    ashleys_qc(input.join(strandseq.out.pass_strandseq_bam_folder), params.min_ashleys_qc_prob)

    //strandseq_phase.R
    strandseq_phase(input.join(ashleys_qc.out.ashleys_pass_folder), params.chromosomes)

    //patmat
    patmat(input.join(strandseq_phase.out.vcf), params.reference_genome)

    //put patmat phased resuts into the bam
    haplotag(input.join(patmat.out.pofo_assignment_vcf), params.reference_genome, params.reference_genome_idx)


}

/*
Takes care of the preliminary strand-seq analysis including:
-making a single pair of fastq files per cell
-align with bowtie2
-sort, dupmark, etc.
-Filter out cells with low fragment counts

returns:
-folders for cells that pass and fail
-lists of cells that pass and fail
*/
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
        //strandseq_mark_duplicates.out.groupTuple().view()
        strand_pass_fail_cells(strandseq_mark_duplicates.out.groupTuple(), params.min_reads_threshold)

    emit:
        passed_strand_seq_list = strand_pass_fail_cells.out.passed_strand_seq_list
        failed_strand_seq_list = strand_pass_fail_cells.out.failed_strand_seq_list
        pass_strandseq_bam_folder = strand_pass_fail_cells.out.pass_strandseq_bams
        fail_strandseq_bam_folder = strand_pass_fail_cells.out.fail_strandseq_bams

}