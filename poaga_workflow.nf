#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    ashleys_qc
    compress_haplotagged_bam
    log_pipeline_version
    patmat
    prepare_whatshap_files
    sort_haplotagged_bam
    strand_pass_fail_cells
    strandseq_bowtie
    strandseq_concat_fastq
    strandseq_cutadapt
    strandseq_mark_duplicates
    strandseq_phase
    strandseq_tidy_reads
    whatshap_haplotag

} from './poaga_processes.nf'

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

workflow {
    input = Channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map{ row-> tuple(row.patient,row.strand_seq_fastq,row.clair3_vcf,row.bam,row.bai) }

    // Log pipeline version:
    ch_versions = log_pipeline_version()

    //implemented below. Does the basic aligning, tidying and read count qc for strand-seq
    strandseq(input)
    ch_versions = ch_versions.mix(strandseq.out.ch_versions)

    //ashley's qc - scores the cells, returns folders with the cells that pass and fail
    ashleys_qc(input.join(strandseq.out.pass_strandseq_bam_folder), params.min_ashleys_qc_prob)
    ch_versions = ch_versions.mix(ashleys_qc.out.versions.take(1))

    //strandseq_phase.R
    strandseq_phase(input.join(ashleys_qc.out.ashleys_pass_folder), params.chromosomes)
    ch_versions = ch_versions.mix(strandseq_phase.out.versions.take(1))

    //patmat
    patmat(input.join(strandseq_phase.out.vcf), params.reference_genome)
    ch_versions = ch_versions.mix(patmat.out.versions.take(1))

    //Haplotagging the BAM file with whatshap
    prepare_whatshap_files(input.join(patmat.out.pofo_assignment_vcf)) 
    ch_versions = ch_versions.mix(prepare_whatshap_files.out.versions.take(1))

    whatshap_haplotag(prepare_whatshap_files.out.prep_files, params.reference_genome, params.reference_genome_idx)
    ch_versions = ch_versions.mix(whatshap_haplotag.out.versions.take(1))

    sort_haplotagged_bam(whatshap_haplotag.out.haplotagged_bam)
    ch_versions = ch_versions.mix(sort_haplotagged_bam.out.versions.take(1))

    //Compress BAM to CRAM
    compress_haplotagged_bam(sort_haplotagged_bam.out.haplotagged_sorted_bam, params.reference_genome)
    ch_versions = ch_versions.mix(sort_haplotagged_bam.out.versions.take(1))

    // Merge and publish logs
    ch_versions.unique().collectFile(name: params.logfile, storeDir: params.output_folder)
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
        ch_versions = Channel.empty()

        //concat_fastq
        strandseq_concat_fastq(input)
        ch_versions = ch_versions.mix(strandseq_concat_fastq.out.versions.take(1))

        //will process each fastq pair as its own job
        strandseq_cutadapt(strandseq_concat_fastq.out.fastqs.transpose())
        ch_versions = ch_versions.mix(strandseq_cutadapt.out.versions.take(1))

        //bowtie2 - outputs sam
        strandseq_bowtie(strandseq_cutadapt.out.fastqs, params.bowtie_index)
        ch_versions = ch_versions.mix(strandseq_bowtie.out.versions.take(1))

        //tidybam
        strandseq_tidy_reads(strandseq_bowtie.out.sam_file)
        ch_versions = ch_versions.mix(strandseq_tidy_reads.out.versions.take(1))
        
        //remove bams with very few reads
        strandseq_mark_duplicates(strandseq_tidy_reads.out.bam_file)
        ch_versions = ch_versions.mix(strandseq_mark_duplicates.out.versions.take(1))

        //split bams into pass/fail sets
        strand_pass_fail_cells(strandseq_mark_duplicates.out.bam_file.groupTuple(), params.min_reads_threshold)
        ch_versions = ch_versions.mix(strand_pass_fail_cells.out.versions.take(1))

    emit:
        passed_strand_seq_list = strand_pass_fail_cells.out.passed_strand_seq_list
        failed_strand_seq_list = strand_pass_fail_cells.out.failed_strand_seq_list
        pass_strandseq_bam_folder = strand_pass_fail_cells.out.pass_strandseq_bams
        fail_strandseq_bam_folder = strand_pass_fail_cells.out.fail_strandseq_bams
        ch_versions

}