#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    extract_variants_of_interest
} from './processes/regions_of_interest'

/*
* Main workflow
*/

workflow {

    // HELP message
    if (params.help) {
        help = """
             |annotate_variants.nf: Extract variants of interest from Parent of Origin (POAGA) pipeline results.
             |
             |Usage: nextflow run extract_variants.nf --pofo_vcf <VCF> --pofo_scores <SCORES> --pofo_dmr <DMR> -variants <CSV> --output_folder <DIR>
             |
             |Required arguments:
             |  --pofo_vcf              Path to the PofO VCF "*PofO_Assignment.vcf" output(s). Use wildcards to input more than one sample.
             |  --pofo_scores           Path to the PofO Score "*PofO_Scores.tsv" output(s). Use wildcards to input more than one sample.
             |  --pofo_dmr              Path to the PofO DMRs "*CpG-Methylation-Status-at-DMRs.tsv" output(s). Use wildcards to input more than one sample.
             |  --input_csv             Location of the input CSV file. Each row corresponds to a sample to be processed.
             |                          The CSV file must have a header with 5 columns:
             |     - patient            Sample label.
             |     - strand_seq_fastq   Folder with the Strand-seq FastQ files for the sample.
             |     - clair3_vcf         Clair3 VCF output for the sample's Nanopore data.
             |     - bam                Nanopore BAM file.
             |     - bai                Index file for the Nanopore BAM file.
             |
             |  --variants             Variant CSV file(s). Each row corresponds to a variant of interest.
             |                          The CSV file must have a header with 2 columns:
             |     - patient            Sample label.
             |     - variant            Variant in HGVS cDNA notation.
             |
             |  --output_folder         Output folder to store results.
             |
             |""".stripMargin()

        // Print the help with the stripped margin and exit
        println(help)
        exit(0)
    }

    if (!params.input_csv) {
        error("Error: Please specify the input CSV file with the --input_csv parameter.")
    }
    if (!params.pofo_vcf) {
        error("Error: Please specify the input VCF file(s) with the --pofo_vcf parameter.")
    }
    if (!params.variants) {
        error("Error: Please specify the variants file with the --variants parameter.")
    }
    if (!params.output_folder) {
        error("Error: Please specify the output location with the --output_folder parameter.")
    }
    if (!params.pofo_scores) {
        error("Error: Please specify the input VCF file(s) with the --pofo_vcf parameter.")
    }
    if (!params.pofo_dmr) {
        error("Error: Please specify the input VCF file(s) with the --pofo_vcf parameter.")
    }
    if (!params.haplotagged_sorted_bam) {
        error("Error: Please specify the input VCF file(s) with the --pofo_vcf parameter.")
    }



    log.info(
        """\
Extract variants of interest and PofO from POAGA pipeline results.
===================================
input_csv           : ${params.input_csv}
pofo_vcf            : ${params.pofo_vcf}
output_folder	    : ${params.output_folder}
variants            : ${params.variants}

""".stripIndent()
    )

    bam_files = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.patient, row.bam, row.bai) }


    pofo_vcf = Channel
        .fromPath(params.pofo_vcf)
        .map { file -> tuple(file.baseName.replace("_PofO_Assignment", ""), file) }

    pofo_scores = Channel
        .fromPath(params.pofo_scores)
        .map { file -> tuple(file.baseName.replace("_PofO_Scores", ""), file) }

    pofo_dmr = Channel
        .fromPath(params.pofo_dmr)
        .map { file -> tuple(file.baseName.replace("_CpG-Methylation-Status-at-DMRs", ""), file) }

    extract_variants_of_interest(
        params.variants,
        pofo_vcf,
        pofo_scores,
        pofo_dmr,
        Channel.fromPath("${projectDir}/assets/Imprinted_DMR_List_V1.tsv"),
        bam_files,
    )
}
