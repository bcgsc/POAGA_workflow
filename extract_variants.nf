#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    extract_variants_of_interest
} from './poaga_processes.nf'

// HELP message
if ( params.help ) {
    help = """
             |annotate_variants.nf: Extract variants of interest from Parent of Origin (POAGA) pipeline results.
             |
             |Usage: nextflow run extract_variants.nf --pofo_vcf <VCF> --variants <CSV> --output_folder <DIR>
             |
             |Required arguments:
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

if (!params.pofo_vcf) {
    error("Error: Please specify the input VCF file(s) with the --pofo_vcf parameter.")
}
if (!params.variants) {
    error("Error: Please specify the input CSV file with the --input_csv parameter.")
}
if (!params.output_folder) {
    error("Error: Please specify the output location with the --output_folder parameter.")
}


log.info """\
Extract variants of interest and PofO from POAGA pipeline results.
===================================
output_folder	    : ${params.output_folder}
variants            : ${params.variants}
"""
.stripIndent()


/*
* Main workflow
*/

workflow {

    pofo_vcf = Channel.fromPath(params.pofo_vcf)
                    .map { vcf -> tuple(vcf.baseName.replace("_PofO_Assignment",""), vcf)}

    extract_variants_of_interest(params.variants, pofo_vcf)
}
