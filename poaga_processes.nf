#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
binaries
*/

process find_fastq_files {
    tag "${library}"
    cpus 2
    memory "8 GB"
  
    input:
        tuple val(flowcell),
            val(library),
            path(fast5_folder),
            path(fastq_folder),
            path(sequencing_summary)

	output:
        tuple val("${flowcell}"), val("${library}"), path("*gz")

    script:
	"""
       ls ${fastq_folder}/*pass*.gz | xargs -i ln -s {}
    """
}

