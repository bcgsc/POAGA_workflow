#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process prepare_whatshap_files {
    tag "${patient}"
    cpus 16
    label "mem64"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"

    input:
        tuple val(patient),
            path(strandseq_fastq),
            path(clair_vcf),
            path(bam),
            path(bai),
            path(pofo_assignment_vcf)

	output:
        tuple val("${patient}"), path("PofO_Assignment_ForWhatsHap.vcf.gz"), path("PofO_Assignment_ForWhatsHap.vcf.gz.tbi"), path("whatshap.bam"), path("whatshap.bam.bai"), emit: prep_files
        path("versions.yml"), emit: versions

    script:
	"""
        prepare_whatshap_files.sh $pofo_assignment_vcf $bam ${task.cpus}

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$( tabix --version 2>&1 | grep "^tabix" )
        END_VERSIONS
    """
}

process whatshap_haplotag {
    tag "${patient}"
    cpus 4
    label "mem12"
    container "quay.io/biocontainers/whatshap:2.2--py310h0dbaff4_0"

    input:
        tuple val(patient),
            path(whatshap_vcf),
            path(whatshap_vcf_idx),
            path(whatshap_bam),
            path(whatshap_bai)
        path(reference_genome)
        path(reference_genome_idx)

	output:
        tuple val("${patient}"), path("PofO_HaploTagged.bam"), emit: haplotagged_bam
        path("versions.yml"), emit: versions

    script:
	"""
        # haplotag with whatshap
        whatshap haplotag --ignore-read-groups --tag-supplementary --output-threads ${task.cpus} --reference $reference_genome -o PofO_HaploTagged.bam $whatshap_vcf $whatshap_bam

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        whatshap: \$( whatshap --version )
        END_VERSIONS
    """
}


process sort_haplotagged_bam {
    tag "${patient}"
    cpus 16
    label "mem64"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"

    input:
        tuple val(patient),
            path(haplotagged_bam)

	output:
        tuple val("${patient}"), path("PofO_HaploTagged.sorted.bam"), path("PofO_HaploTagged.sorted.bam.bai"), emit: haplotagged_sorted_bam
        path("versions.yml"), emit: versions

    script:
	"""
        # Sort and index bam file
        samtools sort -@ ${task.cpus} ${haplotagged_bam} -o PofO_HaploTagged.sorted.bam
        samtools index -@ ${task.cpus} PofO_HaploTagged.sorted.bam

        # Log tool versions and container
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
    """
}
