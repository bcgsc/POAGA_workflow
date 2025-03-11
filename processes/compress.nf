#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process compress_haplotagged_bam {
    tag "${patient}"
    cpus 12
    label "mem64"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy'


    input:
        tuple val(patient),
            path(bam),
            path(bai)
        path reference_genome

	output:
        tuple val("${patient}"), path("PofO_HaploTagged.sorted.cram"), path("PofO_HaploTagged.sorted.cram.crai"), emit: cram
        path("versions.yml"), emit: versions

    script:
	"""
        samtools view -@ ${task.cpus} -O CRAM,store_md=1,store_nm=1 -T ${reference_genome} $bam -o PofO_HaploTagged.sorted.cram
        samtools index -@ ${task.cpus} PofO_HaploTagged.sorted.cram

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
    """
}


process compress_vcf {
    tag "${patient}"
    cpus 12
    memory "24 GB"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    publishDir "${params.output_folder}/${patient}/${publish_dir}", mode: 'copy', overwrite: true, enabled: "${ publish_dir != '' }"

    input:
        tuple val(patient),
            path(vcf)
        val(publish_dir)

    output:
        tuple val(patient),
            path("${vcf}.gz"),
            path("${vcf}.gz.tbi"), emit: compressed_vcf
        path("versions.yml"), emit: versions
    script:
    """
        bgzip $vcf
        tabix ${vcf}.gz

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        bgzip: \$( bgzip --version 2>&1 | grep "^bgzip")
        tabix: \$( tabix --version 2>&1 | grep "^tabix" )
        END_VERSIONS
    """
}

