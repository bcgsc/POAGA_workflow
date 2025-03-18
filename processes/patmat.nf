#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// https://github.com/vahidAK/PatMat
//Run vahid's patmat to predict which variants came from each parent
process patmat {
    tag "${patient}"
    cpus 32
    label "mem400"
    container "${params.singularity_folder}/patmat-1.4.0_dev.sif"
    // Main output
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_PofO_Scores.tsv'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_PofO_Assigned_info.tsv'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_CpG-Methylation-Status-at-DMRs.tsv'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: "*", enabled: params.debug

    input:
        tuple val(patient),
            path(clair_phased_vcf),
            path(bam),
            path(bai),
            path(strandseq_phased_vcf)
        path(reference_genome)
        path(reference_genome_idx)


	output:
        tuple val(patient), path("*"), emit: patmat_all_output
        // VCF below is not published here, because bgzip+tabix later publishes it.
        tuple val(patient), path("patmat_PofO_Assigned.vcf"), emit: pofo_assignment_vcf
        tuple val(patient), path("patmat_PofO_Scores.tsv"), emit: pofo_scores
        tuple val(patient), path("patmat_CpG-Methylation-Status-at-DMRs.tsv"), emit: pofo_dmr_status
        path("versions.yml"), emit: versions

    script:
	"""
        tabix ${clair_phased_vcf}
        python3 /opt/PatMat/patmat/patmat.py \
        --strand_vcf ${strandseq_phased_vcf}  \
        --vcf ${clair_phased_vcf} \
        --bam ${bam} \
        --reference ${reference_genome} \
        --output "patmat" \
        -p ${task.cpus} \
        -dp 6 | tee stdout

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        patmat: \$( python3 /opt/PatMat/patmat/patmat.py --version )
        END_VERSIONS
    """
}

