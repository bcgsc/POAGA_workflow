#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process log_pipeline_version {
    output:
        path("versions.yml"), emit: versions
    script:
    """
        cat <<-END_VERSIONS > versions.yml
        ---
        Pipeline version: ${workflow.manifest.version}
        Nextflow version: ${nextflow.version}
        END_VERSIONS
    """
}

process calc_het_hom_discordance {
    tag "${patient}"
    cpus 2
    label "mem12"
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    publishDir "${params.output_folder}/${patient}/qc_metrics", mode: 'copy'

    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair_phased_vcf),
            path(bam),
            path(bai),
            path(strandseq_phased_vcf)

	output:
        tuple val("${patient}"), path("het_hom_discordance.txt"), emit: het_hom_discordance
        path("versions.yml"), emit: versions

    script:
    """
        bgzip ${strandseq_phased_vcf}
        tabix ${strandseq_phased_vcf}.gz
        tabix ${clair_phased_vcf}
        echo "DiscordantSites\tTotalSites\tDiscordance" > het_hom_discordance.txt
        bcftools gtcheck -e 0 -g ${clair_phased_vcf} ${strandseq_phased_vcf}.gz \
            | grep "^DC" \
            | awk  '{ printf "%d\t%d\t%.4f", \$4, \$6, \$4 / \$6  }' \
            >> het_hom_discordance.txt

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        bctools: \$( bcftools --version 2>&1 | grep "^bcftools" )
        END_VERSIONS
   """

}
