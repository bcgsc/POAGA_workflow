process mpileup_strandseq {
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    input:
    tuple val(patient), path(clair3_vcf), path(strandseq_pass_bam_folder)
    path region

    output:
    tuple val(patient), path(NANO_OUT), path(SSEQ_OUT), emit: vcf_files
    path ("versions.yml"), emit: versions

    script:
    SSEQ_OUT = "${patient}.sseq.subset.vcf.gz"
    NANO_OUT = "${patient}.nano.subset.vcf.gz"
    """
        bcftools index ${clair3_vcf}
        bcftools view -Oz -R ${region} -o ${NANO_OUT} ${clair3_vcf}
        bcftools mpileup --no-reference -Q 30 -R ${region} ${strandseq_pass_bam_folder}/*bam \
            | bcftools view -G -AA -H -a \
            | awk '{ if (\$5 !="." ) print \$1"\t"\$2"\t"\$5}' \
            | gzip > ${SSEQ_OUT}

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        bctools: \$( bcftools --version 2>&1 | grep "^bcftools" )
        END_VERSIONS
        """
}

process calc_discordance {
    label "mem32"
    container "quay.io/biocontainers/pandas:2.2.1"
    publishDir "${params.output_folder}/${patient}/qc_metrics", mode: 'copy'

    input:
    tuple val(patient), path(nano_vcf), path(sseq_vcf)

    output:
    tuple val(patient), path(OUT), emit: snp_discordance
    path ("versions.yml"), emit: versions

    script:
    OUT = "${patient}.snp_discordance.txt"
    """
        compare_gt.py -n ${nano_vcf} -s ${sseq_vcf} -o ${OUT}

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        python: \$( python --version 2>&1)
        pandas: \$( pip list | grep pandas)
        END_VERSIONS

        """
}
