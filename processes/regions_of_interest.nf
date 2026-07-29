#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    compress_vcf
} from '../processes/compress'

workflow extract_variants_of_interest {
    take:
    variants
    pofo_vcf
    pofo_scores
    pofo_dmr_status
    imprinted_regions
    haplotagged_bam

    main:
    // parse variants with TransVar
    variants_channel = Channel
        .fromPath(variants)
        .splitCsv(header: true)
        .map { row -> tuple(row.patient, row.variant) }
        .groupTuple()
        .transpose()

    annotate_variants(variants_channel)

    // compress and index vcf to search variants
    compress_vcf(pofo_vcf, "")

    // extract
    annotate_variants.out.variant_annotation.groupTuple().join(compress_vcf.out.compressed_vcf)
        | extract_variants

    // get Genes:
    parse_regions_channel = annotate_variants.out.variant_annotation
        .groupTuple()
        .combine(imprinted_regions)

    parse_regions_of_interest(parse_regions_channel)

    regions_and_bam_channel = parse_regions_of_interest.out.regions_of_interest.join(haplotagged_bam)

    subset_regions_of_interest_bam(regions_and_bam_channel)

    regions_and_methylation_channel = parse_regions_of_interest.out.regions_of_interest
        .join(pofo_scores)
        .join(pofo_dmr_status)

    subset_pofo_methylation_of_interest(regions_and_methylation_channel)

    // build tool versions channel
    ch_versions = annotate_variants.out.versions
        .take(1)
        .mix(compress_vcf.out.versions.take(1))
        .mix(extract_variants.out.versions.take(1))

    emit:
    ch_versions
    regions_of_interest = subset_pofo_methylation_of_interest.out.regions_of_interest_methylation
}

process annotate_variants {
    label "variants_of_interest"
    tag "${patient}"
    cpus 12
    debug true
    memory "24 GB"
    container "pedrofeijao/transvar:hg38"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
    tuple val(patient), val(variant)

    output:
    tuple val(patient), path("${variant.replace(">", "-").replace(":", ".")}.annotation.txt"), emit: variant_annotation
    path ("versions.yml"), emit: versions

    script:
    """
        transvar canno --refseq -i '${variant}' --longest > ${variant.replace(">", "-").replace(":", ".")}.annotation.txt

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        TransVar: \$( transvar --version 2>&1)
        END_VERSIONS
    """
}


process extract_variants {
    label "variants_of_interest"
    tag "${patient}"
    cpus 12
    memory "24 GB"
    container "pedrofeijao/vcfpy:0.13.8"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
    tuple val(patient), path(annotations), path(vcf), path(vcf_idx)

    output:
    tuple val(patient), path("variants_of_interest.vcf"), emit: variants_interest_vcf
    tuple val(patient), path("variants_of_interest_status.txt"), emit: variants_interest_txt
    path ("versions.yml"), emit: versions

    script:
    """
        extract_variants.py -v ${vcf} -a ${annotations} -o variants_of_interest.vcf -s variants_of_interest_status.txt

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        python: \$( python --version 2>&1)
        vcfpy: \$( pip list | grep vcfpy)
        END_VERSIONS
    """
}

process parse_regions_of_interest {
    label "variants_of_interest"
    tag "${patient}"
    cpus 12
    memory "24 GB"
    container "pedrofeijao/transvar:hg38"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
    tuple val(patient), path(annotations), path(imprinted_regions)

    output:
    tuple val(patient), path("regions_of_interest.bed"), emit: regions_of_interest
    path ("versions.yml"), emit: versions

    script:
    """
        parse_regions_of_interest.sh ${params.flanking_region_size} ${imprinted_regions} ${annotations} > regions_of_interest.bed

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        TransVar: \$( transvar --version 2>&1)
        END_VERSIONS
    """
}

process subset_regions_of_interest_bam {
    label "variants_of_interest"
    tag "${patient}"
    cpus 24
    memory "24 GB"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
    tuple val(patient), path(bed_file), path(bam_file), path(bai_file)

    output:
    tuple val(patient), path("regions_of_interest.bam"), path("regions_of_interest.bam.bai"), emit: regions_of_interest_bam
    path ("versions.yml"), emit: versions

    script:
    """
        samtools view -b -h -L ${bed_file} ${bam_file} > "regions_of_interest.bam"
        samtools index -@ ${task.cpus} regions_of_interest.bam

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
    """
}

process subset_pofo_methylation_of_interest {
    label "variants_of_interest"
    tag "${patient}"
    cpus 24
    memory "24 GB"
    container "pedrofeijao/vcfpy:0.13.8"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
    tuple val(patient), path(bed_file), path(pofo_scores), path(pofo_dmr_status)

    output:
    tuple val(patient), path("pofo_scores_regions_of_interest.tsv"), path("CpG-methylation-status-at-DMRs-regions_of_interest.tsv"), emit: regions_of_interest_methylation
    path ("versions.yml"), emit: versions

    script:
    """
        subset_pofo_methylation_of_interest.py \
            -b ${bed_file} \
            -s ${pofo_scores} \
            -d ${pofo_dmr_status} \
            -so pofo_scores_regions_of_interest.tsv \
            -do CpG-methylation-status-at-DMRs-regions_of_interest.tsv

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        python: \$( python --version 2>&1)
        END_VERSIONS

    """
}
