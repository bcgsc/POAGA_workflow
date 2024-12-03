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

// https://github.com/vahidAK/PatMat
// phased the strandseq data with the clair3 vcf
process strandseq_phase {
    tag "${patient}"
    cpus 16
    label "mem400"
    container "${params.singularity_folder}/patmat-dev.sif"
    publishDir "${params.output_folder}/${patient}/strandseqphase", mode: 'copy', pattern: "${patient}.strandseqphased.inv_aware.vcf"
    publishDir "${params.output_folder}/${patient}/strandseqphase", mode: 'copy', pattern: "*", enabled: params.debug

    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair3_vcf),
            path(bam),
            path(bai),
            path(ashleys_pass_bam_folder)
        val(chromosomes)

	output:
        tuple val("${patient}"), path("${patient}.strandseqphased.inv_aware.vcf"), emit: vcf
        tuple val("${patient}"), path("*"), emit: strandseq_phase_all_output
        path("versions.yml"), emit: versions

    script:
        def STRANDSEQ_PHASE_VERSION="1.3.0"
	"""
        /opt/PatMat/Strand-seq/strandseq_phase.R \
        -p TRUE \
        -i ${ashleys_pass_bam_folder} \
        -o . \
        -t ${task.cpus} \
        -n ${patient} \
        --chromosomes ${chromosomes} \
        ${clair3_vcf}

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        strandseq_phase.R: ${STRANDSEQ_PHASE_VERSION}
        END_VERSIONS
    """
}

//Run vahid's patmat to predict which variants came from each parent
process patmat {
    tag "${patient}"
    cpus 32
    label "mem400"
    container "${params.singularity_folder}/patmat-dev.sif"
    // Main output
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_PofO_Assigned.vcf'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_PofO_Scores.tsv'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_PofO_Assigned_info.tsv'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_CpG-Methylation-Status-at-DMRs.tsv'

    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair_phased_vcf),
            path(bam),
            path(bai),
            path(strandseq_phased_vcf)
        path(reference_genome)
        path(reference_genome_idx)


	output:
        tuple val(patient), path("*"), emit: patmat_all_output
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



process strandseq_concat_fastq {
    cpus 4
    label "mem12"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"

    input:
        tuple val(patient),
            path(strandseq_folder),
            path(clair3_vcf),
            path(bam),
            path(bai)

	output:
        tuple val("${patient}"), path("merged_fastqs/*"), emit:fastqs
        path("versions.yml"), emit: versions

    //Added the symlin to rename the input folder to something without underscores
    script:
	"""
        mkdir merged_fastqs
        for r1 in ${strandseq_folder}/*_R1*.fastq.gz
        do
            r2=\${r1//_R1/_R2}  # get R2 filename
            prefix=\$( basename \$r1 | sed -r 's/_R1.*fastq.gz//')  # prefix (sample name before R{1,2})
            OUT=merged_fastqs/\$prefix
            mkdir -p \$OUT
            cp \$r1 \$OUT/\$(basename \$r1)
            cp \$r2 \$OUT/\$(basename \$r2)
        done

        # Log tool versions and container
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        END_VERSIONS
    """
}

process strandseq_cutadapt {
    cpus 4
    label "mem12"
    container "quay.io/biocontainers/cutadapt:2.6--py36h516909a_0"

    input:
        tuple val(patient),
            path(fastq_folder)

	output:
        tuple val("${patient}"), path("*trimmed.1.fastq.gz"), path("*trimmed.2.fastq.gz"), emit: fastqs
        path("versions.yml"), emit: versions

    script:
	"""
        read1=`ls ${fastq_folder}/*_R1*`
        read2=`ls ${fastq_folder}/*_R2*`
        cell_name=`ls -d ${fastq_folder}`

        cutadapt \
            -g AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
            -G CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
            -o \${cell_name}_trimmed.1.fastq.gz \
            -p \${cell_name}_trimmed.2.fastq.gz \
            \${read1} \
            \${read2} \
            -m 30 \
            -q 15

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        cutadapt: \$( cutadapt --version )
        END_VERSIONS
    """
}


process strandseq_bowtie {
    cpus 16
    label "mem12"
    container "quay.io/biocontainers/bowtie2:2.3.5.1--py27he513fc3_0"

    input:
        tuple val(patient),
            path(fastq1),
            path(fastq2)
        path(bowtie_ref)

	output:
        tuple val("${patient}"), path("${fastq1.simpleName}.sam"), emit: sam_file
        path("versions.yml"), emit: versions

    //had to use a folder for the reference and give the name here to get it available in the container
    script:
	"""
        bowtie2 -x ${bowtie_ref}/hg38_no_alt -p ${task.cpus} -1 ${fastq1} -2 ${fastq2} -S ${fastq1.simpleName}.sam


        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        END_VERSIONS
    """
}


process strandseq_tidy_reads {
    cpus 4
    label "mem12"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"

    input:
        tuple val(patient),
           path(sam_file)

	output:
        tuple val("${patient}"), path("${sam_file.simpleName}.sorted.bam"), emit: bam_file
        path("versions.yml"), emit: versions

    //had to use a folder for the reference and give the name here to get it available in the container
    script:
	"""
        chrs=\$(printf "chr%s " {1..22} X Y)
        samtools view -bS ${sam_file} | samtools sort - -o ${sam_file.simpleName}.bam
        samtools index ${sam_file.simpleName}.bam

        # Remove unmapped reads, secondary alignments, poor-quality reads, unwanted chromsomes, and then clean up any abandoned paired reads
        samtools view -F2052 -q10 -h ${sam_file.simpleName}.bam \$chrs | grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*chrEBV' | samtools sort -n - | samtools fixmate -O bam - - | samtools sort - | samtools view -bh -f1 -o ${sam_file.simpleName}.sorted.bam

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
    """
}



//samtools says we're missing "ms" tags evening though fixmate is run on the bam
process strandseq_mark_duplicates {
    cpus 4
    memory "12 GB"
    container "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"

    input:
        tuple val(patient),
           path(bam_file)

	output:
        tuple val("${patient}"), path("${bam_file.simpleName}.dup.bam"), emit:bam_file
        path("versions.yml"), emit: versions

    //had to use a folder for the reference and give the name here to get it available in the container
    script:
	"""
        picard MarkDuplicates -I ${bam_file} -O ${bam_file.simpleName}.dup.bam -M ${bam_file.simpleName}.metrics

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        picard: \$( picard MarkDuplicates --version 2>&1 | grep Version )
        END_VERSIONS
    """
}

process strand_pass_fail_cells {
    cpus 4
    label "mem12"
    container "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"

    input:
        tuple val(patient),
           path(bam_files)
        val(min_reads_threshold)

	output:
        tuple val("${patient}"), path("pass.strandseq.files.txt"), emit: passed_strand_seq_list
        tuple val("${patient}"), path("fail.strandseq.files.txt"), emit: failed_strand_seq_list
        tuple val("${patient}"), path("pass_strandseq_bams"), emit: pass_strandseq_bams
        tuple val("${patient}"), path("fail_strandseq_bams"), emit: fail_strandseq_bams
        path("versions.yml"), emit: versions

    stub:
    """
        mkdir pass_strandseq_bams
        mkdir fail_strandseq_bams
        touch pass.strandseq.files.txt
        touch fail.strandseq.files.txt
        touch versions.yml
        cp *bam pass_strandseq_bams
        cd pass_strandseq_bams
        ls *bam | xargs -i sambamba index -t ${task.cpus} {}
    """

    //had to use a folder for the reference and give the name here to get it available in the container
    script:

	"""
        mkdir pass_strandseq_bams
        mkdir fail_strandseq_bams

        ls *bam | xargs -i echo "sambamba flagstat {} > {}.flagstat" | bash

        grep -H total *.flagstat | awk '{ print \$1 }' | awk -F ":" '\$2>$min_reads_threshold { print \$1 }' | sed 's/.flagstat//' > pass.strandseq.files.txt
        grep -H total *.flagstat | awk '{ print \$1 }' | awk -F ":" '\$2<=$min_reads_threshold { print \$1 }' | sed 's/.flagstat//' > fail.strandseq.files.txt

        cat pass.strandseq.files.txt | xargs -i cp {} pass_strandseq_bams
        cat fail.strandseq.files.txt | xargs -i cp {} fail_strandseq_bams

        cd pass_strandseq_bams
        ls *bam | xargs -i sambamba index {}
        cd ..

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        sambamba: \$( sambamba 2>&1 | grep "^sambamba" )
        END_VERSIONS
    """
}


process ashleys_qc {
    tag "${patient}"
    cpus 24
    label "mem256"
    container "quay.io/pedrofeijao/ashleys-qc:v0.2.0"
    publishDir "${params.output_folder}/${patient}/ashleys_qc", mode: 'copy', pattern: '{fail,pass,quality.txt}'

    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair3_vcf),
            path(bam),
            path(bai),
            path(strandseq_pass_bam_folder)
        val(pass_probability)

	output:
        tuple val("${patient}"), path("quality.txt")
        tuple val("${patient}"), path("pass"), emit: ashleys_pass_folder
        tuple val("${patient}"), path("fail"), emit: ashleys_fail_folder
        path("versions.yml"), emit: versions

    stub:
    """
        mkdir pass
        mkdir fail
        cp ${strandseq_pass_bam_folder}/*bam ${strandseq_pass_bam_folder}/*bai pass
        echo "Skipped Ashleys QC" > quality.txt
        touch versions.yml
    """
    script:
    def ASHLEYS_VERSION="0.2.0"
	"""
        ashleys_qc.sh ${strandseq_pass_bam_folder} ${pass_probability} ${task.cpus}

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        Ashleys: $ASHLEYS_VERSION
        END_VERSIONS
    """
}

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

process annotate_variants {
    tag "${patient}"
    cpus 12
    debug true
    memory "24 GB"
    container "pedrofeijao/transvar:hg38"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
        tuple val(patient),
            val(variant)

	output:
        tuple val(patient), path("${variant.replace(">","-").replace(":",".")}.annotation.txt"), emit: variant_annotation
        path("versions.yml"), emit: versions

    script:
	"""
        transvar canno --refseq -i '${variant}' --longest > ${variant.replace(">","-").replace(":",".")}.annotation.txt

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        TransVar: \$( transvar --version 2>&1)
        END_VERSIONS
    """
}

process compress_vcf {
    tag "${patient}"
    cpus 12
    debug true
    memory "24 GB"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"

    input:
        tuple val(patient),
            path(vcf)

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


process extract_variants {
    tag "${patient}"
    cpus 12
    memory "24 GB"
    container "pedrofeijao/vcfpy:0.13.8"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
        tuple val(patient),
            path(annotations),
            path(vcf),
            path(vcf_idx)

	output:
        tuple val(patient), path("variants_of_interest.vcf"), emit: variants_interest_vcf
        tuple val(patient), path("variants_of_interest_status.txt"), emit: variants_interest_txt
        path("versions.yml"), emit: versions

    script:
	"""
        extract_variants.py -v $vcf -a $annotations -o variants_of_interest.vcf -s variants_of_interest_status.txt

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
    tag "${patient}"
    cpus 12
    memory "24 GB"
    container "pedrofeijao/transvar:hg38"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
        tuple val(patient),
            path(annotations),
            path(imprinted_regions)

	output:
        tuple val(patient), path("regions_of_interest.bed"), emit: regions_of_interest
        path("versions.yml"), emit: versions

    script:
	"""
        parse_regions_of_interest.sh ${params.flanking_region_size} ${imprinted_regions} $annotations > regions_of_interest.bed

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
    tag "${patient}"
    cpus 24
    memory "24 GB"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
        tuple val(patient),
            path(bed_file),
            path(bam_file),
            path(bai_file)

	output:
        tuple val(patient), path("regions_of_interest.bam"), path("regions_of_interest.bam.bai"), emit: regions_of_interest_bam
        path("versions.yml"), emit: versions

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
    tag "${patient}"
    cpus 24
    memory "24 GB"
    container "pedrofeijao/vcfpy:0.13.8"
    publishDir "${params.output_folder}/${patient}/POAGA/variants_of_interest", mode: 'copy'

    input:
        tuple val(patient),
            path(bed_file),
            path(pofo_scores),
            path(pofo_dmr_status)

	output:
        tuple val(patient),
            path("pofo_scores_regions_of_interest.tsv"),
            path("CpG-methylation-status-at-DMRs-regions_of_interest.tsv"), emit: regions_of_interest_methylation
        path("versions.yml"), emit: versions

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
        variants_channel = Channel.fromPath(variants)
                .splitCsv(header:true)
                .map{ row-> tuple(row.patient, row.variant)}
                .groupTuple()
                .transpose()

        annotate_variants(variants_channel)

        // compress and index vcf to search variants
        compress_vcf(pofo_vcf)

        // extract
        annotate_variants.out.variant_annotation
            .groupTuple()
            .join(compress_vcf.out.compressed_vcf)
            | extract_variants

        // get Genes:
        parse_regions_channel =
            annotate_variants.out.variant_annotation.groupTuple()
            .combine(imprinted_regions)

        parse_regions_of_interest(parse_regions_channel)

        regions_and_bam_channel = parse_regions_of_interest.out.regions_of_interest
            .join(haplotagged_bam)

        subset_regions_of_interest_bam(regions_and_bam_channel)

        regions_and_methylation_channel = parse_regions_of_interest.out.regions_of_interest
            .join(pofo_scores)
            .join(pofo_dmr_status)

        subset_pofo_methylation_of_interest(regions_and_methylation_channel)

        // build tool versions channel
        ch_versions = annotate_variants.out.versions.take(1)
            .mix(compress_vcf.out.versions.take(1))
            .mix(extract_variants.out.versions.take(1))
    emit:
        ch_versions
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
