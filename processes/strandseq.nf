#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ashleys_qc {
    tag "${patient}"
    cpus 24
    label "mem256"
    container "quay.io/pedrofeijao/ashleys-qc:v0.2.0"
    publishDir "${params.output_folder}/${patient}/ashleys_qc", mode: 'copy', pattern: '{fail,pass,quality.txt}'

    input:
        tuple val(patient),
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

// https://github.com/vahidAK/PatMat
// phased the strandseq data with the clair3 vcf
process strandseq_phase {
    tag "${patient}"
    cpus 16
    label "mem400"
    container "${params.singularity_folder}/patmat-1.4.0_dev.sif"

    publishDir "${params.output_folder}/${patient}/strandseqphase", mode: 'copy', pattern: "${patient}.chr-states.txt"
    publishDir "${params.output_folder}/${patient}/strandseqphase", mode: 'copy', pattern: "invertyper"
    publishDir "${params.output_folder}/${patient}/strandseqphase", mode: 'copy', pattern: "${patient}*composite.bam*"
    publishDir "${params.output_folder}/${patient}/strandseqphase", mode: 'copy', pattern: "*", enabled: params.debug

    input:
        tuple val(patient),
            path(clair3_vcf),
            path(ashleys_pass_bam_folder)
        val(chromosomes)

	output:
        tuple val("${patient}"), path("${patient}.strandseqphased.inv_aware.vcf"), emit: vcf
        tuple val("${patient}"), path("${patient}.chr-states.txt"), emit: chr_states
        tuple val("${patient}"), path("invertyper"), emit: invertyper_output
        tuple val("${patient}"), path("*"), emit: strandseq_phase_all_output
        tuple val("${patient}"),
            path(WC_OUT),
            path("${WC_OUT}.bai"),
            path(WW_OUT),
            path("${WW_OUT}.bai"), emit: composite_bams
        path("versions.yml"), emit: versions

    script:
        def STRANDSEQ_PHASE_VERSION="1.4.0dev"
        WW_OUT="${patient}_WW_composite.bam"
        WC_OUT="${patient}_WC_composite.bam"
	"""
        strandseq_phase.R ${params.extra_strandseq_options} \
        -p TRUE \
        -i ${ashleys_pass_bam_folder} \
        -o . \
        -t ${task.cpus} \
        -n ${patient} \
        --chromosomes ${chromosomes} \
        --inversion_list /opt/PatMat/Strand-seq/hanlon_2021_BMCgenomics_augmented.Merged.GRCh38.bed \
        --hard_mask /opt/PatMat/Strand-seq/hard_mask.GRCh38.bed \
        --soft_mask /opt/PatMat/Strand-seq/soft_mask.GRCh38.bed \
        ${clair3_vcf}

        report-chr-states.R ${patient}.chr-states.txt

        generate-composite-bam.R invertyper ${WW_OUT} ${WC_OUT}

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        strandseq_phase.R: ${STRANDSEQ_PHASE_VERSION}
        END_VERSIONS
    """
}

process strandseq_concat_fastq {
    cpus 4
    label "mem12"

    input:
        tuple val(patient),
            path(strandseq_folder)

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



process write_chr_phasing_summary {
    tag "${patient}"
    cpus 4
    label "mem12"
    container "pedrofeijao/vcfpy:0.13.8"
    publishDir "${params.output_folder}/${patient}/strandseqphase", mode: 'copy'

    input:
        tuple val(patient),
            path(chr_states)

	output:
        tuple val(patient), path(OUT), emit: chr_phasing_summary
        path("versions.yml"), emit: versions

    script:
    OUT="${patient}.chr_phasing_summary.txt"
	"""
        write_chr_states_summary.py ${chr_states} $OUT

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        python: \$( python --version 2>&1)
        END_VERSIONS
    """
}