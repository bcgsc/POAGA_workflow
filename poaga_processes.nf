#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
binaries
*/
// try to use public containers where available

// https://github.com/vahidAK/PatMat
// phased the strandseq data with the clair3 vcf
process strandseq_phase {
    tag "${patient}"
    cpus 12
    memory "400 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
    publishDir "${params.output_folder}/${patient}", mode: 'copy', pattern: "strandseq_phase"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair3_vcf),
            path(bam),
            path(bai),
            path(ashleys_pass_bam_folder)
        val(chromosomes)

	output:
        tuple val("${patient}"), path("strandseq_phase/${patient}.strandseqphased.inv_aware.vcf"), emit: vcf
        tuple val("${patient}"), path("strandseq_phase"), emit: strandseq_phase_all_output
        path("${task.process}.versions.log")


    script:
        def STRANDSEQ_PHASE_VERSION="1.3.0"
	"""
        /opt/PatMat-1.3.0/Strand-seq/strandseq_phase.R \
        -p TRUE \
        -i ${ashleys_pass_bam_folder} \
        -o strandseq_phase \
        -t ${task.cpus} \
        -n ${patient} \
        --chromosomes ${chromosomes} \
        ${clair3_vcf}

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        strandseq_phase.R: $STRANDSEQ_PHASE_VERSION
        END_VERSIONS
    """
}

//Run vahid's patmat to predict which variants came from each parent
process patmat {
    tag "${patient}"
    cpus 48
    memory "360 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
    publishDir "${params.output_folder}/${patient}/", mode: 'copy', pattern: 'patmat'
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"

    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair_phased_vcf),
            path(bam),
            path(bai),
            path(strandseq_phased_vcf)
        path(reference_genome)

	output:
        tuple val("${patient}"), path("patmat"), path("snp_identity_check.txt"), emit: patmat_all_output
        tuple val("${patient}"), path("patmat/patmat_PofO_Assignment.vcf"), emit: pofo_assignment_vcf
        path("${task.process}.versions.log")

    script:
	"""
        mkdir -p patmat
         python3 /opt/patmat.py \
         -sv ${strandseq_phased_vcf}  \
         -v ${clair_phased_vcf} \
         -mc ${bam} \
         -b  ${bam}  \
         -ref ${reference_genome} \
         -o "patmat/patmat" \
         -p ${task.cpus} \
         -dp 6 | tee stdout

        # Parse the SNP identity line:
        grep "reference heterozygous phased" stdout > snp_identity_check.txt

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        patmat: \$( /opt/patmat.py --version )
        END_VERSIONS
    """
}


process prepare_whatshap_files {
    tag "${patient}"
    cpus 48
    memory "48 GB"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"

    input:
        tuple val(patient),
            path(strandseq_fastq),
            path(clair_vcf),
            path(bam),
            path(bai),
            path(pofo_assignment_vcf)

	output:
        tuple val("${patient}"), path("PofO_Assignment_ForWhatsHap.vcf.gz"), path("PofO_Assignment_ForWhatsHap.vcf.gz.tbi"), path("whatshap.bam"), path("whatshap.bam.bai"), emit: prep_files
        path("${task.process}.versions.log")

    script:
	"""
        prepare_whatshap_files.sh $pofo_assignment_vcf $bam ${task.cpus}

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$( samtools --version | head -n 2 )
        tabix: \$( samtools --version | head -n2 )
        END_VERSIONS
    """
}

process whatshap_haplotag {
    tag "${patient}"
    cpus 48
    memory "48 GB"
    container "quay.io/biocontainers/whatshap:2.2--py310h0dbaff4_0"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"

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
        path("${task.process}.versions.log")

    script:
	"""
        # haplotag with whatshap
        whatshap haplotag --ignore-read-groups --tag-supplementary --output-threads ${task.cpus} --reference $reference_genome -o PofO_HaploTagged.bam $whatshap_vcf $whatshap_bam

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        whatshap: \$( whatshap --version )
        END_VERSIONS
    """
}


process sort_haplotagged_bam {
    tag "${patient}"
    cpus 48
    memory "48 GB"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    publishDir "${params.output_folder}/${patient}/", mode: 'copy', pattern: 'PofO_HaploTagged.sorted.bam*'
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"

    input:
        tuple val(patient),
            path(haplotagged_bam)

	output:
        tuple val("${patient}"), path("PofO_HaploTagged.sorted.bam"), path("PofO_HaploTagged.sorted.bam.bai"), emit: haplotagged_sorted_bam
        path("${task.process}.versions.log")

    script:
	"""
        # Sort and index bam file
        samtools sort -@ ${task.cpus} ${haplotagged_bam} -o PofO_HaploTagged.sorted.bam
        samtools index -@ ${task.cpus} PofO_HaploTagged.sorted.bam

        # Log tool versions and container
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$( samtools --version | head -n 2 )
        END_VERSIONS
    """
}



process strandseq_concat_fastq {
    cpus 12
    memory "48 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
    input:
        tuple val(patient),
            path(strandseq_folder),
            path(clair3_vcf),
            path(bam),
            path(bai)

	output:
        tuple val("${patient}"), path("merged_fastqs/*"), emit:fastqs
        path("${task.process}.versions.log")


    //Added the symlin to rename the input folder to something without underscores
    script:
	"""
        mkdir merged_fastqs 
        for i in \$(ls ${strandseq_folder}/*.fastq.gz | sed 's/${strandseq_folder}\\///' | cut -f1-2 -d_ | sort -u)
        do
            mkdir -p merged_fastqs/"\$i"
            cat ${strandseq_folder}/"\$i"_*R1*.fastq.gz > merged_fastqs/"\$i"/"\$i"_R1_001.fastq.gz
            cat ${strandseq_folder}/"\$i"_*R2*.fastq.gz > merged_fastqs/"\$i"/"\$i"_R2_001.fastq.gz
        done

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        END_VERSIONS
    """
}

process strandseq_cutadapt {
    cpus 12
    memory "48 GB"
    container "quay.io/biocontainers/cutadapt:2.6--py36h516909a_0"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
    input:
        tuple val(patient),
            path(fastq_folder)

	output:
        tuple val("${patient}"), path("*trimmed.1.fastq.gz"), path("*trimmed.2.fastq.gz"), emit: fastqs
        path("${task.process}.versions.log")

    script:
	"""
        read1=`ls ${fastq_folder}/*R1*`
        read2=`ls ${fastq_folder}/*R2*`
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
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        cutadapt: \$( cutadapt --version )
        END_VERSIONS
    """
}


process strandseq_bowtie {
    cpus 12
    memory "12 GB"
    container "quay.io/biocontainers/bowtie2:2.3.5.1--py27he513fc3_0"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
    input:
        tuple val(patient),
            path(fastq1),
            path(fastq2)
        path(bowtie_ref)

	output:
        tuple val("${patient}"), path("${fastq1.simpleName}.sam"), emit: sam_file
        path("${task.process}.versions.log")

    //had to use a folder for the reference and give the name here to get it available in the container 
    script:
	"""	
        bowtie2 -x ${bowtie_ref}/hg38_no_alt -p ${task.cpus} -1 ${fastq1} -2 ${fastq2} -S ${fastq1.simpleName}.sam


        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        bowtie2: \$( bowtie2 --version 2>&1 | grep bowtie2 )
        END_VERSIONS
    """
}


process strandseq_tidy_reads {
    cpus 2
    memory "12 GB"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
    input:
        tuple val(patient),
           path(sam_file)

	output:
        tuple val("${patient}"), path("${sam_file.simpleName}.sorted.bam"), emit: bam_file
        path("${task.process}.versions.log")

    //had to use a folder for the reference and give the name here to get it available in the container 
    script:
	"""	
        chrs=\$(printf "chr%s " {1..22} X Y)
        samtools view -bS ${sam_file} | samtools sort - -o ${sam_file.simpleName}.bam
        samtools index ${sam_file.simpleName}.bam 
        
        # Remove unmapped reads, secondary alignments, poor-quality reads, unwanted chromsomes, and then clean up any abandoned paired reads
        samtools view -F2052 -q10 -h ${sam_file.simpleName}.bam \$chrs | grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*chrEBV' | samtools sort -n - | samtools fixmate -O bam - - | samtools sort - | samtools view -bh -f1 -o ${sam_file.simpleName}.sorted.bam

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        samtools: \$( samtools --version | head -n 2 )
        END_VERSIONS
    """
}



//samtools says we're missing "ms" tags evening though fixmate is run on the bam
process strandseq_mark_duplicates {
    cpus 12
    memory "12 GB"
    container "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
    input:
        tuple val(patient),
           path(bam_file)

	output:
        tuple val("${patient}"), path("${bam_file.simpleName}.dup.bam"), emit:bam_file
        path("${task.process}.versions.log")

    //had to use a folder for the reference and give the name here to get it available in the container 
    script:
	"""	
       picard MarkDuplicates -I ${bam_file} -O ${bam_file.simpleName}.dup.bam -M ${bam_file.simpleName}.metrics

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        picard: \$( picard MarkDuplicates --version 2>&1 | grep Version )
        END_VERSIONS
    """
}

process strand_pass_fail_cells {
    cpus 12
    memory "12 GB"
    container "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
    publishDir "${params.output_folder}/${patient}/", mode: 'copy', pattern: "{pass_strandseq_bams,fail_strandseq_bams}"
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
    input:
        tuple val(patient),
           path(bam_files)
        val(min_reads_threshold)

	output:
        tuple val("${patient}"), path("pass.strandseq.files.txt"), emit: passed_strand_seq_list
        tuple val("${patient}"), path("fail.strandseq.files.txt"), emit: failed_strand_seq_list
        tuple val("${patient}"), path("pass_strandseq_bams"), emit: pass_strandseq_bams
        tuple val("${patient}"), path("fail_strandseq_bams"), emit: fail_strandseq_bams
        path("${task.process}.versions.log")
        
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
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        sambamba: \$( sambamba flagstat 2>&1 | head -n4 )
        END_VERSIONS
    """
}


process ashleys_qc {
    tag "${patient}"
    cpus 48
    memory "48 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
    publishDir "${params.output_folder}/${patient}/ashleys_qc", mode: 'copy', pattern: '{ashleys_fail,ashleys_pass,quality.txt}'
    publishDir "${params.output_folder}/${patient}/logs/", mode: 'copy', pattern: "${task.process}.versions.log"
  
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
        tuple val("${patient}"), path("ashleys_pass"), emit: ashleys_pass_folder
        tuple val("${patient}"), path("ashleys_fail"), emit: ashleys_fail_folder
        path("${task.process}.versions.log")
        
    script:
    def ASHLEYS_VERSION="0.2.0"
	"""
        ashleys_qc.sh ${strandseq_pass_bam_folder} ${pass_probability} ${task.cpus}

        # Log tool versions and containers
        cat <<-END_VERSIONS > ${task.process}.versions.log
        Process:"${task.process}"
        Container: "${task.container}"
        Ashleys: $ASHLEYS_VERSION
        END_VERSIONS
    """
}
