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
    cpus 48
    memory "200 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
    publishDir "${params.output_folder}/${patient}/", mode: 'copy'
  
    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair3_vcf),
            path(bam),
            path(bai),
            path(ashleys_pass_bam_folder)

	output:
        tuple val("${patient}"), path("strandseq_phase/${patient}.strandseqphased.inv_aware.vcf"), emit: vcf
        tuple val("${patient}"), path("strandseq_phase"), emit: strandseq_phase_all_output


    script:
	"""
        /opt/PatMat-1.3.0/Strand-seq/strandseq_phase.R \
        -p TRUE \
        -i ${ashleys_pass_bam_folder} \
        -o strandseq_phase \
        -t 24 \
        -n ${patient} \
        ${clair3_vcf}
    """
}

//Run vahid's patmat to predict which variants came from each parent
process patmat {
    tag "${patient}"
    cpus 48
    memory "48 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
    publishDir "${params.output_folder}/${patient}/", mode: 'copy'
  
    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair_phased_vcf),
            path(bam),
            path(bai),
            path(strandseq_phased_vcf)
        path(reference_genome)

	output:
        tuple val("${patient}"), path("patmat_out*")

    script:
	"""
         python3 /opt/patmat.py \
         -sv ${strandseq_phased_vcf}  \
         -v ${clair_phased_vcf} \
         -mc ${bam} \
         -b  ${bam}  \
         -ref ${reference_genome} \
         -o patmat_out \
         -p 48 \
         -dp 10
    """
}



process strandseq_concat_fastq {
    cpus 12
    memory "48 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
  
    input:
        tuple val(patient),
            path(strandseq_folder),
            path(clair3_vcf),
            path(bam),
            path(bai)

	output:
        tuple val("${patient}"), path("merged_fastqs/*")


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
    """
}

process strandseq_cutadapt {
    cpus 12
    memory "48 GB"
    container "quay.io/biocontainers/cutadapt:2.6--py36h516909a_0"
  
    input:
        tuple val(patient),
            path(fastq_folder)

	output:
        tuple val("${patient}"), path("*trimmed.1.fastq.gz"), path("*trimmed.2.fastq.gz")

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
    """
}


process strandseq_bowtie {
    cpus 12
    memory "12 GB"
    container "quay.io/biocontainers/bowtie2:2.3.5.1--py27he513fc3_0"
  
    input:
        tuple val(patient),
            path(fastq1),
            path(fastq2)
        path(bowtie_ref)

	output:
        tuple val("${patient}"), path("${fastq1.simpleName}.sam")

    //had to use a folder for the reference and give the name here to get it available in the container 
    script:
	"""	
        bowtie2 -x ${bowtie_ref}/hg38_no_alt -p 12 -1 ${fastq1} -2 ${fastq2} -S ${fastq1.simpleName}.sam
    """
}


process strandseq_tidy_reads {
    cpus 2
    memory "12 GB"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
  
    input:
        tuple val(patient),
           path(sam_file)

	output:
        tuple val("${patient}"), path("${sam_file.simpleName}.sorted.bam")

    //had to use a folder for the reference and give the name here to get it available in the container 
    script:
	"""	
        chrs=\$(printf "chr%s " {1..22} X Y)
        samtools view -bS ${sam_file} | samtools sort - -o ${sam_file.simpleName}.bam
        samtools index ${sam_file.simpleName}.bam 
        
        # Remove unmapped reads, secondary alignments, poor-quality reads, unwanted chromsomes, and then clean up any abandoned paired reads
        samtools view -F2052 -q10 -h ${sam_file.simpleName}.bam \$chrs | grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*chrEBV' | samtools sort -n - | samtools fixmate -O bam - - | samtools sort - | samtools view -bh -f1 -o ${sam_file.simpleName}.sorted.bam
    """
}



//samtools says we're missing "ms" tags evening though fixmate is run on the bam
process strandseq_mark_duplicates {
    cpus 12
    memory "12 GB"
    container "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
  
    input:
        tuple val(patient),
           path(bam_file)

	output:
        tuple val("${patient}"), path("${bam_file.simpleName}.dup.bam")

    //had to use a folder for the reference and give the name here to get it available in the container 
    script:
	"""	
       picard MarkDuplicates -I ${bam_file} -O ${bam_file.simpleName}.dup.bam -M ${bam_file.simpleName}.metrics
    """
}

process strand_pass_fail_cells {
    cpus 12
    memory "12 GB"
    container "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
    publishDir "${params.output_folder}/${patient}/", mode: 'copy'
  
    input:
        tuple val(patient),
           path(bam_files)

	output:
        tuple val("${patient}"), path("pass.strandseq.files.txt"), emit: passed_strand_seq_list
        tuple val("${patient}"), path("fail.strandseq.files.txt"), emit: failed_strand_seq_list
        tuple val("${patient}"), path("pass_strandseq_bams"), emit: pass_strandseq_bams
        tuple val("${patient}"), path("fail_strandseq_bams"), emit: fail_strandseq_bams
        
    //had to use a folder for the reference and give the name here to get it available in the container 
    script:
	"""	
        mkdir pass_strandseq_bams
        mkdir fail_strandseq_bams

        ls *bam | xargs -i echo "sambamba flagstat {} > {}.flagstat" | bash 

        grep -H total *.flagstat | awk '{ print \$1 }' | awk -F ":" '\$2>10000 { print \$1 }' | sed 's/.flagstat//' > pass.strandseq.files.txt
        grep -H total *.flagstat | awk '{ print \$1 }' | awk -F ":" '\$2<=10000 { print \$1 }' | sed 's/.flagstat//' > fail.strandseq.files.txt

        cat pass.strandseq.files.txt | xargs -i cp {} pass_strandseq_bams
        cat fail.strandseq.files.txt | xargs -i cp {} fail_strandseq_bams

        cd pass_strandseq_bams
        ls *bam | xargs -i sambamba index {}

    """
}


process ashleys_qc {
    tag "${patient}"
    cpus 48
    memory "48 GB"
    container "/gsc/www/bcgsc.ca/downloads/rcorbett/poaga.sif"
    publishDir "${params.output_folder}/${patient}/", mode: 'copy'
  
    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair3_vcf),
            path(bam),
            path(bai),
            path(strandseq_pass_bam_folder)

	output:
        tuple val("${patient}"), path("quality.txt")
        tuple val("${patient}"), path("ashleys_pass"), emit: ashleys_pass_folder
        tuple val("${patient}"), path("ashleys_fail"), emit: ashleys_fail_folder
        
    script:
	"""
        ashleys.py -j 48 \
        features \
        -f ${strandseq_pass_bam_folder} \
        -w 5000000 2000000 1000000 800000 600000 400000 200000 \
        -o ./features.tsv

        ashleys.py predict -p ./features.tsv -o ./quality.txt -m /opt/ashleys-qc/models/svc_default.pkl

        mkdir ashleys_pass
        mkdir ashleys_fail

        awk 'NR>1 && \$2==1{ print \$1 }' ./quality.txt | xargs -i cp ${strandseq_pass_bam_folder}/{} ${strandseq_pass_bam_folder}/{}.bai ashleys_pass
        awk 'NR>1 && \$2==0{ print \$1 }' ./quality.txt | xargs -i cp ${strandseq_pass_bam_folder}/{} ${strandseq_pass_bam_folder}/{}.bai ashleys_fail
        
    """
}
