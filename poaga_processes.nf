#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process log_pipeline_version {
    output:
        path("versions.yml"), emit: versions
    script:
    """
        cat <<-END_VERSIONS > versions.yml
        ---
        Pipeline version:"${workflow.manifest.version}"
        Nextflow version: "${nextflow.version}"
        END_VERSIONS
    """
}

// https://github.com/vahidAK/PatMat
// phased the strandseq data with the clair3 vcf
process strandseq_phase {
    tag "${patient}"
    cpus 12
    memory "400 GB"
    container "quay.io/pedrofeijao/patmat:1.3"
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
        /opt/PatMat-1.3.0/Strand-seq/strandseq_phase.R \
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
        strandseq_phase.R: "${STRANDSEQ_PHASE_VERSION}"
        END_VERSIONS
    """
}

//Run vahid's patmat to predict which variants came from each parent
process patmat {
    tag "${patient}"
    cpus 48
    memory "360 GB"
    container "quay.io/pedrofeijao/patmat:1.3"
    // Main output
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_PofO_Assignment.vcf'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'patmat_PofO_Scores.tsv'
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'snp_identity_check.txt'
    // Debug output
    publishDir "${params.output_folder}/${patient}/patmat", mode: 'copy', pattern: 'patmat_callDML.tsv.gz*', enabled: params.debug
    publishDir "${params.output_folder}/${patient}/patmat", mode: 'copy', pattern: 'patmat_CpG-Methylation-Status-at-DMRs.tsv', enabled: params.debug
    publishDir "${params.output_folder}/${patient}/patmat", mode: 'copy', pattern: 'patmat_HP1_HP2_PerReadInfo.tsv', enabled: params.debug
    publishDir "${params.output_folder}/${patient}/patmat", mode: 'copy', pattern: 'patmat_NonPofO*', enabled: params.debug
    publishDir "${params.output_folder}/${patient}/patmat", mode: 'copy', pattern: 'patmat_PofO_Assignment_', enabled: params.debug

    input:
        tuple val(patient),
            path(strandseq_bam_folder),
            path(clair_phased_vcf),
            path(bam),
            path(bai),
            path(strandseq_phased_vcf)
        path(reference_genome)

	output:
        tuple val("${patient}"), path("*"), emit: patmat_all_output
        tuple val("${patient}"), path("patmat_PofO_Assignment.vcf"), emit: pofo_assignment_vcf
        path("versions.yml"), emit: versions

    script:
	"""
         python3 /opt/PatMat-1.3.0/patmat/patmat.py \
         -sv ${strandseq_phased_vcf}  \
         -v ${clair_phased_vcf} \
         -mc ${bam} \
         -b  ${bam}  \
         -ref ${reference_genome} \
         -o "patmat" \
         -p ${task.cpus} \
         -dp 6 | tee stdout

        # Parse the SNP identity line:
        grep "reference heterozygous phased" stdout > snp_identity_check.txt

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
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
    cpus 48
    memory "48 GB"
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
    cpus 48
    memory "48 GB"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    publishDir "${params.output_folder}/${patient}/POAGA", mode: 'copy', pattern: 'PofO_HaploTagged.sorted.bam*'

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
    cpus 12
    memory "48 GB"
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
    cpus 12
    memory "12 GB"
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
    cpus 2
    memory "12 GB"
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
    cpus 12
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
    cpus 12
    memory "12 GB"
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
    cpus 48
    memory "48 GB"
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
