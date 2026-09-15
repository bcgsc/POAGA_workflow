#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    compress_haplotagged_bam ;
    compress_vcf as compress_patmat_vcf ;
    compress_vcf as compress_strandseq_phase_vcf ;
    filter_vcf_annotations
} from './processes/compress'

include {
    calc_het_hom_discordance ;
    log_pipeline_version
} from './processes/qc'

include {
    extract_variants_of_interest
} from './processes/regions_of_interest'

include {
    patmat
} from './processes/patmat'

include {
    prepare_whatshap_files ;
    sort_haplotagged_bam ;
    whatshap_haplotag
} from './processes/whatshap'

include {
    ashleys_qc ;
    convert_fastq_pair_to_single ;
    strand_pass_fail_cells ;
    strandseq_bowtie ;
    strandseq_cutadapt ;
    strandseq_mark_duplicates ;
    strandseq_phase ;
    strandseq_tidy_reads ;
    write_chr_phasing_summary
} from './processes/strandseq'

include {
    mpileup_strandseq ;
    calc_discordance
} from './processes/snp_concordance'

/*
* Main workflow
*/

workflow {


    // HELP message
    if (params.help) {
        help = """
             |poaga_processes.nf: Parent of Origin (POAGA) pipeline using strand-seq and nanopore data.
             |
             |Usage: nextflow run poaga_processes.nf --input_csv <CSV> --output_folder <DIR>
             |
             |Required arguments:
             |  --input_csv             Location of the input CSV file. Each row corresponds to a sample to be processed.
             |                          The CSV file must have a header with either 6 or 7 columns:
             |     - patient            Sample label.
             |     - strand_seq_fastq   Folder containing Strand-seq FastQ files. Optional if strand_seq_fastq_list is given.
             |     - strand_seq_fastq_list Optional semicolon-separated list of individual Strand-seq FastQ file paths.
             |                            Takes precedence when both are supplied.
             |     - clair3_vcf         Clair3 VCF output for the sample's Nanopore data.
             |     - sniffles2_vcf      sniffles2 SV VCF output for the sample's Nanopore data.
             |     - bam                Nanopore BAM file.
             |     - bai                Index file for the Nanopore BAM file.
             |
             |  --output_folder         Output folder to store POAGA results.
             |
             |Optional arguments:
             |  --variants              Variant CSV file(s). Each row corresponds to a variant of interest.
             |                          The CSV file must have a header with 2 columns:
             |     - patient            Sample label, matching the 'input_csv' label(s).
             |     - variant            Variant in HGVS cDNA notation.
             |  --idmr_regions_file     Path to iDMR list TSV file used for variant-of-interest extraction.
             |                [default: ${params.idmr_regions_file}]
             |  --min_reads_threshold   Minimum number of high quality reads to pass a StrandSeq library file.
             |                [default: ${params.min_reads_threshold}]
             |  --chromosomes           Comma-separated list of chromosomes to process.
             |                [default: ${params.chromosomes}]
             |  --min_ashleys_qc_prob   Minimum probability score for a Strand-seq library to pass the ASHLEYS QC step.
             |                [default: ${params.min_ashleys_qc_prob}]
             |  -w                      The NextFlow work directory.
             |                [default: ${workDir}]
             |  -profile                Nextflow config profile(s) with pre-defined parameters (comma-separated profiles to use more than one).
             |                Available profiles include 'hg38_no_alt' (hg38 reference genome) and 'test' (run minimal dataset test)
             |
             |
             |Useful Nextflow arguments:
             |  -resume                Does not execute steps that have already been executed with the same parameters. Useful for resuming
             |                         failed runs or re-runs with some parameter changes.
             |  -with-report           Creates an HTML report with many nextflow execution metrics.
             |  -with-trace            Creates a tab-separated file with time, cpu and memory stats of each process.
             |  -with-timeline         Creates an HTML execution timeline of the pipeline
             |""".stripMargin()

        // Print the help with the stripped margin and exit
        log.info(help)
        exit(0)
    }

    if (!params.input_csv) {
        error("Error: Please specify the input CSV file with the --input_csv parameter.")
    }
    if (!params.output_folder) {
        error("Error: Please specify the output location with the --output_folder parameter.")
    }


    log.info(
        """\
Run the POAGA pipeline using strand-seq and nanopore data
===================================
output_folder	    : ${params.output_folder}
samples_file        : ${params.input_csv}
""".stripIndent()
    )

    input = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .multiMap { row ->
            def strand_fastq = row.get('strand_seq_fastq')
            def strand_fastq_list = row.get('strand_seq_fastq_list')

            strandseq: tuple(row.patient, strand_fastq, strand_fastq_list)
            nanopore_vcf: tuple(row.patient, row.clair3_vcf)
            nanopore_sv_vcf: tuple(row.patient, row.sniffles2_vcf)
            nanopore_bam: tuple(row.patient, row.bam, row.bai)
        }

    // "transform" the folder input in the channel to a fastq list input.
    strandseq_reads = build_strandseq_fastq_input(input.strandseq)

    // Log pipeline version:
    ch_versions = log_pipeline_version()

    //implemented below. Does the basic aligning, tidying and read count qc for strand-seq
    strandseq(strandseq_reads)
    ch_versions = ch_versions.mix(strandseq.out.ch_versions)

    // Filter clair3 VCF annotations, might break some tools:
    nanopore_vcf = filter_vcf_annotations(input.nanopore_vcf).filtered_vcf

    // SNP concordance:
    snp_concordance_channel = nanopore_vcf.join(strandseq.out.pass_strandseq_bam_folder)
    snp_region_file = channel.fromPath(params.snp_region, checkIfExists: true)
    mpileup_vcf_files = mpileup_strandseq(snp_concordance_channel, snp_region_file).vcf_files
    ch_versions = ch_versions.mix(mpileup_strandseq.out.versions.take(1))

    calc_discordance(mpileup_vcf_files)
    ch_versions = ch_versions.mix(calc_discordance.out.versions.take(1))

    // Filter samples that SNP concordance check
    snp_concordant_samples = calc_discordance.out.snp_discordance
        .filter { patient, discordance_file ->
            return check_snp_discordance(patient, discordance_file)
        }
        .map { patient, discordance_file -> patient }

    strandseq_bam_folder = strandseq.out.pass_strandseq_bam_folder.join(snp_concordant_samples, by: 0)
    // join to get snp concordant only

    //ashley's qc - scores the cells, returns folders with the cells that pass and fail
    pass_bams = ashleys_qc(strandseq_bam_folder, params.min_ashleys_qc_prob).ashleys_pass_folder
    ch_versions = ch_versions.mix(ashleys_qc.out.versions.take(1))

    //strandseq_phase.R
    strandseq_phase_channel = nanopore_vcf
        .join(pass_bams)
        .join(strandseq.out.is_paired_end)
    strandseq_phase(strandseq_phase_channel, params.chromosomes)
    ch_versions = ch_versions.mix(strandseq_phase.out.versions.take(1))

    //chromosome phasing summary
    write_chr_phasing_summary(strandseq_phase.out.chr_states)
    ch_versions = ch_versions.mix(write_chr_phasing_summary.out.versions.take(1))

    //het-hom-discordance
    calc_het_hom_discordance(nanopore_vcf.join(strandseq_phase.out.vcf))
    ch_versions = ch_versions.mix(calc_het_hom_discordance.out.versions.take(1))

    // Filter samples that failed het-hom check
    clean_samples = calc_het_hom_discordance.out.het_hom_discordance
        .filter { patient, results ->
            return check_het_hom_discordance(patient, results)
        }
        .map { patient, results -> patient }

    filtered_vcfs = strandseq_phase.out.vcf.join(clean_samples, by: 0)
    // join by 1st element, get only clean samples.

    // compress strandseq VCF
    compress_strandseq_phase_vcf(filtered_vcfs, "strandseqphase")

    //patmat
    patmat_channel = nanopore_vcf
        .join(input.nanopore_sv_vcf)
        .join(input.nanopore_bam)
        .join(filtered_vcfs)
    pofo_assignment_vcf = patmat(
        patmat_channel,
        params.reference_genome,
        params.reference_genome_idx,
        channel.fromPath(params.idmr_regions_file, checkIfExists: true),
    ).pofo_assignment_vcf
    ch_versions = ch_versions.mix(patmat.out.versions.take(1))

    vcf_to_compress = pofo_assignment_vcf.mix(patmat.out.pofo_assignment_sv_vcf)
    compress_patmat_vcf(vcf_to_compress, "POAGA")
    ch_versions = ch_versions.mix(compress_patmat_vcf.out.versions.take(1))

    //Haplotagging the BAM file with whatshap
    prepare_whatshap_channel = input.nanopore_bam.join(pofo_assignment_vcf)
    prep_files = prepare_whatshap_files(prepare_whatshap_channel).prep_files
    ch_versions = ch_versions.mix(prepare_whatshap_files.out.versions.take(1))

    haplotagged_bam = whatshap_haplotag(prep_files, params.reference_genome, params.reference_genome_idx).haplotagged_bam
    ch_versions = ch_versions.mix(whatshap_haplotag.out.versions.take(1))

    sort_haplotagged_bam(haplotagged_bam)
    ch_versions = ch_versions.mix(sort_haplotagged_bam.out.versions.take(1))

    //Compress BAM to CRAM
    compress_haplotagged_bam(sort_haplotagged_bam.out.haplotagged_sorted_bam, params.reference_genome)
    ch_versions = ch_versions.mix(sort_haplotagged_bam.out.versions.take(1))

    // Variants of interest:Not
    if (params.variants) {
        extract_variants_of_interest(
            params.variants,
            patmat.out.pofo_assignment_vcf,
            patmat.out.pofo_scores,
            patmat.out.pofo_dmr_status,
            channel.fromPath(params.idmr_regions_file, checkIfExists: true),
            sort_haplotagged_bam.out.haplotagged_sorted_bam,
        )
        ch_versions = ch_versions.mix(extract_variants_of_interest.out.ch_versions)
    }

    // Merge and publish logs
    ch_versions.unique().collectFile(name: params.logfile, storeDir: params.output_folder)
}

/*
Takes care of the preliminary strand-seq analysis including:
-making a single pair of fastq files per cell
-align with bowtie2
-sort, dupmark, etc.
-Filter out cells with low fragment counts

returns:
-folders for cells that pass and fail
-lists of cells that pass and fail
*/
workflow strandseq {
    take:
    strandseq_reads

    main:
    ch_versions = Channel.empty()

    //will process each fastq pair as its own job
    strandseq_cutadapt(strandseq_reads)
    ch_versions = ch_versions.mix(strandseq_cutadapt.out.versions.take(1))

    // build is_paired channel, per fastq per patient.
    is_paired_fastq = strandseq_cutadapt.out.fastqs.map { patient, fastq, is_paired -> tuple(patient, is_paired) }

    // log PE vs SE counts:
    print_pe_vs_counts(is_paired_fastq)

    // Determine full SE/PE vs mixed:
    is_mixed_ch = determine_if_analysis_is_mixed(is_paired_fastq)

    // Combine is_mixed channel with fastqs, to check if we need "conversion" from PE to SE
    strandseq_cutadapt.out.fastqs
        .combine(is_mixed_ch, by: 0)
        .branch {
            needs_conversion: it[2] == true && it[3] == true
            already_single: it[2] == false
            all_paired: it[2] == true && it[3] == false
        }
        .set { fastq_branches }

    // Convert paired files to single when needed
    convert_fastq_pair_to_single(fastq_branches.needs_conversion)

    // remove 'is_mixed' from non processing channel and merge all channels back together
    already_single_standardized = fastq_branches.already_single.map { patient, fastq, is_paired, is_mixed -> tuple(patient, fastq, is_paired) }
    all_paired_standardized = fastq_branches.all_paired.map { patient, fastq, is_paired, is_mixed -> tuple(patient, fastq, is_paired) }

    // Merge all channels back together
    already_single_standardized
        .mix(all_paired_standardized)
        .mix(convert_fastq_pair_to_single.out.fastq)
        .set { unified_fastqs_ch }

    //bowtie2 - outputs sam
    strandseq_bowtie(unified_fastqs_ch, params.bowtie_index)
    ch_versions = ch_versions.mix(strandseq_bowtie.out.versions.take(1))

    // Channel with paired-end boolean per patient, to pass downstream
    // If data had mixed SE/PE, each PE fastq pair was "merged" into a single SE fastq.
    is_paired_end = unified_fastqs_ch
        .map { tuple(it[0], it[2]) }
        .unique()
    // should be tuple of (patient, true/false).


    //tidybam
    strandseq_tidy_reads(strandseq_bowtie.out.sam_file)
    ch_versions = ch_versions.mix(strandseq_tidy_reads.out.versions.take(1))

    //remove bams with very few reads
    strandseq_mark_duplicates(strandseq_tidy_reads.out.bam_file)
    ch_versions = ch_versions.mix(strandseq_mark_duplicates.out.versions.take(1))

    //split bams into pass/fail sets
    strand_pass_fail_cells(strandseq_mark_duplicates.out.bam_file.groupTuple(), params.min_reads_threshold)
    ch_versions = ch_versions.mix(strand_pass_fail_cells.out.versions.take(1))

    emit:
    passed_strand_seq_list = strand_pass_fail_cells.out.passed_strand_seq_list
    failed_strand_seq_list = strand_pass_fail_cells.out.failed_strand_seq_list
    pass_strandseq_bam_folder = strand_pass_fail_cells.out.pass_strandseq_bams
    fail_strandseq_bam_folder = strand_pass_fail_cells.out.fail_strandseq_bams
    is_paired_end = is_paired_end
    ch_versions
}

def determine_if_analysis_is_mixed(is_paired_fastq) {
    def is_mixed_ch = is_paired_fastq
        .groupTuple()
        .map { patient, is_paired_list ->
            def paired_count = is_paired_list.count(true)
            def single_count = is_paired_list.count(false)
            def is_mixed = (paired_count > 0 && single_count > 0)
            return tuple(patient, is_mixed)
        }
    // Log if mixed:
    is_mixed_ch
        .filter { it[1] == true }
        .view { "INFO: Analysis for patient ${it[0]} has mixed SE/PE libraries, treating all as SE." }

    return is_mixed_ch
}

def print_pe_vs_counts(is_paired_fastq) {
    is_paired_fastq
        .groupTuple()
        .map { patient, is_paired_end_list ->
            def pe_count = is_paired_end_list.count(true)
            def se_count = is_paired_end_list.count(false)
            tuple(patient, pe_count, se_count)
        }
        .view { patient, pe_count, se_count ->
            "Patient: ${patient}, Paired End libraries: ${pe_count}, Single End libraries: ${se_count}"
        }
}

def extract_column(file, column_name) {
    // Read file and skip header
    def lines = file.readLines()

    // Get header to find Discordance column index
    def header = lines[0].split('\t')
    def valueColIdx = header.findIndexOf { it == column_name }

    if (valueColIdx == -1) {
        log.info("ERROR: No '${column_name}' column found in file ${file}")
        return -1
    }
    // Get value from second line
    return lines[1].split('\t')[valueColIdx].toFloat()
}

def check_het_hom_discordance(patient, results) {
    def value = extract_column(results, "Discordance")
    log.info("Patient ${patient} HET-HOM discordance: ${value}")

    // Compare with threshold (check for -1, meaning did not parse the file correctly)
    def passes = (value != -1 && value <= params.het_hom_threshold)

    // Log filtered samples
    if (!passes) {
        log.info("Patient ${patient} filtered out due to HET-HOM discordance value: ${value}")
    }
    return passes
}

def check_snp_discordance(patient, discordance_file) {
    def value = extract_column(discordance_file, "DiscordantP")
    log.info("Patient ${patient} SNP discordance: ${value}")

    // Compare with threshold (check for -1, meaning did not parse the file correctly)
    def passes = (value != -1 && value <= params.snp_discordance_threshold)

    // Log filtered samples
    if (!passes) {
        log.info("Patient ${patient} filtered out due to SNP discordance value: ${value}")
    }
    return passes
}

def build_strandseq_fastq_input(strandseq_input) {
    // Given the strandseq input with tuples (patient, fastq_folder, fastq_list), find all cell fastq files
    // and return a tuple (patient, cell_name, r1_file, r2_file, is_paired) where:
    // - cell name: parsed cell name from R1 file
    // - r1_file and r2_file: R1/R2 fastQ files, with R2 possibly missing, for SE libraries.
    // - is_paired: boolean indicating if R2 file exists, which means a PE library.
    // The fastq files are sourced from either a folder or a ';'-separated list of files, whichever is provided.
    return strandseq_input.flatMap { patient, fastq_folder, fastq_list ->
        def r1_files
        if (fastq_list?.toString()?.trim()) {
            // A ';'-separated list of individual fastq file paths was provided.
            // This branch takes precedence when both inputs are present.
            def all_files = fastq_list.split(';').collect { it.trim() }.findAll { it }.collect { file(it) }
            r1_files = all_files.findAll { it.name ==~ /.*_R1.*\.fastq\.gz$/ }
        }
        else if (fastq_folder?.toString()?.trim()) {
            // A folder was provided; glob all R1 files inside it.
            r1_files = file("${fastq_folder}/*_R1*.fastq.gz")
        }
        else {
            error("No Strand-seq FASTQ input found for patient ${patient}. Provide either strand_seq_fastq or strand_seq_fastq_list.")
        }

        // Create a list of pairs, handling optional R2
        return r1_files.collect { r1_file ->
            // Get the cell id:
            def cell_id = r1_file.name.toString().replaceAll(/_R1.+$/, '')
            // Get the matching R2 file
            def r2_file = file(r1_file.toString().replace('_R1', '_R2'))
            // Check if R2 exists to determine if it is PE or SE:
            def is_paired = r2_file.exists()

            return tuple(patient, cell_id, r1_file, r2_file, is_paired)
        }
    }
}
