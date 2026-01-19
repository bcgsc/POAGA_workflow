# POAGA Pipeline - Developer Guide

This guide is intended for developers who will be maintaining and updating the Parent of Origin-Aware Genomic Analysis (POAGA) pipeline. It provides an architectural overview of the pipeline, key components, customization points, and best practices.

## Pipeline Architecture

The POAGA pipeline is built using Nextflow DSL2 and consists of several modules that perform the analysis of strand-seq and nanopore data to determine parent-of-origin information. The pipeline architecture follows a modular approach with the following key components:

### Main Workflow Files

- `poaga_processes.nf`: Parent workflow that includes all processes and defines the main workflow
- `extract_variants.nf`: Supplementary workflow for extracting variants of interest
- `processes/` directory: Contains individual process files organized by functionality

### Modular Process Organization

Processes are grouped into modules based on functionality:

- `processes/compress.nf`: File compression and manipulation
- `processes/qc.nf`: Quality control processes
- `processes/regions_of_interest.nf`: Handling regions of interest
- `processes/patmat.nf`: Parent-of-origin assignment
- `processes/whatshap.nf`: Phasing/haplotagging
- `processes/strandseq.nf`: Strand-seq processing
- `processes/snp_concordance.nf`: SNP concordance checks

## Key Pipeline Components and Flow

The pipeline follows this general workflow:

1. **Input Processing**:
   - Parses input CSV with patient, strand-seq, and nanopore data locations
   - Builds channel structures for downstream processing

2. **StrandSeq Processing**:
   - Cutadapt trimming
   - Paired/single-end processing logic
   - Bowtie2 alignment
   - Read filtering/tidying
   - Duplicate marking
   - QC-based filtering

3. **QC Steps**:
   - SNP concordance checks
   - Ashley's QC filtering
   - Het/hom discordance filtering

4. **Phasing and Parent-of-Origin Assignment**:
   - StrandSeq phasing
   - PatMat tool processing
   - Haplotagging with WhatsHap

5. **Output Generation**:
   - Compressed BAM/CRAM files
   - VCF generation and compression
   - Variant of interest extraction (optional)

## Channel Design and Management

The pipeline uses Nextflow channels extensively for data flow control. Key channel design patterns include:

### Multi-Map Pattern for the Input CSV

Used to split the input CSV file into multiple channels with the FastQ and VCF files:

```groovy
input = Channel
    .fromPath(params.input_csv)
    .splitCsv(header:true)
    .multiMap { row ->
        strandseq: tuple(row.patient, row.strand_seq_fastq)
        nanopore_vcf: tuple(row.patient, row.clair3_vcf)
        nanopore_sv_vcf: tuple(row.patient, row.sniffles2_vcf)
        nanopore_bam: tuple(row.patient, row.bam, row.bai)
    }
```

### Join/Combine Patterns

Used to synchronize multiple parallel processes:

```groovy
patmat_channel = nanopore_vcf
    .join(input.nanopore_sv_vcf)
    .join(input.nanopore_bam)
    .join(filtered_vcfs)
```

### Filtering and Transformation

Used to conditionally process samples based on QC metrics:

```groovy
clean_samples = calc_het_hom_discordance.out.het_hom_discordance
    .filter { patient, results ->
        return check_het_hom_discordance(patient, results)
    }.map { patient, results -> patient }
```

## Helper Functions

The pipeline contains several helper functions to manage workflow logic:

- `determine_if_analysis_is_mixed()`: Detects mixed single/paired-end data
- `print_pe_vs_counts()`: Logs paired vs single-end library counts
- `check_het_hom_discordance()`: Evaluates heterozygous vs homozygous discordance metric against a given threshold
- `check_snp_discordance()`: Evaluates SNP discordance against a given threshold
- `build_strandseq_fastq_input()`: Builds standardized strandseq input channels

## Configuration System

The pipeline configuration is managed through:

### Core Parameters

The main configurable parameters are defined with defaults:

```groovy
// Mandatory
params.input_csv = null
params.output_folder = null

// Optional with defaults
params.min_reads_threshold = 10000
params.chromosomes = "chr1,chr2,...,chr22"
params.min_ashleys_qc_prob = 0.5
```

### Resource Profiles

Process resources are managed with labels:

```groovy
withLabel: mem256 {
  memory = {256.GB * task.attempt}
  maxRetries = 2
}
```

### Reference Profiles

Genome references are managed through profiles:

```groovy
profiles {
  hg38_no_alt {
    params.bowtie_index = "/projects/poaga/reference/"
    params.reference_genome = "/projects/poaga/reference/hg38_no_alt.fa"
    params.reference_genome_idx = "/projects/poaga/reference/hg38_no_alt.fa.fai"
  }
}
```

## Container Management

The pipeline uses Singularity containers for reproducibility:

- Each process specifies a container, e.g., `container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"`
- The `singularity` section in the config controls mount behavior and cache location
- Custom containers are stored in `params.singularity_folder`

## Error Handling and Retry Logic

The pipeline implements robust error handling:

- Default `errorStrategy = 'retry'` for all processes
- Memory-based retries with increasing allocation: `memory = {64.GB * task.attempt}`
- Maximum retry limits set per resource class: `maxRetries = 3`

## Adding New Functionality

When extending the pipeline, follow these guidelines:

### Adding a New Process

1. Create a new process in the appropriate process file based on functionality
2. Define inputs, outputs, and container requirements
3. Include version logging in the script section
4. Emit version information for central version tracking

Example template:

```groovy
process new_process {
    tag "${patient}"
    cpus 4
    label "mem12"
    container "container:version"
    publishDir "${params.output_folder}/${patient}/output_dir", mode: 'copy'

    input:
        tuple val(patient), path(input_file)

    output:
        tuple val(patient), path(output_file), emit: result
        path("versions.yml"), emit: versions

    script:
    """
        command_to_run ${input_file} > ${output_file}

        # Log tool versions and containers
        cat <<-END_VERSIONS > versions.yml
        ---
        Process:"${task.process}"
        Container: "${task.container}"
        tool: \$( tool --version )
        END_VERSIONS
    """
}
```

### Integrating a New Process

1. Include the process in the main workflow file:
   ```groovy
   include { new_process } from './processes/module'
   ```

2. Add the process to the workflow:
   ```groovy
   new_process(input_channel)
   ch_versions = ch_versions.mix(new_process.out.versions.take(1))
   ```

3. Connect the process outputs to downstream processes through channels

### Adding Parameters

1. Add parameter in the main config with a sensible default
2. Document the parameter in the help message
3. Use the parameter in processes with `params.parameter_name`

## Testing

The pipeline includes a `test` profile for testing that runs much faster than a full sample.

```groovy
test {
  params.input_csv = "example_minimal_samples.csv"
  params.min_reads_threshold = 500
  params.chromosomes = "chr19"
  params.min_ashleys_qc_prob = 0
}
```

## Version Management

The pipeline tracks tool versions using the a common Nextflow tracking approach:

1. Each process generates a `versions.yml` file
2. The main workflow collects and merges these files
3. The result is published to `params.logfile`

