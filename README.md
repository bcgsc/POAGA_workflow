# POAGA: Parent of Origin-Aware Genomic Analysis Pipeline

## Overview

POAGA is a Nextflow-based bioinformatics pipeline that combines Strand-seq and Oxford Nanopore Technologies (ONT) data to determine parent-of-origin information for genetic variants. The pipeline processes Strand-seq data for chromosomal phasing and integrates this with Nanopore long-read data to assign maternal and paternal origin to variants.

## Features

- Processing and QC of Strand-seq libraries
- Integration of Nanopore sequencing data
- SNP concordance analysis between short and long reads
- Chromosomal phasing using Strand-seq data
- Parent-of-origin assignment for SNVs and SVs
- BAM haplotagging
- Analysis of methylation patterns in imprinted regions
- Extraction of variants of interest

## Requirements

- Nextflow (≥23.04.4)
- Singularity/Apptainer
- Access to a SLURM-based HPC environment (or modify executor settings for your environment)
- Reference genome files

## Input Requirements

### Main Input CSV

The pipeline requires an input CSV file with the following columns:

- `patient`: Sample label/identifier
- `strand_seq_fastq`: Path to folder with Strand-seq FastQ files
- `clair3_vcf`: Path to Clair3 VCF output for Nanopore data
- `sniffles2_vcf`: Path to Sniffles2 SV VCF output for Nanopore data
- `bam`: Path to Nanopore BAM file
- `bai`: Path to index file for the Nanopore BAM

Example input CSV format:
```
patient,strand_seq_fastq,clair3_vcf,sniffles2_vcf,bam,bai
sample1,/path/to/strandseq/sample1,/path/to/clair3/sample1.vcf.gz,/path/to/sniffles2/sample1.vcf.gz,/path/to/nanopore/sample1.bam,/path/to/nanopore/sample1.bam.bai
sample2,/path/to/strandseq/sample2,/path/to/clair3/sample2.vcf.gz,/path/to/sniffles2/sample2.vcf.gz,/path/to/nanopore/sample2.bam,/path/to/nanopore/sample2.bam.bai
```

**IMPORTANT**: The StrandSeq read names must have the format \*_*R{1,2}*.fastq.gz. Single-ended libraries are
also supported, meaning that the R2 files are optional. Also, the  minimum fragment count for usable strand-seq cells is 5000.

### Variants of Interest CSV input (Optional)

When using the `--variants` parameter, you need to provide a CSV file with the following columns:

- `patient`: Sample label/identifier (must match those in the main input CSV)
- `variant`: Variant in HGVS cDNA notation (e.g., `NM_000492:c.1521_1523delCTT`)

Example variants CSV format:
```
patient,variant
sample1,NM_000492:c.1521_1523delCTT
sample1,NM_000492:c.3909C>G
sample2,NM_000059:c.6275_6276delTT
```

The pipeline also requires a list of imprinted DMR (Differentially Methylated Regions), configured with `--idmr_regions_file` (default: `${projectDir}/assets/Imprinted_DMR_List_V1.tsv`).

## Usage

### Basic Usage

```bash
nextflow run poaga_workflow.nf \
  --input_csv /path/to/samples.csv \
  --output_folder /path/to/results \
  -profile hg38_no_alt
```

### With Variants of Interest

```bash
nextflow run poaga_workflow.nf \
  --input_csv /path/to/samples.csv \
  --output_folder /path/to/results \
  --variants /path/to/variants_of_interest.csv \
  --flanking_region_size 20000 \
  -profile hg38_no_alt
```

The above example includes a larger flanking region (20kb) around each variant of interest.

### Using a subset of Chromosomes

```bash
nextflow run poaga_workflow.nf \
  --input_csv /path/to/samples.csv \
  --output_folder /path/to/results \
  --chromosomes "chr1,chr2,chr3,chr4,chr5" \
  -profile hg38_no_alt
```

### With custom QC Parameters

```bash
nextflow run poaga_workflow.nf \
  --input_csv /path/to/samples.csv \
  --output_folder /path/to/results \
  --min_reads_threshold 15000 \
  --min_ashleys_qc_prob 0.6 \
  --het_hom_threshold 0.02 \
  --snp_discordance_threshold 0.02 \
  -profile hg38_no_alt
```

### Running with Clinical Priority (SLURM priority)

```bash
nextflow run poaga_workflow.nf \
  --input_csv /path/to/samples.csv \
  --output_folder /path/to/results \
  -profile hg38_no_alt,clinical_priority
```

### Running in Test Mode
The `test` profile loads a "minimal example" dataset that should run in about 15-20 minutes.
```bash
nextflow run poaga_workflow.nf \
  --output_folder ./test_output \
  -profile hg38_no_alt,test
```

## Parameter description

### Required Parameters

- `--input_csv`: Location of the input CSV file
- `--output_folder`: Output folder to store results

### Optional Parameters

- `--variants`: CSV file(s) with variants of interest in HGVS cDNA notation
- `--idmr_regions_file`: Path to iDMR list TSV used during variants-of-interest extraction (default: `${projectDir}/assets/Imprinted_DMR_List_V1.tsv`)
- `--min_reads_threshold`: Minimum number of high-quality reads to pass a StrandSeq library (default: 10000)
- `--chromosomes`: Comma-separated list of chromosomes to process (default: all autosomes)
- `--min_ashleys_qc_prob`: Minimum probability score for Strand-seq library QC (default: 0.5)
- `--het_hom_threshold`: Maximum allowed het-hom discordance value (default: 0.03)
- `--snp_discordance_threshold`: Maximum allowed SNP discordance value (default: 0.03)
- `--snp_region`: Region/BED file for SNP concordance analysis (default: `${projectDir}/assets/snp_regions.bed`)
- `--flanking_region_size`: Size of flanking region for variant analysis in base pairs (default: 10000)
- `--extra_strandseq_options`: Additional options for Strand-seq processing

### Profiles

- `hg38_no_alt`: Uses hg38 reference genome without alternative contigs
- `clinical_priority`: Uses higher priority in job scheduler for clinical samples
- `test`: Runs minimal dataset for testing/development

## Pipeline Stages

1. **Strand-seq processing**
   - Adapter trimming (cutadapt)
   - Alignment (Bowtie2)
   - Read filtering and QC

2. **SNP concordance**
   - Mpileup of Strand-seq data
   - Discordance calculation with Nanopore data
   - Sample filtering based on discordance threshold

3. **Strand-seq QC**
   - ASHLEYS QC for cell scoring
   - Filtering of low-quality cells

4. **Phasing**
   - Strand-seq phasing
   - Chromosome state determination
   - Het-hom discordance calculation

5. **Parent-of-origin assignment**
   - Patmat processing for SNVs and SVs
   - VCF generation with parent-of-origin information

6. **Haplotagging**
   - WhatsHap processing
   - BAM sorting and compression
   - CRAM generation

7. **Variant analysis** (if specified)
   - Annotation of variants using TransVar
   - Extraction of variants from POAGA VCF files
   - Generation of regions of interest BED files around variants
   - Subsetting of BAM files to regions of interest
   - Analysis of methylation patterns in imprinted regions

## Output Files

The pipeline generates multiple output files in the specified output folder:
- Phased VCF files with parent-of-origin information
- Haplotagged BAM/CRAM files
- QC reports for Strand-seq data
- Chromosome phasing summaries
- Tool version logs

When using the `--variants` parameter, the pipeline creates a `variants_of_interest` folder for each patient containing
subsets of the main output files targeting the given variants.

### POAGA Pipeline Output Structure

This section details the final output structure of the POAGA pipeline, organized by directories and files.

#### Output Directory Structure

For each patient processed by the pipeline, a directory is created under the specified `--output_folder` with the following structure:

```
output_folder/
└── <patient>/
    ├── ashleys_qc/
    ├── POAGA/
    │   └── variants_of_interest/
    ├── qc_metrics/
    └── strandseqphase/
```

Below is a detailed description of the files contained in each directory.

#### 1. POAGA Directory

The main output directory containing the core results.

##### Core Output Files:

- `PofO_HaploTagged.sorted.cram` + `.crai`: Compressed haplotagged BAM file and its index
- `patmat_PofO_Scores.tsv`: Parent-of-origin methylation scores
- `patmat_PofO_Assigned_info.tsv`: Information about parent-of-origin assignments
- `*_PofO_Assignment_SVs_info.tsv`: Information about structural variant parent-of-origin assignments
- `patmat_CpG-Methylation-Status-at-DMRs.tsv`: Methylation status at differentially methylated regions

##### VCF Files (when processed by compress_patmat_vcf):

- `*.POAGA.vcf.gz` + `.tbi`: Compressed and indexed VCF file with parent-of-origin assignments

#### 2. variants_of_interest Directory

Created when `--variants` parameter is provided. Contains analysis of specific variants of interest.

##### Files:

- `*.annotation.txt`: TransVar annotations for each variant
- `variants_of_interest.vcf`: Subset of VCF containing only variants of interest
- `variants_of_interest_status.txt`: Status report of each variant of interest
- `regions_of_interest.bed`: BED file of regions around variants with flanking regions
- `regions_of_interest.bam` + `.bai`: Subset of haplotagged BAM containing only reads in regions of interest
- `pofo_scores_regions_of_interest.tsv`: Subset of parent-of-origin scores for regions of interest
- `CpG-methylation-status-at-DMRs-regions_of_interest.tsv`: Methylation status at DMRs in regions of interest

#### 3. strandseqphase Directory

Contains outputs from the Strand-seq phasing process.

##### Files:

- `<patient>.chr-states.txt`: Chromosome state information from Strand-seq phasing
- `<patient>.chr_phasing_summary.txt`: Summary of chromosome phasing results
- `<patient>_WW_composite.bam` + `.bai`: Composite BAM of Watson-Watson strand states
- `<patient>_WC_composite.bam` + `.bai`: Composite BAM of Watson-Crick strand states
- `invertyper/`: Directory with inversion typing results

##### Additional Files (when processed by compress_strandseq_phase_vcf):

- `*.strandseqphase.vcf.gz` + `.tbi`: Compressed and indexed phased VCF from Strand-seq

#### 4. ashleys_qc Directory

Contains quality control results from ASHLEYS.

##### Files:

- `pass/`: Directory containing BAM files that passed ASHLEYS QC
- `fail/`: Directory containing BAM files that failed ASHLEYS QC
- `quality.txt`: QC scores for each cell

#### 5. qc_metrics Directory

Contains quality metrics for the analysis.

##### Files:

- `het_hom_discordance.txt`: Heterozygous-homozygous discordance metrics

#### 6. root_directory File

The pipeline also outputs a version log file in the root output directory:

- `tool_versions.yml`: Version information for all tools used in the pipeline


#### Debug Mode

When the pipeline is run with `params.debug=true`, additional intermediate files may be published:

- Complete `patmat` output in the POAGA directory
- Complete `strandseq_phase` output in the strandseqphase directory

## Troubleshooting and Other details

### Common Issues

1. **Insufficient resources**: The patmat and StrandSeq phasing steps can use a lot of memory, exceeding the memory specification which can cause them to fail. A common workaround is running specific high-memory cases in a high memory `gphost`  machine.

2. **Failed QC**: Check QC logs and adjust thresholds if necessary. For instance, some samples might require a higher threshold for some of the concordance checks.
   ```bash
   nextflow run poaga_workflow.nf -resume --input_csv samples.csv --output_folder results --het_hom_threshold 0.05
   ```

### Resuming Failed Runs

To resume a failed run add the `-resume` flag:
```bash
nextflow run poaga_workflow.nf -resume --input_csv samples.csv --output_folder results
```

### Variant Analysis Workflow details

The variant extraction workflow consists of several steps:

1. **Variant Annotation**: Uses TransVar to annotate variants specified in HGVS notation
   ```bash
   # Example of a standalone TransVar command (for reference)
   transvar canno --refseq -i 'NM_000492:c.1521_1523delCTT' --longest
   ```

2. **Variant Extraction**: Identifies the variants in the parent-of-origin-assigned VCF files

3. **Region of Interest Generation**: Creates BED files for regions around variants, considering:
   - The variant position
   - A flanking region (controlled by `--flanking_region_size`)
   - Nearby imprinted regions from the DMR list

4. **BAM Subsetting**: Extracts reads from the specified regions for focused analysis
   ```bash
   # Example of a standalone samtools command (for reference)
   samtools view -b -h -L regions_of_interest.bed input.bam > regions_of_interest.bam
   ```

5. **Methylation Analysis**: Focuses methylation analysis on the regions of interest

### Custom Reference Genome

To use a custom reference genome, modify the `nextflow.config` file or create a new profile:

```
profiles {
  custom_reference {
    params.bowtie_index = "/path/to/custom/bowtie/index"
    params.reference_genome = "/path/to/custom/reference.fa"
    params.reference_genome_idx = "/path/to/custom/reference.fa.fai"
  }
}
```

## Maintenance and Development

This pipeline is designed to be maintainable and extensible. If you're interested in contributing to the development, understanding the pipeline architecture, or extending its functionality, please refer to our [Development Documentation](./DEVELOPMENT.md).

The developer documentation includes:

- Pipeline architecture and workflow overview
- Process organization and data flow
- Guide for adding new functionality
- Configuration system details
