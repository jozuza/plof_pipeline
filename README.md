# LoF Pipeline (PLOF) - Comprehensive Documentation

**Author:** JCMM  
**Date:** February 2026  
**Purpose:** End-to-end variant annotation pipeline for identifying and annotating Loss-of-Function (LoF) variants using VEP and LoFTEE

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Prerequisites & Dependencies](#prerequisites--dependencies)
4. [Installation & Setup](#installation--setup)
5. [Directory Structure](#directory-structure)
6. [Configuration](#configuration)
7. [Running the Pipeline](#running-the-pipeline)
8. [Detailed Module Documentation](#detailed-module-documentation)
9. [Input/Output Formats](#inputoutput-formats)
10. [Environment Variables](#environment-variables)
11. [Troubleshooting](#troubleshooting)
12. [Performance Considerations](#performance-considerations)

---

## Overview

The **PLOF Pipeline** is a comprehensive variant analysis workflow designed to:

- **Download** VCF files from Google Cloud Storage (GCS)
- **Filter** variants with PASS status and quality control criteria
- **Normalize** variants and reduce to relevant genomic regions
- **Annotate** variants for functional impact using Ensembl VEP and LoFTEE
- **Parse** results into structured tabular format for downstream analysis

This pipeline is optimized for identifying **Loss-of-Function variants** that are likely to have high functional impact on protein-coding genes.

### Key Features

- ✅ Supports both **site-only** and **cohort** variant datasets
- ✅ Automatic quality control with Hardy-Weinberg equilibrium testing
- ✅ Integration with Google Cloud Storage for large-scale data access
- ✅ Chromosome-level parallelization for efficiency
- ✅ Detailed logging and error handling
- ✅ Intermediate file caching to skip completed steps

---

## Pipeline Architecture

The pipeline follows a modular, sequential design with these main stages:

```
INPUT VCF (from GCS)
       ↓
[STEP 0] Download
       ↓
[STEP 1] PASS Filter
       ↓
[STEP 2A] QC Filter (HWE, Missingness) [cohort mode only]
       ↓
[STEP 2B] Normalization + Genotype Reduction + Protein Coding Filter
       ↓
[STEP 3] VEP + LoFTEE Annotation
       ↓
OUTPUT VCF (Annotated)
       ↓
[OPTIONAL] Parse to TSV
```

### Execution Model

- **Chromosome-based processing**: Each chromosome is processed independently in a serial loop (can be parallelized externally)
- **Step caching**: If output files exist, steps are automatically skipped
- **Comprehensive logging**: All operations are logged to `<BASENAME>.<CHR>.pipeline.log`

---

## Prerequisites & Dependencies

### System Requirements

- **OS**: Linux (tested on Ubuntu 20.04+)
- **Disk Space**: Minimum 500 GB for VCF files and temporary processing
- **RAM**: 32–64 GB recommended
- **Threads/CPU**: 7+ cores (configurable via `THREADS` variable)

### Required Software

1. **bcftools** (>=1.14) - VCF manipulation and filtering
2. **vcftools** (>=1.16) - VCF analysis and QC
3. **Ensembl VEP** (>=106) - Variant Effect Predictor
4. **LoFTEE plugin** - LoF annotation for VEP
5. **Google Cloud SDK** (`gsutil`) - Cloud file access
6. **Reference files**:
   - Human genome FASTA (GRCh38/hg38)
   - VEP cache (offline mode)
   - LoFTEE plugin files
   - Protein coding genes BED file

### Optional

- **Conda/Conda-pack**: For environment reproducibility

---

## Installation & Setup

### 1. Clone/Download the Repository

```bash
git clone <repository-url>
cd plof_pipeline
```

### 2. Prepare Directories

```bash
mkdir -p /data/vcf_files /data/temp /data/results/vep_loftee /data/references/human_genome
```

### 3. Download Reference Files

#### Human Genome (GRCh38)

```bash
# Download Ensembl GRCh38 reference
wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    -O /data/references/human_genome/Homo_sapiens_assembly38.fasta.gz

gunzip /data/references/human_genome/Homo_sapiens_assembly38.fasta.gz
samtools faidx /data/references/human_genome/Homo_sapiens_assembly38.fasta
```

#### Protein Coding Genes (BED)

```bash
# From Ensembl or GENCODE
wget https://www.gencodegenes.org/releases/current.html  # Download GRCh38 comprehensive gene annotation

# Extract protein-coding gene regions
zcat gencode.v43.basic.annotation.gff3.gz | \
    awk '$3 == "exon" && $9 ~/gene_type "protein_coding"/' | \
    awk '{print $1"\t"$4-1"\t"$5}' | \
    sort -u > /data/references/human_genome/coding_genes.bed
```

#### VEP Cache & LoFTEE Plugin

```bash
# Download VEP cache (requires space - ~10GB)
mkdir -p /opt/vep/vep_cache
# Use VEP's installer or manually download from Ensembl

# Download LoFTEE plugin
mkdir -p /opt/vep_plugins/loftee
wget https://github.com/konradjk/loftee/archive/refs/heads/master.zip -O /tmp/loftee.zip
unzip /tmp/loftee.zip -d /opt/vep_plugins/
mv /opt/vep_plugins/loftee-master /opt/vep_plugins/loftee
```

### 4. Install System Dependencies

```bash
apt-get update && apt-get install -y \
    bcftools \
    vcftools \
    samtools \
    curl

./scripts/prep_cloud.sh  # Install Google Cloud SDK
```

### 5. Verify Installation

```bash
bcftools --version
vcftools --version
gsutil --version
```

---

## Directory Structure

```
plof_pipeline/
├── README.md                          # This file
├── scripts/
│   ├── plof_end_to_end.sh            # Main pipeline controller
│   ├── prep_cloud.sh                 # Google Cloud SDK setup
│   ├── pass_filter.sh                # PASS filter module
│   ├── filterQC.sh                   # QC filter module
│   ├── normalize.sh                  # Normalization module
│   ├── vep_loftee.sh                 # VEP+LoFTEE annotation module
│   └── parse.sh                      # TSV parsing module (optional)
│
├── /data/
│   ├── vcf_files/                    # Input VCF files (from GCS or local)
│   ├── temp/                         # Temporary processing files
│   ├── results/
│   │   └── vep_loftee/              # Final annotated VCFs and logs
│   └── references/
│       └── human_genome/
│           ├── Homo_sapiens_assembly38.fasta
│           ├── Homo_sapiens_assembly38.fasta.fai
│           └── coding_genes.bed
```

---

## Configuration

All user-configurable parameters are defined in [plof_end_to_end.sh](scripts/plof_end_to_end.sh) under the **USER CONFIGURATION** section.

### Key Configuration Parameters

```bash
# Data Type
export INPUTTYPE="site_only"          # Options: "site_only" | "cohort"
                                      # site_only = no genotype data
                                      # cohort = contains genotype samples

# File Naming
export VCFBASENAME="bravo.pub"        # Base name for all output files
export EXTENSIONNAME="vcf.gz"         # VCF file extension

# Cloud Storage (GCS)
export BUCKET="gs://projects-usp/public-datasets/TOPMed/references"
                                      # GCS bucket containing input VCFs

# Computational Resources
export THREADS=7                      # Number of parallel threads for tools
                                      # Recommend: (CPU cores - 1)
```

### Fixed Paths

```bash
export INPUTDIR="/data/vcf_files"              # Downloaded VCF input
export TMPDIR="/data/temp"                     # Temporary working files
export OUTPUTDIR="/data/results/vep_loftee"    # Final output directory
export REFFASTA="/data/references/human_genome/Homo_sapiens_assembly38.fasta"
export REFCODING="/data/references/human_genome/coding_genes.bed"
```

### Chromosome Range

Currently set to chromosomes 22 and 21 for testing:

```bash
for chr in {22..21}; do
```

**To process all chromosomes**, change to:

```bash
for chr in {1..22}; do
```

---

## Running the Pipeline

### Basic Execution

```bash
bash scripts/plof_end_to_end.sh
```

### With Logging to File

```bash
bash scripts/plof_end_to_end.sh 2>&1 | tee pipeline_run.log
```

### Docker Execution (if containerized)

```bash
docker run \
  -v /data:/data \
  -v /opt/vep:/opt/vep \
  -v /opt/vep_plugins:/opt/vep_plugins \
  -e BUCKET="gs://your-bucket" \
  plof-pipeline bash /scripts/plof_end_to_end.sh
```

### Processing Specific Chromosomes

To process only chromosome 10:

```bash
export CHR="chr10"
bash scripts/pass_filter.sh
bash scripts/filterQC.sh
bash scripts/normalize.sh
bash scripts/vep_loftee.sh
```

### Resume After Interruption

Since each step checks for existing output files, simply re-run the main script:

```bash
bash scripts/plof_end_to_end.sh  # Will skip completed steps
```

To force reprocessing of a specific step, delete the output file:

```bash
rm /data/temp/bravo.pub.clean.chr21.vcf.gz
bash scripts/plof_end_to_end.sh  # Will reprocess QC filter
```

---

## Detailed Module Documentation

### STEP 0: Prep Cloud Setup

**File**: `prep_cloud.sh`

**Purpose**: Install Google Cloud SDK and dependencies for GCS file access.

**Operations**:
- Updates APT package lists
- Installs `curl`, `lsb-release`, `gnupg`
- Adds Google Cloud SDK repository
- Installs `google-cloud-sdk` (includes `gsutil`)

**When to Use**: Only needs to run once before first pipeline execution.

---

### STEP 1: PASS Filter

**File**: `pass_filter.sh`

**Purpose**: Extract only variants with FILTER=PASS status.

**Input**: 
- `${INPUTDIR}/${INPUTVCF}` (raw VCF from GCS)

**Output**:
- `${PASSVCF}` = `/data/temp/bravo.pub.pass.chr${N}.vcf.gz`

**Operations**:
1. Use `bcftools view -f PASS` to filter variants
2. Compress with gzip (`-Oz`)
3. Create tabix index (`.tbi`)

**Example**:
```bash
bcftools view chr21.bravo.pub.vcf.gz -f PASS -Oz -o pass.vcf.gz
bcftools index -t pass.vcf.gz
```

---

### STEP 2A: QC Filter

**File**: `filterQC.sh`

**Purpose**: Apply quality control filters for cohort data (skipped for site-only data).

**Input**: 
- `${PASSVCF}` (from Step 1)

**Output**:
- `${CLEANVCF}` = `/data/temp/bravo.pub.clean.chr${N}.vcf.gz`

**QC Metrics** (for cohort mode):

| Metric | Tool | Criterion | Description |
|--------|------|-----------|-------------|
| Hardy-Weinberg Equilibrium | vcftools | p ≤ 1e-8 | Variants deviating from HWE |
| Site Missingness | vcftools | ≥ 5% | Sites with >5% missing genotypes |
| Individual Missingness | vcftools | ≥ 5% | Samples with >5% missing genotypes |

**Operations**:
1. Detect sample count from VCF
2. If site-only or no samples: skip QC, pass through PASSVCF → CLEANVCF
3. If cohort mode:
   - Run HWE test
   - Calculate site-level missingness
   - Calculate individual-level missingness
   - Exclude low-quality sites and individuals
4. Output filtered VCF

---

### STEP 2B: Normalization

**File**: `normalize.sh`

**Purpose**: 
1. Normalize variant representation (split multiallelic, left-align indels)
2. Remove homozygous-reference genotypes (if present)
3. Filter to protein-coding gene regions

**Input**: 
- `${PASSVCF}` (site-only mode) or `${CLEANVCF}` (cohort mode)

**Output**:
- `${NORMVCF}` = `/data/results/vep_loftee/bravo.pub.norm.chr${N}.vcf.gz` (normalized)
- `${PROTCOD}` = `/data/temp/bravo.pub.protcod.chr${N}.vcf.gz` (protein-coding only)

**Operations**:
1. **Normalization**: `bcftools norm -m -any` with reference genome
2. **Genotype Reduction**: Remove `GT="0/0"` (hom-ref) if GT field exists
3. **Protein Coding Filter**: `bcftools view -T coding_genes.bed`

**Key Parameters**:
- `-m -any`: Split any multiallelic variants
- `NUM_SAMPLES_NORM`: Check if genotypes present before reduction

---

### STEP 3: VEP + LoFTEE Annotation

**File**: `vep_loftee.sh`

**Purpose**: Annotate variants for functional impact using Ensembl VEP and LoFTEE.

**Input**: 
- `${PROTCOD}` (normalized, protein-coding variants)

**Output**:
- `${VEPVCF}` = `/data/temp/bravo.pub.vep.chr${N}.vcf`
- `${FINALVCF}` = `/data/results/vep_loftee/bravo.pub.veploftee.chr${N}.vcf` (currently commented out)

**VEP Configuration**:

```bash
/opt/vep/src/ensembl-vep/vep \
  --assembly GRCh38              # Human reference assembly
  --cache                        # Use offline cache
  --offline                      # No internet access
  --format vcf                   # VCF input format
  --vcf                          # VCF output format
  --plugin LoF                   # Enable LoFTEE plugin
  --pick                         # Select canonical transcript
  --minimal                      # Minimal output
  --fork 7                       # Parallel processing
```

**LoFTEE Output Classes**:

| Class | Meaning | Description |
|-------|---------|-------------|
| HC | High Confidence | Likely true Loss-of-Function variant |
| LC | Low Confidence | Possible but uncertain LoF variant |
| OS | Other Splice | Splice site variant, not clearly LoF |
| NMD | Nonsense Mediated Decay | LoF but subject to NMD escape |

**Note**: Full annotation on HC variants is currently commented out in the script.

---

### OPTIONAL: Parse to TSV

**File**: `parse.sh`

**Purpose**: Convert annotated VEF from VCF format to structured tab-separated values (TSV).

**Input**: 
- VEP annotated VCFs from `/data/results/vep_loftee/`

**Output**:
- TSV files with format: `bravo.pub.vep.parsed.chr${N}.tsv`

**Operations**:
1. Filter VCF for HC and LC LoF variants
2. Extract CSQ field and expand into columns
3. Create standardized TSV with fields from VEP annotation

**TSV Columns** (example):
```
CHROM	POS	ID	REF	ALT	QUAL	FILTER	Het	Hom	Consequence	Gene	HGNC	...
```

---

## Input/Output Formats

### Input VCF

**Source**: Google Cloud Storage  
**Format**: Compressed VCF (`.vcf.gz`)  
**Example filename**: `chr21.bravo.pub.vcf.gz`

**Required VCF components**:
- CHROM, POS, REF, ALT fields
- FILTER field (with at least PASS entries)
- Optional: genotype data (GT field) for cohort mode

### Intermediate Files

| Stage | Filename Pattern | Description |
|-------|-----------------|-------------|
| Pass | `bravo.pub.pass.chr${N}.vcf.gz` | PASS-filtered variants |
| QC | `bravo.pub.clean.chr${N}.vcf.gz` | QC-filtered (if cohort) |
| Normalized | `bravo.pub.norm.chr${N}.vcf.gz` | Normalized variants |
| Protein Coding | `bravo.pub.protcod.chr${N}.vcf.gz` | Protein-coding only |

### Output VCF

**Format**: VCF with VEP CSQ annotation  
**Filename**: `bravo.pub.veploftee.chr${N}.vcf`  
**Key Fields**:
- `CSQ`: VEP consequence information
- For LoFTEE: `CSQ|LoF=HC|LC|OS|NMD` indicators

### Log Files

**Filename**: `bravo.pub.chr${N}.pipeline.log`  
**Location**: `/data/results/vep_loftee/`  
**Contents**: Timestamped messages from each step

---

## Environment Variables

### User-Defined (in main script)

| Variable | Default | Type | Description |
|----------|---------|------|-------------|
| `INPUTTYPE` | `site_only` | string | `site_only` \| `cohort` |
| `VCFBASENAME` | `bravo.pub` | string | Base name for outputs |
| `EXTENSIONNAME` | `vcf.gz` | string | VCF file extension |
| `BUCKET` | GCS path | string | Google Cloud Storage bucket |
| `THREADS` | 7 | int | Number of parallel threads |

### System-Defined (derived)

| Variable | Scope | Description |
|----------|-------|-------------|
| `CHR` | chromosome loop | Current chromosome (e.g., `chr21`) |
| `LOG` | chromosome loop | Log file path for current chromosome |
| `PASSVCF` | step 1 | PASS-filtered VCF path |
| `CLEANVCF` | step 2A | QC-filtered VCF path |
| `NORMVCF` | step 2B | Normalized VCF path |
| `PROTCOD` | step 2B | Protein-coding filtered VCF path |
| `VEPVCF` | step 3 | VEP-annotated VCF path |
| `FINALVCF` | step 3 | Final output VCF path |

### Exporting to Child Scripts

All variables are `export`ed to make them available to sourced scripts:

```bash
export CHR="chr21"
bash /scripts/pass_filter.sh  # Has access to $CHR, $LOG, etc.
```

---

## Troubleshooting

### Common Issues

#### 1. "bcftools: command not found"

**Solution**:
```bash
apt-get install -y bcftools
# or: conda install -c bioconda bcftools
```

#### 2. "Reference file not found"

**Solution**: Ensure reference files exist:
```bash
ls -lh /data/references/human_genome/
# Should show: Homo_sapiens_assembly38.fasta, .fasta.fai, coding_genes.bed
```

#### 3. "gsutil: command not found" (cloud access)

**Solution**: Re-run prep_cloud.sh or manually install:
```bash
bash scripts/prep_cloud.sh
# or: apt-get install -y google-cloud-sdk
```

#### 4. VEP plugin errors

**Solution**:
- Verify VEP installation: `/opt/vep/src/ensembl-vep/vep --version`
- Check LoFTEE plugin files:
  ```bash
  ls -lh /opt/vep_plugins/loftee/
  ```
- Re-download plugin if missing

#### 5. "Out of memory" during VEP

**Solution**:
- Reduce `THREADS` variable
- Process fewer chromosomes per run
- Filter VCF to smaller regions before VEP

#### 6. Pipeline skips steps unexpectedly

**Reason**: Output files exist from a previous run  
**Solution**: Delete intermediate files to reprocess:
```bash
rm -f /data/temp/bravo.pub.*.vcf.gz
rm -f /data/results/vep_loftee/bravo.pub.*.vcf.gz
bash scripts/plof_end_to_end.sh
```

#### 7. "Chromosome loop variables undefined in child scripts"

**Reason**: Variables not exported properly  
**Solution**: Ensure `export` keyword in parent script (should be present)

### Debug Mode

Add verbose logging to identify issues:

```bash
set -x  # Enable command tracing
bash scripts/plof_end_to_end.sh 2>&1 | tee debug.log
```

### Log File Analysis

```bash
# View latest errors
tail -50 /data/results/vep_loftee/bravo.pub.chr21.pipeline.log

# Search for specific step
grep "PASS FILTER" /data/results/vep_loftee/bravo.pub.chr21.pipeline.log

# Check for exit codes
grep -i "error\|failed" /data/results/vep_loftee/bravo.pub.chr21.pipeline.log
```

---

## Performance Considerations

### Runtime Estimates

| Step | 1M variants | 10M variants | Remarks |
|------|-------------|--------------|---------|
| PASS Filter | ~2 min | ~20 min | Linear with variant count |
| QC Filter | ~5 min | ~50 min | Slower for cohort mode |
| Normalization | ~3 min | ~30 min | Depends on complexity |
| VEP+LoFTEE | ~10 min | ~100+ min | Most time-consuming step |
| **Total** | **~20 min** | **~3+ hours** | Per chromosome |

### Optimization Tips

1. **Increase THREADS**:
   ```bash
   export THREADS=16  # If 16+ cores available
   ```

2. **Process Chromosomes in Parallel** (external script):
   ```bash
   for chr in {1..22}; do
       bash plof_end_to_end.sh &
   done
   wait
   ```

3. **Filter by Region** (pre-step):
   ```bash
   bcftools view -T regions.bed input.vcf.gz > filtered.vcf.gz
   ```

4. **Reduce VEP Annotation** (if full annotation not needed):
   - Use `--pick` (already in use) for canonical transcript only
   - Disable unused plugins
   - Use `--minimal` (already in use)

### Resource Requirements

- **Disk Space**: ~300 GB per chromosome (VCF + indexed files)
- **RAM**: 32 GB minimum (especially for VEP)
- **Network**: High bandwidth if accessing GCS frequently
- **CPU**: 7+ cores recommended for parallelization

---

## Citation

If you use this pipeline, please cite:

```
McLaren et al. (2016) "The Ensembl Variant Effect Predictor"
Genome Biology, 17:122. https://doi.org/10.1186/s13059-016-0974-4

Karczewski et al. (2020) "The ExAC Browser and ExAC constrainedProjection of human genetic variation"
Nature Reviews Genetics, 21(4):244-256. https://doi.org/10.1038/s41576-020-0210-7
```

---

## Contact & Support

For issues, questions, or contributions:
- **Code**: See script headers for author information
- **Documentation**: Refer to individual script comments
- **Bug Reports**: File an issue on the repository

---

## License

[Specify your license here - e.g., MIT, GPL, etc.]

---

## Change Log

| Date | Version | Changes |
|------|---------|---------|
| Feb 16, 2026 | 1.0 | Initial pipeline release |
| Feb 25, 2026 | 1.0.1 | Comprehensive documentation |

