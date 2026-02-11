# COSIGT on UKB-RAP

COSIGT bash pipeline for genotyping from pangenome graphs. This implementation is designed to run on the UK Biobank Research Analysis Platform (UKB-RAP) and starts from pre-built variation graphs.

## Overview

This pipeline genotypes samples against pangenome graphs in three main stages:

1. **Setup**: Organize input files and download containers
2. **Preprocessing**: Build reference k-mer database and extract unmapped reads per sample
3. **Genotyping**: Process each genomic region independently across all samples

The pipeline is designed for parallel execution on HPC systems using SLURM.

## Prerequisites

- Singularity
- SLURM job scheduler
- Python 3.6+

## Input Files

Prepare three tab-separated files mapping your data:

1. **Graph map** (`graph.map.tsv`): region to `.og` graph file
   ```
   chr1_100000_200000    /path/to/chr1_100000_200000.og
   chr2_300000_400000    /path/to/chr2_300000_400000.og
   ```

2. **Alignment map** (`aln.map.tsv`): alignment file to sample ID
   ```
   /path/to/HG00096.cram    HG00096
   /path/to/HG00171.cram    HG00171
   ```

3. **Regions BED** (`roi.bed`): genomic regions to genotype (3-5 columns)
   ```
   chr1    100000    200000    region1
   chr2    300000    400000    region2
   ```

You'll also need an indexed reference genome (FASTA + FAI).

## Quick Start

### 1. Organize inputs

This creates the working directory structure and symlinks all inputs:

```bash
python src/organize.py \
    -g test/graph.map.tsv \
    -r test/aln.map.tsv \
    -b test/roi.bed \
    -f test/GRCh38.primary_assembly.genome.fa \
    -o output_directory
```

This generates:
- `resources/`: symlinks to graphs, alignments, regions, and reference
- `config/config.yaml`: pipeline configuration file

### 2. Download containers

Pre-download all Singularity containers to avoid conflicts when running parallel jobs:

```bash
module load singularity/3.8.5
bash src/download_containers.sh config/config.yaml
```

This only needs to be done once.

### 3. Preprocess reference

Build the reference k-mer database (run once per reference genome):

```bash
bash src/preprocess_reference.sh config/config.yaml 8
```

Uses meryl to count 31-mers in the reference. Adjust threads based on available resources.

### 4. Preprocess samples

Extract unmapped reads from each sample's alignment. This submits one SLURM job per sample:

```bash
bash src/submit_sample_preprocessing.sh \
    -c config/config.yaml \
    -t 8 \
    -m 4G \
    --time 00:30:00
```

Check options with `bash src/submit_sample_preprocessing.sh --help`

Jobs can run in parallel since each processes a different sample.

### 5. Process regions

Run COSIGT genotyping on each region. This submits one SLURM job per region:

```bash
bash src/submit_region_processing.sh \
    -c config/config.yaml \
    -t 8 \
    -m 20G \
    --time 00:30:00
```

Each region job processes all samples for that region. Regions are independent and run in parallel.

## Small Test Example

Use the included test data to verify the pipeline works:

```bash
# Setup
python src/organize.py \
    -g test/graph.map.small.tsv \
    -r test/aln.map.small.tsv \
    -b test/roi.small.bed \
    -f test/GRCh38.primary_assembly.genome.fa #this is not included in the folder \
    -o cosigt_test

module load singularity/3.8.5
bash src/download_containers.sh config/config.yaml
bash src/preprocess_reference.sh config/config.yaml 8
bash src/submit_sample_preprocessing.sh -t 8 -m 5G
bash src/submit_region_processing.sh -t 8 -m 10G
```

## Output Structure

Results are organized by sample and region:

```
output_directory/
├── containers/           # Singularity .sif files
├── meryl/
│   └── reference/       # Reference k-mer database
├── samtools/fasta/
│   └── {sample}/
│       └── unmapped.fasta.gz  # Per-sample unmapped reads
├──...
│
└── cosigt/
    └── {sample}/{chrom}/{region}/
        └── {region}.cosigt_genotype.tsv  # Final genotypes
```

The main output is `cosigt/{sample}/{chrom}/{region}/{region}.cosigt_genotype.tsv` containing genotype calls.

## Pipeline Steps (Per Region)

Each region job performs these steps for all samples:

1. Extract assemblies from graph
2. Build BWA-MEM2 index for assemblies
3. Convert graph to GFA and extract node information
4. Build k-mer filter index (graph-specific kmers)
5. Filter low-complexity/coverage-outlier nodes
6. Cluster haplotypes by similarity
7. For each sample:
   - Extract mapped reads from region
   - Combine with k-mer-filtered unmapped reads
   - Realign to graph assemblies
   - Project alignments to graph (GAF format)
   - Calculate node coverage
   - Genotype using COSIGT algorithm

## Citation

If you use this pipeline, please cite the COSIGT paper [XXX].

