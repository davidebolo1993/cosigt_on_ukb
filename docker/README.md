# COSIGT Pipeline - Dockerized Version

COSIGT pipeline for genotyping from pangenome graphs using Docker/Singularity containers. This implementation is designed to run on HPC systems starting from pre-built variation graphs.

## Overview

This pipeline genotypes samples against pangenome graphs in four main stages:

1. **Organize**: Set up directory structure and configuration
2. **Preprocess Reference**: Build reference k-mer database (once per reference)
3. **Preprocess Samples**: Extract unmapped reads (once per sample)
4. **Process Regions**: Genotype all samples for each region

All pipeline steps run through Singularity containers, ensuring reproducibility and portability.

## Prerequisites

- Singularity (≥3.5)

## Container Images

Four containers are available on Docker Hub (auto-pulled by Singularity):

- `davidebolo1993/cosigt-organize:latest`
- `davidebolo1993/cosigt-preprocess-reference:latest`
- `davidebolo1993/cosigt-preprocess-sample:latest`
- `davidebolo1993/cosigt-process-region:latest`

## Input Files

Prepare three tab-separated files mapping your data:

### 1. Graph Map (`graph.map.tsv`)
Maps region identifiers to `.og` graph files:
```
chr1_100000_200000	/path/to/chr1_100000_200000.og
chr2_300000_400000	/path/to/chr2_300000_400000.og
```

### 2. Alignment Map (`aln.map.tsv`)
Maps alignment files to sample IDs:
```
/path/to/HG00096.cram	HG00096
/path/to/HG00171.cram	HG00171
```
**Format:** `alignment_path<TAB>sample_id`

Both BAM and CRAM formats are supported. Index files (`.bai`, `.csi`, or `.crai`) must be present.

### 3. Regions BED (`roi.bed`)
Genomic regions to genotype (standard BED format):
```
chr1	100000	200000	region1
chr2	300000	400000	region2
```

### 4. Reference Genome
Indexed reference genome in FASTA format with `.fai` index.

## Pipeline Execution

### Step 1: Organize Inputs

Create directory structure and configuration:

```bash
singularity run \
    -B "$PWD:/work,/path/to/data1,/path/to/data2" \
    docker://davidebolo1993/cosigt-organize:latest \
    -g test/graph.map.tsv \
    -r test/aln.map.tsv \
    -b test/roi.bed \
    -f test/GRCh38.primary_assembly.genome.fa \
    -o output_directory
```

**Bind mounts (`-B`):** Include all directories containing your input files. For example:
- `-B "$PWD:/work"` - current directory
- Add paths to your graph files (e.g., `/scratch/user`)
- Add paths to your alignment files (e.g., `/project/data`)

**Generated:**
- `resources/` - symlinks to all inputs organized by chromosome
- `config/config.yaml` - pipeline configuration

### Step 2: Preprocess Reference

Build reference k-mer database (run once per reference):

```bash
singularity run \
    -B "$PWD:/work" \
    docker://davidebolo1993/cosigt-preprocess-reference:latest \
    config/config.yaml 8
```

Replace `8` with number of threads to use.

**Output:** `output_directory/meryl/reference/` - 31-mer database

### Step 3: Preprocess Samples

Extract unmapped reads for each sample (parallelizable):

```bash
# Get list of samples
SAMPLES=$(grep -A 100 "^samples:" config/config.yaml | grep "^- " | sed 's/^- //')

# Process each sample
for SAMPLE in $SAMPLES; do
    singularity run \
        -B "$PWD:/work,/path/to/alignments" \
        docker://davidebolo1993/cosigt-preprocess-sample:latest \
        $SAMPLE config/config.yaml 8
done
```

**Tip:** Submit these as separate jobs on a cluster for parallel execution.

**Output:** `output_directory/samtools/fasta/{sample}/unmapped.fasta.gz`

### Step 4: Process Regions

Run genotyping for each region (parallelizable):

```bash
# Get list of regions
REGIONS=$(grep -A 100 "^regions:" config/config.yaml | grep "^- " | sed 's/^- //')

# Process each region
for REGION in $REGIONS; do
    singularity run \
        -B "$PWD:/work,/path/to/graphs,/path/to/alignments,/localscratch" \
        docker://davidebolo1993/cosigt-process-region:latest \
        $REGION config/config.yaml 12
done
```

**Important bind mounts:**
- Include all paths from organize step
- Add `/localscratch` or temp directory for intermediate files

**Tip:** Each region is independent - submit as separate cluster jobs.

**Output:** `output_directory/cosigt/{sample}/{chrom}/{region}/{region}.cosigt_genotype.tsv`

## Quick Test Example

Test with small dataset (2 samples, 1 region):

```bash
# Ensure test data paths are correct in TSV files
cat test/graph.map.small.tsv
cat test/aln.map.small.tsv

# Step 1: Organize
singularity run \
    -B "$PWD:/work,/scratch/davide.bolognini,/project/ham" \
    docker://davidebolo1993/cosigt-organize:latest \
    -g test/graph.map.small.tsv \
    -r test/aln.map.small.tsv \
    -b test/roi.small.bed \
    -f test/GRCh38.primary_assembly.genome.fa \
    -o cosigt_test

# Step 2: Preprocess reference
singularity run \
    -B "$PWD:/work" \
    docker://davidebolo1993/cosigt-preprocess-reference:latest \
    config/config.yaml 8

# Step 3: Preprocess samples
singularity run \
    -B "$PWD:/work,/project/ham" \
    docker://davidebolo1993/cosigt-preprocess-sample:latest \
    HG00096 config/config.yaml 8

singularity run \
    -B "$PWD:/work,/project/ham" \
    docker://davidebolo1993/cosigt-preprocess-sample:latest \
    HG00171 config/config.yaml 8

# Step 4: Process region
singularity run \
    -B "$PWD:/work,/scratch/davide.bolognini,/project/ham,/localscratch" \
    docker://davidebolo1993/cosigt-process-region:latest \
    chr1_103304997_103901127 config/config.yaml 12
```

## Output Structure

```
output_directory/
├── meryl/
│   └── reference/              # Reference k-mer database
├── samtools/fasta/
│   ├── HG00096/
│   │   └── unmapped.fasta.gz   # Unmapped reads
│   └── HG00171/
│       └── unmapped.fasta.gz
├── cosigt/
│   ├── HG00096/
│   │   └── chr1/
│   │       └── chr1_103304997_103901127/
│   │           └── chr1_103304997_103901127.cosigt_genotype.tsv
│   └── HG00171/
│       └── chr1/
│           └── chr1_103304997_103901127/
│               └── chr1_103304997_103901127.cosigt_genotype.tsv
└── [intermediate files in other directories]
```

The main output files are:
```
cosigt/{sample}/{chrom}/{region}/{region}.cosigt_genotype.tsv
```

## Parallelization Strategy

### For Production Runs

**Samples (Step 3):**
Submit one job per sample. All samples can run in parallel.

**Regions (Step 4):**
Submit one job per region. All regions can run in parallel. Each region job processes all samples.

### Example SLURM Script

```bash
#!/bin/bash
#SBATCH --job-name=cosigt_region
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=2:00:00

REGION=$1
CONFIG=config/config.yaml
THREADS=12

singularity run \
    -B "$PWD:/work,/scratch/user,/project/data,/localscratch" \
    docker://davidebolo1993/cosigt-process-region:latest \
    $REGION $CONFIG $THREADS
```

Submit for each region:
```bash
for REGION in $(grep -A 100 "^regions:" config/config.yaml | grep "^- " | sed 's/^- //'); do
    sbatch process_region.sh $REGION
done
```

## Per-Region Pipeline Steps

Each region job performs:

1. Extract assemblies from graph
2. Build BWA-MEM2 index for assemblies
3. Convert graph to GFA and extract paths/nodes
4. Create BED files for graph features
5. Build k-mer filter index (graph-specific k-mers)
6. Filter low-complexity nodes (coverage outliers)
7. Cluster haplotypes by sequence similarity
8. For each sample:
   - Extract mapped reads from region
   - Filter unmapped reads using graph k-mers
   - Realign reads to graph assemblies (BWA-MEM2)
   - Project alignments to graph (GAF format via gafpack)
   - Calculate node coverage
   - Genotype using COSIGT algorithm

## License

This project is licensed under the MIT License.

## Contact

For issues, questions, or contributions, please open an issue on GitHub or contact the developers.

