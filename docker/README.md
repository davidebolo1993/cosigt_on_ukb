# COSIGT Pipeline - Dockerized Version

COSIGT pipeline for genotyping from pangenome graphs using Docker/Singularity containers. This implementation is designed to run on HPC systems starting from pre-built variation graphs.

## Overview

This pipeline genotypes samples against pangenome graphs in three main stages:

1. **Preprocess Reference**: Build reference k-mer database (once per reference)
2. **Preprocess Samples**: Extract unmapped reads (once per sample)
3. **Process Regions**: Genotype all samples for each region

All pipeline steps run through Singularity containers, ensuring reproducibility and portability.

## Prerequisites

- Singularity (≥3.5)

## Container Images

Three containers are available on Docker Hub (auto-pulled by Singularity):

- `davidebolo1993/cosigt-preprocess-reference:latest`
- `davidebolo1993/cosigt-preprocess-sample:latest`
- `davidebolo1993/cosigt-process-region:latest`

## Input Files

Prepare the following files:

### 1. Configuration File (`config.yaml`)

```yaml
output: /path/to/output_directory
reference: /path/to/reference.fa
alignment_map: /path/to/aln.map.tsv
graph_map: /path/to/graph.map.tsv

samples:
  - HG00096
  - HG00171

regions:
  - chr1_103304997_103901127
  - chr2_300000_400000
```

### 2. Alignment Map (`aln.map.tsv`)

Tab-separated file mapping samples to alignment files:

```tsv
HG00096	/path/to/HG00096.cram
HG00171	/path/to/HG00171.cram
```

**Format:** `sample<TAB>alignment_path`

Both BAM and CRAM formats are supported. Index files (`.bai`, `.csi`, or `.crai`) must be present in the same directory.

### 3. Graph Map (`graph.map.tsv`)

Tab-separated file mapping regions to graph files:

```tsv
chr1_103304997_103901127	/path/to/chr1_103304997_103901127.og
chr2_300000_400000	/path/to/chr2_300000_400000.og
```

**Format:** `region<TAB>graph_path`

Graph files must be in ODGI format (`.og`).

### 4. Reference Genome

Indexed reference genome in FASTA format with `.fai` index.

## Pipeline Execution

### Step 1: Preprocess Reference

Build reference k-mer database (run once per reference):

```bash
singularity run \
    -B "$PWD:/work" \
    docker://davidebolo1993/cosigt-preprocess-reference:latest \
    config.yaml 8
```

Replace `8` with number of threads to use.

**Output:** `output_directory/meryl/reference/` - 31-mer database

**Note:** This step only needs to be run once. The meryl database can be reused for all samples and regions.

### Step 2: Preprocess Samples

Extract unmapped reads for each sample (parallelizable):

```bash
# Get list of samples from config
SAMPLES=$(grep -A 10000 "^samples:" config.yaml | grep "^- " | sed 's/^- //')

# Process each sample
for SAMPLE in $SAMPLES; do
    singularity run \
        -B "$PWD:/work" \
        docker://davidebolo1993/cosigt-preprocess-sample:latest \
        $SAMPLE config.yaml 8
done
```

**Important:** Make sure all paths in `aln.map.tsv` are accessible. Add bind mounts with `-B` if needed:

```bash
singularity run \
    -B "$PWD:/work,/path/to/alignments" \
    docker://davidebolo1993/cosigt-preprocess-sample:latest \
    $SAMPLE config.yaml 8
```

**Tip:** Submit these as separate jobs on a cluster for parallel execution.

**Output:** `output_directory/samtools/fasta/{sample}/unmapped.fasta.gz`

### Step 3: Process Regions

Run genotyping for each region:

```bash
# Get list of regions from config
REGIONS=$(grep -A 100 "^regions:" config.yaml | grep "^- " | sed 's/^- //')

# Process each region
for REGION in $REGIONS; do
    singularity run \
        -B "$PWD:/work,/path/to/graphs,/path/to/alignments,/localscratch" \
        docker://davidebolo1993/cosigt-process-region:latest \
        $REGION config.yaml 12
done
```

**Important bind mounts:**
- Include all directories containing your input files
- Add `/localscratch` or temp directory for intermediate files
- Ensure graph files (`.og`) and alignments are accessible

**Tip:** Each region is independent - submit as separate cluster jobs.

**Output:** `output_directory/cosigt/{sample}/{chrom}/{region}/{region}.cosigt_genotype.tsv`

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