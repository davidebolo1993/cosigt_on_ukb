#!/usr/bin/env bash

set -euo pipefail

################################################################################
# COSIGT Pipeline - Region Processing (Docker/native version)
# 
# No Singularity - all tools available natively in container
################################################################################

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <region> <config_file> <threads>"
    echo ""
    echo "Example: $0 chr1_100000_200000 config/config.yaml 12"
    exit 1
fi

REGION="$1"
CONFIG="$2"
THREADS="$3"

################################################################################
# Simple YAML Parser
################################################################################

get_yaml_value() {
    local key="$1"
    grep "^${key}:" "$CONFIG" | sed "s/^${key}:[[:space:]]*//" | tr -d "'" | tr -d '"'
}

get_yaml_list() {
    local key="$1"
    awk "/^${key}:/{flag=1;next}/^[a-zA-Z]/{flag=0}flag" "$CONFIG" | \
        sed 's/^[[:space:]]*-[[:space:]]*//' | tr -d "'" | tr -d '"' | grep -v '^[[:space:]]*$'
}

# Load configuration
OUTPUT_DIR=$(get_yaml_value "output")
REFERENCE=$(get_yaml_value "reference")
PANSN_PREFIX=$(get_yaml_value "pansn_prefix")

# Load samples
SAMPLES=()
while IFS= read -r line; do
    [ -n "$line" ] && SAMPLES+=("$line")
done < <(get_yaml_list "samples")

# Extract chromosome from region
CHROM=$(echo "$REGION" | rev | cut -d'_' -f3- | rev)

echo "==================================================================="
echo "COSIGT Pipeline - Region: $REGION (Chromosome: $CHROM)"
echo "==================================================================="
echo "Samples: ${#SAMPLES[@]} | Threads: $THREADS"
echo "==================================================================="
echo ""

################################################################################
# Validate Prerequisites
################################################################################

MERYL_REF_DB="${OUTPUT_DIR}/meryl/reference"
if [ ! -d "$MERYL_REF_DB" ] || [ ! -f "${MERYL_REF_DB}/merylIndex" ]; then
    echo "ERROR: Reference k-mer database not found!"
    echo "  Expected: $MERYL_REF_DB"
    echo ""
    echo "Please run preprocessing first"
    exit 1
fi

################################################################################
# Step 1: Extract Assemblies from Graph
################################################################################

echo "[Step 1] Processing graph and extracting assemblies"

ALLELES_DIR="${OUTPUT_DIR}/alleles/${CHROM}/${REGION}"
mkdir -p "$ALLELES_DIR"

INPUT_GRAPH="resources/assemblies/${CHROM}/${REGION}.og"
if [ ! -f "$INPUT_GRAPH" ]; then
    echo "ERROR: Graph not found: $INPUT_GRAPH"
    exit 1
fi

ALLELES_FASTA="${ALLELES_DIR}/${REGION}.fasta.gz"
if [ ! -f "$ALLELES_FASTA" ]; then
    echo "  Extracting assemblies from graph..."
    odgi paths -i "$INPUT_GRAPH" -f | bgzip -c > "$ALLELES_FASTA"
    samtools faidx "$ALLELES_FASTA"
fi

# Build BWA-MEM2 index
echo "  Building BWA-MEM2 index..."
if [ ! -f "${ALLELES_FASTA}.bwt.2bit.64" ]; then
    bwa-mem2 index "$ALLELES_FASTA"
fi

################################################################################
# Step 2: ODGI Graph Utilities
################################################################################

echo "[Step 2] ODGI utilities and graph analysis"

ODGI_DIR="${OUTPUT_DIR}/odgi"
mkdir -p "${ODGI_DIR}/view/${CHROM}/${REGION}"
mkdir -p "${ODGI_DIR}/paths/${CHROM}/${REGION}"

GFA_FILE="${ODGI_DIR}/view/${CHROM}/${REGION}/${REGION}.gfa.gz"
PATHS_FILE="${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}.tsv.gz"
LENGTH_FILE="${ODGI_DIR}/view/${CHROM}/${REGION}/${REGION}.node.length.tsv"

if [ ! -f "$GFA_FILE" ]; then
    echo "  Converting to GFA and extracting paths..."
    odgi view -i "$INPUT_GRAPH" -g --threads "$THREADS" | gzip > "$GFA_FILE"
    zgrep '^S' "$GFA_FILE" | awk '{print "node."$2"\t"length($3)}' > "$LENGTH_FILE"
    odgi paths -i "$INPUT_GRAPH" -H | cut -f 1,4- | gzip > "$PATHS_FILE"
fi

################################################################################
# Step 3: Prepare BED Files
################################################################################

echo "[Step 3] Preparing BED files"

BEDTOOLS_DIR="${OUTPUT_DIR}/bedtools"
mkdir -p "${BEDTOOLS_DIR}/reference_bed/${CHROM}/${REGION}"
mkdir -p "${BEDTOOLS_DIR}/alignment_bed/${CHROM}/${REGION}"

REF_BED="${BEDTOOLS_DIR}/reference_bed/${CHROM}/${REGION}/${REGION}.bed.gz"
ALIGN_BED="${BEDTOOLS_DIR}/alignment_bed/${CHROM}/${REGION}/${REGION}.bed.gz"

INPUT_BED="resources/regions/${CHROM}/${REGION}.bed"

if [ ! -f "$REF_BED" ]; then
    grep -w "$CHROM" "$INPUT_BED" | gzip > "$REF_BED"
fi

if [ ! -f "$ALIGN_BED" ]; then
    TMP_SORTED=$(mktemp)
    bedtools sort -i "$INPUT_BED" > "$TMP_SORTED"
    bedtools intersect -a "$REF_BED" -b "$TMP_SORTED" -nonamecheck -u | gzip > "$ALIGN_BED"
    bedtools intersect -a "$TMP_SORTED" -b "$REF_BED" -nonamecheck -v | gzip >> "$ALIGN_BED"
    rm "$TMP_SORTED"
fi

################################################################################
# Step 4: Build K-mer Filter Index
################################################################################

echo "[Step 4] Building k-mer filter index"

KFILT_DIR="${OUTPUT_DIR}/kfilt/index/${CHROM}/${REGION}"
mkdir -p "$KFILT_DIR"
KFILT_IDX="${KFILT_DIR}/${REGION}.kfilt.idx"

if [ ! -f "$KFILT_IDX" ]; then
    MERYL_ALLELES="${OUTPUT_DIR}/meryl/${CHROM}/${REGION}/${REGION}"
    mkdir -p "$MERYL_ALLELES"
    MERYL_DIFF="${OUTPUT_DIR}/meryl/${CHROM}/${REGION}/${REGION}_unique"
    UNIQUE_KMERS="${OUTPUT_DIR}/meryl/${CHROM}/${REGION}/${REGION}.unique_kmers.txt"

    echo "  Computing allele-specific k-mers..."
    MEM_GB=$((THREADS < 8 ? THREADS * 2 : 16))

    meryl count k=31 threads="$THREADS" memory="$MEM_GB" "$ALLELES_FASTA" output "$MERYL_ALLELES"
    meryl difference "$MERYL_ALLELES" "$MERYL_REF_DB" output "$MERYL_DIFF"
    meryl print "$MERYL_DIFF" > "$UNIQUE_KMERS"
    kfilt build -k "$UNIQUE_KMERS" -K 31 -o "$KFILT_IDX"

    rm -rf "$MERYL_ALLELES" "$MERYL_DIFF"
fi

################################################################################
# Step 5: Node Filtering
################################################################################

echo "[Step 5] Filtering low-complexity nodes"

PANPLEXITY_DIR="${OUTPUT_DIR}/panplexity/${CHROM}/${REGION}"
mkdir -p "$PANPLEXITY_DIR"
PANPLEXITY_MASK="${PANPLEXITY_DIR}/${REGION}.mask.tsv"

if [ ! -f "$PANPLEXITY_MASK" ]; then
    panplexity \
        --input-gfa "$GFA_FILE" \
        -t auto -k 16 -w 100 -d 100 \
        --complexity linguistic \
        -m "$PANPLEXITY_MASK" \
        --threads "$THREADS"
fi

COMBINED_MASK="${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}.mask.tsv"

if [ ! -f "$COMBINED_MASK" ]; then
    echo "  Filtering coverage outliers..."
    Rscript /usr/local/bin/coverage_outliers.r \
        "$PATHS_FILE" \
        "${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}" \
        "$LENGTH_FILE" \
        "$PANPLEXITY_MASK"
fi

################################################################################
# Step 6: Clustering
################################################################################

echo "[Step 6] Computing path dissimilarities and clustering"

DISSIM_DIR="${ODGI_DIR}/dissimilarity/${CHROM}/${REGION}"
mkdir -p "$DISSIM_DIR"
DISSIM_FILE="${DISSIM_DIR}/${REGION}.tsv.gz"

if [ ! -f "$DISSIM_FILE" ]; then
    odgi similarity -i "$INPUT_GRAPH" --all --distances --threads "$THREADS" | gzip > "$DISSIM_FILE"
fi

CLUSTER_DIR="${OUTPUT_DIR}/cluster/${CHROM}/${REGION}"
mkdir -p "$CLUSTER_DIR"
CLUSTER_JSON="${CLUSTER_DIR}/${REGION}.clusters.json"

if [ ! -f "$CLUSTER_JSON" ]; then
    Rscript /usr/local/bin/cluster.r \
        "$DISSIM_FILE" "$CLUSTER_JSON" "automatic" "100.0" "1"
fi

################################################################################
# Step 7: Process Samples
################################################################################

echo "[Step 7] Processing ${#SAMPLES[@]} samples"
echo ""

process_sample() {
    local SAMPLE="$1"
    echo "  [Sample: $SAMPLE]"

    SAMPLE_ALIGNMENT=$(find resources/alignments -name "${SAMPLE}.*am" 2>/dev/null | head -1)
    if [ -z "$SAMPLE_ALIGNMENT" ]; then
        echo "    WARNING: No alignment found for $SAMPLE, skipping"
        return
    fi

    # 7.1: Extract mapped reads
    SAMTOOLS_DIR="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$SAMTOOLS_DIR"
    MAPPED_FASTA="${SAMTOOLS_DIR}/${REGION}.mapped.fasta.gz"

    if [ ! -f "$MAPPED_FASTA" ]; then
        echo "    Extracting mapped reads..."
        samtools view -T "$REFERENCE" -@ "$THREADS" -L "$ALIGN_BED" -M -b "$SAMPLE_ALIGNMENT" | \
        samtools sort -n -@ "$THREADS" -T "${SAMTOOLS_DIR}/${REGION}" - | \
        samtools fasta -@ "$THREADS" - | gzip > "$MAPPED_FASTA"
    fi

    # 7.2: Use pre-extracted unmapped reads
    UNMAPPED_FASTA="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/unmapped.fasta.gz"
    if [ ! -f "$UNMAPPED_FASTA" ]; then
        echo "    ERROR: Unmapped reads not found for $SAMPLE"
        echo "    Expected: $UNMAPPED_FASTA"
        return 1
    fi

    # 7.3: Filter unmapped reads
    KFILT_SAMPLE_DIR="${OUTPUT_DIR}/kfilt/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$KFILT_SAMPLE_DIR"
    FILTERED_UNMAPPED="${KFILT_SAMPLE_DIR}/${REGION}.unmapped.fasta.gz"

    if [ ! -f "$FILTERED_UNMAPPED" ]; then
        echo "    Filtering unmapped reads..."
        kfilt filter \
            -I "$UNMAPPED_FASTA" -o "$FILTERED_UNMAPPED" \
            -f fasta -z -i "$KFILT_IDX" -n 1 -m 0 -t "$THREADS"
    fi

    # 7.4: Combine reads
    COMBINE_DIR="${OUTPUT_DIR}/combine/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$COMBINE_DIR"
    COMBINED_FASTA="${COMBINE_DIR}/${REGION}.fasta.gz"

    if [ ! -f "$COMBINED_FASTA" ]; then
        cat "$MAPPED_FASTA" "$FILTERED_UNMAPPED" > "$COMBINED_FASTA"
    fi

    # 7.5: Realign and project
    BWA_DIR="${OUTPUT_DIR}/bwa-mem2/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$BWA_DIR"
    REALIGNED_SAM="${BWA_DIR}/${REGION}.realigned.sam"

    GFAINJECT_DIR="${OUTPUT_DIR}/gfainject/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$GFAINJECT_DIR"
    GAF_FILE="${GFAINJECT_DIR}/${REGION}.gaf.gz"

    if [ ! -f "$GAF_FILE" ]; then
        echo "    Realigning reads..."
        bwa-mem2 mem -t "$THREADS" -p -h 10000 -o "$REALIGNED_SAM" "$ALLELES_FASTA" "$COMBINED_FASTA"

        echo "    Projecting to graph..."
        gfainject --gfa "$GFA_FILE" --sam "$REALIGNED_SAM" --alt-hits 10000 | gzip > "$GAF_FILE"
        rm "$REALIGNED_SAM"
    fi

    # 7.6: Calculate coverage
    GAFPACK_DIR="${OUTPUT_DIR}/gafpack/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$GAFPACK_DIR"
    COVERAGE_FILE="${GAFPACK_DIR}/${REGION}.gafpack.gz"

    if [ ! -f "$COVERAGE_FILE" ]; then
        echo "    Computing node coverage..."
        gafpack \
            --gfa "$GFA_FILE" --gaf "$GAF_FILE" \
            --len-scale --weight-queries | gzip > "$COVERAGE_FILE"
    fi

    # 7.7: Genotype
    COSIGT_DIR="${OUTPUT_DIR}/cosigt/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$COSIGT_DIR"
    GENOTYPE_FILE="${COSIGT_DIR}/${REGION}.cosigt_genotype.tsv"

    if [ ! -f "$GENOTYPE_FILE" ]; then
        echo "    Genotyping..."
        cosigt \
            -p "$PATHS_FILE" \
            -g "$COVERAGE_FILE" \
            -c "$CLUSTER_JSON" \
            -o "$COSIGT_DIR" \
            -i "$SAMPLE" \
            -m "$COMBINED_MASK"
    fi

    echo "    ✓ Completed"
}

# Process all samples
for SAMPLE in "${SAMPLES[@]}"; do
    process_sample "$SAMPLE"
done

################################################################################
# Completion
################################################################################

echo ""
echo "==================================================================="
echo "✓ Region $REGION completed successfully!"
echo "Results: ${OUTPUT_DIR}/cosigt/*/${CHROM}/${REGION}/"
echo "==================================================================="

