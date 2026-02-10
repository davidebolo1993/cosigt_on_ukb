#!/usr/bin/env bash
set -euo pipefail

################################################################################
# COSIGT Pipeline - Preprocessing Sample
# 
# Usage: ./preprocess_sample.sh <sample> <config_file> <threads> [cache_dir] [tmp_dir]
#
# Extracts unmapped reads for a single sample
################################################################################

if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 <sample> <config_file> <threads> [cache_dir] [tmp_dir]"
    echo ""
    echo "Arguments:"
    echo "  sample        Sample ID to process"
    echo "  config_file   Path to config.yaml"
    echo "  threads       Number of threads to use"
    echo "  cache_dir     Singularity cache directory [default: \$PWD/singularity_cache]"
    echo "  tmp_dir       Singularity tmp directory [default: /tmp]"
    echo ""
    echo "Example: $0 HG00438 config/config.yaml 8 ./cache /scratch/tmp"
    exit 1
fi

SAMPLE="$1"
CONFIG="$2"
THREADS="$3"
CACHE_DIR="${4:-$PWD/singularity_cache}"
TMP_DIR="${5:-/tmp}"

################################################################################
# Setup Singularity Environment
################################################################################

# Create and set cache directory
CACHE_DIR=$(readlink -f "$CACHE_DIR" 2>/dev/null || echo "$CACHE_DIR")
mkdir -p "$CACHE_DIR"
export SINGULARITY_CACHEDIR="$CACHE_DIR"

# Set tmp directory
TMP_DIR=$(readlink -f "$TMP_DIR" 2>/dev/null || echo "$TMP_DIR")
if [ ! -d "$TMP_DIR" ]; then
    echo "ERROR: tmp directory does not exist: $TMP_DIR"
    exit 1
fi
if [ ! -w "$TMP_DIR" ]; then
    echo "ERROR: tmp directory is not writable: $TMP_DIR"
    exit 1
fi
export SINGULARITY_TMPDIR="$TMP_DIR"

################################################################################
# Simple YAML Parser
################################################################################

get_yaml_value() {
    local key="$1"
    grep "^${key}:" "$CONFIG" | sed "s/^${key}:[[:space:]]*//" | tr -d "'" | tr -d '"'
}

# Load configuration
OUTPUT_DIR=$(get_yaml_value "output")
REFERENCE=$(get_yaml_value "reference")

echo "==================================================================="
echo "COSIGT Pipeline - Preprocessing Sample: $SAMPLE"
echo "==================================================================="
echo "Configuration file: $CONFIG"
echo "Threads: $THREADS"
echo "Singularity cache: $CACHE_DIR"
echo "Singularity tmp: $TMP_DIR"
echo "==================================================================="
echo ""

################################################################################
# Setup Bind Paths for Singularity
################################################################################

# Collect all directories that need to be bound
BIND_PATHS="$PWD"

# Add reference directory
REFERENCE_DIR=$(dirname "$(readlink -f "$REFERENCE")")
BIND_PATHS="${BIND_PATHS},${REFERENCE_DIR}"

# Add output directory if it's not under PWD
OUTPUT_DIR_REAL=$(readlink -f "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
if [[ "$OUTPUT_DIR_REAL" != "$PWD"* ]]; then
    BIND_PATHS="${BIND_PATHS},${OUTPUT_DIR_REAL}"
fi

# Add tmp directory
BIND_PATHS="${BIND_PATHS},${TMP_DIR}"

# Add resources directory for alignments
if [ -d "resources" ]; then
    RESOURCES_REAL=$(readlink -f "resources")
    BIND_PATHS="${BIND_PATHS},${RESOURCES_REAL}"
    
    # Add symlinked alignment directories
    if [ -d "resources/alignments" ]; then
        for item in resources/alignments/*; do
            if [ -L "$item" ]; then
                REAL_PATH=$(readlink -f "$item")
                REAL_DIR=$(dirname "$REAL_PATH")
                BIND_PATHS="${BIND_PATHS},${REAL_DIR}"
            fi
        done
    fi
fi

# Remove duplicates
BIND_PATHS=$(echo "$BIND_PATHS" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')

echo "Singularity bind paths: $BIND_PATHS"
echo ""

################################################################################
# Container Management
################################################################################

CONTAINER_DIR="${OUTPUT_DIR}/containers"
mkdir -p "$CONTAINER_DIR"

get_container() {
    local docker_uri="$1"
    local name=$(echo "$docker_uri" | sed 's|docker://||' | tr '/:' '_').sif
    local path="${CONTAINER_DIR}/${name}"
    
    if [ ! -f "$path" ]; then
        echo "[Container] Downloading: $docker_uri" >&2
        singularity pull "$path" "$docker_uri"
    fi
    echo "$path"
}

run_container() {
    local container="$1"
    shift
    singularity exec --cleanenv --bind "$BIND_PATHS" "$container" "$@"
}

echo "[Step 1] Loading samtools container..."
SAMTOOLS_C=$(get_container "docker://davidebolo1993/samtools:1.22")

################################################################################
# Find Sample Alignment
################################################################################

echo "[Step 2] Finding alignment for sample: $SAMPLE"

SAMPLE_ALIGNMENT=$(find resources/alignments -name "${SAMPLE}.*am" 2>/dev/null | head -1)

if [ -z "$SAMPLE_ALIGNMENT" ]; then
    echo "ERROR: No alignment found for sample: $SAMPLE"
    echo "       Looked in: resources/alignments/"
    exit 1
fi

echo "  Found: $SAMPLE_ALIGNMENT"

################################################################################
# Extract Unmapped Reads
################################################################################

UNMAPPED_DIR="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}"
mkdir -p "$UNMAPPED_DIR"
UNMAPPED_FASTA="${UNMAPPED_DIR}/unmapped.fasta.gz"

if [ -f "$UNMAPPED_FASTA" ]; then
    echo "[Step 3] Unmapped reads already extracted: $UNMAPPED_FASTA"
    echo "         Remove it if you want to re-extract."
else
    echo "[Step 3] Extracting unmapped reads..."
    echo "         Input: $SAMPLE_ALIGNMENT"
    echo "         Output: $UNMAPPED_FASTA"
    echo "         Threads: $THREADS"
    echo ""
    
    run_container "$SAMTOOLS_C" samtools view -u -f 4 -@ "$THREADS" -T "$REFERENCE" "$SAMPLE_ALIGNMENT" | \
    run_container "$SAMTOOLS_C" samtools sort -n -@ "$THREADS" -T "${UNMAPPED_DIR}/unmapped_tmp" - | \
    run_container "$SAMTOOLS_C" samtools fasta -0 /dev/null -@ "$THREADS" - | gzip > "$UNMAPPED_FASTA"
    
    echo ""
    echo "✓ Unmapped reads extracted successfully!"
fi


echo "==================================================================="
echo "✓ Sample preprocessing complete for: $SAMPLE"
echo ""
echo "Unmapped reads: $UNMAPPED_FASTA"
echo ""
echo "==================================================================="

