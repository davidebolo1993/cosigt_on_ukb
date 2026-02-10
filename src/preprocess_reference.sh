#!/usr/bin/env bash
set -euo pipefail

################################################################################
# COSIGT Pipeline - Preprocessing Reference (Run Once)
# 
# Usage: src/preprocess_reference.sh <config_file> <threads> [cache_dir] [tmp_dir]
#
# Builds shared reference k-mer database before region/sample processing
################################################################################

if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <config_file> <threads> [cache_dir] [tmp_dir]"
    echo ""
    echo "Arguments:"
    echo "  config_file   Path to config.yaml"
    echo "  threads       Number of threads to use"
    echo "  cache_dir     Singularity cache directory [default: \$PWD/singularity_cache]"
    echo "  tmp_dir       Singularity tmp directory [default: /tmp]"
    echo ""
    echo "Example: $0 config/config.yaml 8 ./cache /scratch/tmp"
    exit 1
fi

CONFIG="$1"
THREADS="$2"
CACHE_DIR="${3:-$PWD/singularity_cache}"
TMP_DIR="${4:-/tmp}"

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
echo "COSIGT Pipeline - Preprocessing Reference"
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

echo "[Step 1] Downloading kfilt container..."
KFILT_C=$(get_container "docker://davidebolo1993/kfilt:0.1.1")

################################################################################
# Build Reference K-mer Database
################################################################################

MERYL_REF_DB="${OUTPUT_DIR}/meryl/reference"

if [ -d "$MERYL_REF_DB" ]; then
    echo "[Step 2] Reference k-mer database already exists: $MERYL_REF_DB"
    echo "         Remove it if you want to rebuild."
else
    echo "[Step 2] Building reference k-mer database..."
    echo "         Reference: $REFERENCE"
    echo "         Threads: $THREADS"
    echo "         Output: $MERYL_REF_DB"
    
    mkdir -p "${OUTPUT_DIR}/meryl"
    
    # Calculate memory (2GB per thread is reasonable for meryl)
    MEM_GB=$((${THREADS} * 2))
    
    echo "         Memory: ${MEM_GB}GB"
    echo ""
    
    run_container "$KFILT_C" meryl count \
        k=31 \
        threads="$THREADS" \
        memory="$MEM_GB" \
        "$REFERENCE" \
        output "$MERYL_REF_DB"
    
    echo ""
    echo "✓ Reference k-mer database built successfully!"
fi


echo "==================================================================="
echo "✓ Preprocessing complete!"
echo ""
echo "Database location: $MERYL_REF_DB"
echo "Container cache: $CACHE_DIR"
echo ""
echo "==================================================================="

