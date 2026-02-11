#!/usr/bin/env bash

################################################################################
# COSIGT Pipeline - Pre-download All Containers
#
# Usage: ./download_containers.sh <config_file> [cache_dir]
#
# Downloads all required Singularity containers before running the pipeline
################################################################################

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <config_file> [cache_dir]"
    echo ""
    echo "Arguments:"
    echo "  config_file   Path to config.yaml"
    echo "  cache_dir     Singularity cache directory [default: \$PWD/singularity_cache]"
    echo ""
    echo "Example: $0 config/config.yaml ./singularity_cache"
    exit 1
fi

CONFIG="$1"
CACHE_DIR="${2:-$PWD/singularity_cache}"

################################################################################
# Setup Singularity Environment
################################################################################

CACHE_DIR=$(readlink -f "$CACHE_DIR" 2>/dev/null || echo "$CACHE_DIR")
mkdir -p "$CACHE_DIR"
export SINGULARITY_CACHEDIR="$CACHE_DIR"

################################################################################
# Simple YAML Parser
################################################################################

get_yaml_value() {
    local key="$1"
    grep "^${key}:" "$CONFIG" | sed "s/^${key}:[[:space:]]*//" | tr -d "'" | tr -d '"'
}

# Load configuration
OUTPUT_DIR=$(get_yaml_value "output")

echo "==================================================================="
echo "COSIGT Pipeline - Pre-downloading Containers"
echo "==================================================================="
echo "Config:         $CONFIG"
echo "Cache dir:      $CACHE_DIR"
echo "Container dir:  ${OUTPUT_DIR}/containers"
echo "==================================================================="
echo ""

################################################################################
# Container Management
################################################################################

CONTAINER_DIR="${OUTPUT_DIR}/containers"
mkdir -p "$CONTAINER_DIR"

# List of all required containers
declare -A CONTAINERS=(
    ["samtools"]="docker://davidebolo1993/samtools:1.22"
    ["bwa-mem2"]="docker://davidebolo1993/bwa-mem2:2.2.1"
    ["bedtools"]="docker://davidebolo1993/bedtools:2.31.0"
    ["kfilt"]="docker://davidebolo1993/kfilt:0.1.1"
    ["odgi"]="docker://pangenome/odgi:1753347183"
    ["gfainject"]="docker://davidebolo1993/gfainject:0.2.1"
    ["gafpack"]="docker://davidebolo1993/gafpack:0.1.3"
    ["panplexity"]="docker://davidebolo1993/panplexity:0.1.1"
    ["cosigt"]="docker://davidebolo1993/cosigt:0.1.7"
    ["renv"]="docker://davidebolo1993/renv:4.3.3"
)

TOTAL=${#CONTAINERS[@]}
CURRENT=0
DOWNLOADED=0
SKIPPED=0

for name in "${!CONTAINERS[@]}"; do
    ((CURRENT++))
    docker_uri="${CONTAINERS[$name]}"
    sif_name=$(echo "$docker_uri" | sed 's|docker://||' | tr '/:' '_').sif
    sif_path="${CONTAINER_DIR}/${sif_name}"
    
    echo "[${CURRENT}/${TOTAL}] ${name}: ${docker_uri}"
    
    if [ -f "$sif_path" ]; then
        echo "  ✓ Already downloaded: $sif_path"
        ((SKIPPED++))
    else
        echo "  → Downloading to: $sif_path"
        if singularity pull "$sif_path" "$docker_uri"; then
            echo "  ✓ Downloaded successfully"
            ((DOWNLOADED++))
        else
            echo "  ✗ ERROR: Failed to download $docker_uri"
            exit 1
        fi
    fi
    echo ""
done

################################################################################
# Summary
################################################################################

echo "==================================================================="
echo "✓ Container download complete!"
echo "==================================================================="
echo "Total containers:  $TOTAL"
echo "Downloaded:        $DOWNLOADED"
echo "Already present:   $SKIPPED"
echo ""
echo "Container location: $CONTAINER_DIR"
echo ""
echo "==================================================================="

