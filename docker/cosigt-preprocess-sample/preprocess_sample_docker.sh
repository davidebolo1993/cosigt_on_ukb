#!/usr/bin/env bash

set -euo pipefail

################################################################################
# COSIGT Pipeline - Preprocessing Sample (Docker/native version)
################################################################################

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <sample> <config_file> <threads>"
    exit 1
fi

SAMPLE="$1"
CONFIG="$2"
THREADS="$3"

# Simple YAML Parser
get_yaml_value() {
    local key="$1"
    grep "^${key}:" "$CONFIG" | sed "s/^${key}:[[:space:]]*//" | tr -d "'" | tr -d '"'
}

OUTPUT_DIR=$(get_yaml_value "output")
REFERENCE=$(get_yaml_value "reference")

echo "==================================================================="
echo "COSIGT Pipeline - Preprocessing Sample: $SAMPLE"
echo "==================================================================="

# Find alignment
SAMPLE_ALIGNMENT=$(find resources/alignments -name "${SAMPLE}.*am" 2>/dev/null | head -1)

if [ -z "$SAMPLE_ALIGNMENT" ]; then
    echo "ERROR: No alignment found for sample: $SAMPLE"
    exit 1
fi

echo "Found alignment: $SAMPLE_ALIGNMENT"

# Extract unmapped reads
UNMAPPED_DIR="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}"
mkdir -p "$UNMAPPED_DIR"
UNMAPPED_FASTA="${UNMAPPED_DIR}/unmapped.fasta.gz"

if [ -f "$UNMAPPED_FASTA" ]; then
    echo "Unmapped reads already extracted"
else
    echo "Extracting unmapped reads..."
    samtools view -u -f 4 -@ "$THREADS" -T "$REFERENCE" "$SAMPLE_ALIGNMENT" | \
    samtools sort -n -@ "$THREADS" -T "${UNMAPPED_DIR}/unmapped_tmp" - | \
    samtools fasta -0 /dev/null -@ "$THREADS" - | gzip > "$UNMAPPED_FASTA"

    echo "✓ Unmapped reads extracted successfully!"
fi

echo "==================================================================="
echo "✓ Sample preprocessing complete"
echo "Unmapped reads: $UNMAPPED_FASTA"
echo "==================================================================="

