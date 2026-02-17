#!/usr/bin/env bash

################################################################################
# COSIGT Pipeline - Preprocessing Reference (Docker/native version)
################################################################################

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <config_file> <threads>"
    exit 1
fi

CONFIG="$1"
THREADS="$2"

# Simple YAML Parser
get_yaml_value() {
    local key="$1"
    grep "^${key}:" "$CONFIG" | sed "s/^${key}:[[:space:]]*//" | tr -d "'" | tr -d '"'
}

INPUT_DIR=$(get_yaml_value "input")
OUTPUT_DIR=$(get_yaml_value "output")
REFERENCE=$(get_yaml_value "reference")

echo "==================================================================="
echo "COSIGT Pipeline - Preprocessing Reference"
echo "==================================================================="
echo "Reference: $REFERENCE"
echo "Threads: $THREADS"
echo "Output: $OUTPUT_DIR"
echo "==================================================================="

MERYL_REF_DB="${OUTPUT_DIR}/meryl/reference"

if [ -d "${INPUT_DIR}/meryl/reference" ]; then
    echo "Reference k-mer database already exists"
else
    mkdir -p "${OUTPUT_DIR}/meryl"
    MEM_GB=$((THREADS * 2))

    echo "Building reference k-mer database..."
    meryl count \
        k=31 \
        threads="$THREADS" \
        memory="$MEM_GB" \
        "$REFERENCE" \
        output "$MERYL_REF_DB"

    echo "✓ Reference k-mer database built successfully!"
fi

echo "==================================================================="
echo "✓ Preprocessing complete!"
echo "Database location: $MERYL_REF_DB"
echo "==================================================================="
