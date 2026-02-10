#!/bin/bash
set -euo pipefail

################################################################################
# Submit Region Processing Jobs
#
# Usage: ./submit_region_processing.sh [options]
#
# Submits one processing job per region listed in config.yaml
################################################################################

# Default values
CONFIG="config/config.yaml"
CACHE_DIR="$PWD/singularity_cache"
TMP_DIR="/tmp"
THREADS=12
MEMORY="32G"
TIME="12:00:00"
PARTITION="cpuq"

# Usage function
show_usage() {
    cat << EOF
Usage: $0 [options]

Options:
    -c, --config FILE       Config file [default: config/config.yaml]
    --cache-dir DIR         Singularity cache directory [default: \$PWD/singularity_cache]
    --tmp-dir DIR           Singularity tmp directory [default: /tmp]
    -t, --threads N         Threads per job [default: 12]
    -m, --memory MEM        Memory per job [default: 32G]
    --time TIME             Time limit per job [default: 12:00:00]
    -p, --partition PART    SLURM partition [default: cpuq]
    -h, --help              Show this help message

EOF
}

# Parse command-line arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --tmp-dir)
            TMP_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        -p|--partition)
            PARTITION="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

################################################################################
# Validate inputs
################################################################################

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

mkdir -p logs

echo "==================================================================="
echo "Submitting Region Processing Jobs"
echo "==================================================================="
echo "Config:     $CONFIG"
echo "Cache dir:  $CACHE_DIR"
echo "Tmp dir:    $TMP_DIR"
echo "Threads:    $THREADS"
echo "Memory:     $MEMORY"
echo "Time limit: $TIME"
echo "Partition:  $PARTITION"
echo "==================================================================="
echo ""

################################################################################
# Validate prerequisites
################################################################################

OUTPUT_DIR=$(grep "^output:" "$CONFIG" | sed "s/^output:[[:space:]]*//" | tr -d "'" | tr -d '"')

# Check reference k-mer database
MERYL_REF_DB="${OUTPUT_DIR}/meryl/reference"

if [ ! -d "$MERYL_REF_DB" ] || [ ! -f "${MERYL_REF_DB}/merylIndex" ]; then
    echo "ERROR: Reference k-mer database not found: $MERYL_REF_DB"
    echo "       Please run: bash src/preprocess_reference.sh first"
    exit 1
fi

echo "✓ Reference k-mer database found"

# Check unmapped reads for all samples
MISSING_SAMPLES=()
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    UNMAPPED_FILE="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/unmapped.fasta.gz"
    if [ ! -f "$UNMAPPED_FILE" ]; then
        MISSING_SAMPLES+=("$SAMPLE")
    fi
done < <(awk '/^samples:/{flag=1;next}/^[a-zA-Z]/{flag=0}flag' "$CONFIG" | \
    sed 's/^[[:space:]]*-[[:space:]]*//' | tr -d "'" | tr -d '"' | grep -v '^[[:space:]]*$')

if [ ${#MISSING_SAMPLES[@]} -gt 0 ]; then
    echo ""
    echo "ERROR: Unmapped reads not found for ${#MISSING_SAMPLES[@]} sample(s):"
    for SAMPLE in "${MISSING_SAMPLES[@]}"; do
        echo "  - $SAMPLE"
    done
    echo ""
    echo "Please run sample preprocessing first:"
    echo "  bash submit_sample_preprocessing.sh"
    exit 1
fi

echo "✓ Unmapped reads found for all samples"
echo ""

################################################################################
# Extract regions and submit jobs
################################################################################

REGION_JOBS=()
SUBMITTED=0
SKIPPED=0

# Extract regions from config
while read -r REGION; do
    # Skip empty lines
    [ -z "$REGION" ] && continue
    
    # Check if already processed (check for final genotype files)
    CHROM=$(echo "$REGION" | rev | cut -d'_' -f3- | rev)
    
    # Check if any sample has completed this region
    COMPLETED=true
    while read -r SAMPLE; do
        [ -z "$SAMPLE" ] && continue
        GENOTYPE_FILE="${OUTPUT_DIR}/cosigt/${SAMPLE}/${CHROM}/${REGION}/${REGION}.cosigt_genotype.tsv"
        if [ ! -f "$GENOTYPE_FILE" ]; then
            COMPLETED=false
            break
        fi
    done < <(awk '/^samples:/{flag=1;next}/^[a-zA-Z]/{flag=0}flag' "$CONFIG" | \
        sed 's/^[[:space:]]*-[[:space:]]*//' | tr -d "'" | tr -d '"' | grep -v '^[[:space:]]*$')
    
    if [ "$COMPLETED" = true ]; then
        echo "[${REGION}] Already processed, skipping"
        ((SKIPPED++))
        continue
    fi
    
    # Submit job
    JOB_ID=$(sbatch --parsable \
        --job-name=cosigt_region_${REGION} \
        --output=logs/region_${REGION}_%j.out \
        --error=logs/region_${REGION}_%j.err \
        --time=${TIME} \
        --cpus-per-task=${THREADS} \
        --mem=${MEMORY} \
        --partition=${PARTITION} \
        --wrap="module load singularity 2>/dev/null || true; bash src/run_cosigt_region.sh ${REGION} ${CONFIG} ${THREADS} ${CACHE_DIR} ${TMP_DIR}")
    
    REGION_JOBS+=("$JOB_ID")
    echo "[${REGION}] Submitted (Job ID: ${JOB_ID})"
    ((SUBMITTED++))
    
done < <(awk '/^regions:/{flag=1;next}/^[a-zA-Z]/{flag=0}flag' "$CONFIG" | \
    sed 's/^[[:space:]]*-[[:space:]]*//' | tr -d "'" | tr -d '"' | grep -v '^[[:space:]]*$')

################################################################################
# Summary
################################################################################

echo ""
echo "==================================================================="
echo "✓ Submission complete!"
echo "==================================================================="
echo "Jobs submitted: ${SUBMITTED}"
echo "Jobs skipped:   ${SKIPPED}"
echo ""

