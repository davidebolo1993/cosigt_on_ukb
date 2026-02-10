#!/bin/bash

################################################################################
# Submit Sample Preprocessing Jobs
#
# Usage: ./submit_sample_preprocessing.sh [options]
#
# Submits one preprocessing job per sample listed in config.yaml
################################################################################

# Default values
CONFIG="config/config.yaml"
CACHE_DIR="$PWD/singularity_cache"
TMP_DIR="/tmp"
THREADS=4
MEMORY="4G"
TIME="00:30:00"
PARTITION="cpuq"

# Usage function
show_usage() {
    cat << EOF
Usage: $0 [options]

Options:
    -c, --config FILE       Config file [default: config/config.yaml]
    --cache-dir DIR         Singularity cache directory [default: \$PWD/singularity_cache]
    --tmp-dir DIR           Singularity tmp directory [default: /tmp]
    -t, --threads N         Threads per job [default: 4]
    -m, --memory MEM        Memory per job [default: 4G]
    --time TIME             Time limit per job [default: 00:30:00]
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
echo "Submitting Sample Preprocessing Jobs"
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
# Extract samples and submit jobs
################################################################################

SAMPLE_JOBS=()
SUBMITTED=0
SKIPPED=0

# Extract samples from config
while read -r SAMPLE; do
    # Skip empty lines
    [ -z "$SAMPLE" ] && continue
    
    # Check if already processed
    OUTPUT_DIR=$(grep "^output:" "$CONFIG" | sed "s/^output:[[:space:]]*//" | tr -d "'" | tr -d '"')
    UNMAPPED_FILE="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/unmapped.fasta.gz"
    
    if [ -f "$UNMAPPED_FILE" ]; then
        echo "[${SAMPLE}] Already processed, skipping"
        ((SKIPPED++))
        continue
    fi
    
    # Submit job
    JOB_ID=$(sbatch --parsable \
        --job-name=cosigt_sample_${SAMPLE} \
        --output=logs/preprocess_${SAMPLE}_%j.out \
        --error=logs/preprocess_${SAMPLE}_%j.err \
        --time=${TIME} \
        --cpus-per-task=${THREADS} \
        --mem=${MEMORY} \
        --partition=${PARTITION} \
        --wrap="module load singularity 2>/dev/null || true; bash src/preprocess_sample.sh ${SAMPLE} ${CONFIG} ${THREADS} ${CACHE_DIR} ${TMP_DIR}")
    
    SAMPLE_JOBS+=("$JOB_ID")
    echo "[${SAMPLE}] Submitted (Job ID: ${JOB_ID})"
    ((SUBMITTED++))
    
done < <(awk '/^samples:/{flag=1;next}/^[a-zA-Z]/{flag=0}flag' "$CONFIG" | \
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

