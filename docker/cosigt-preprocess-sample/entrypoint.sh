#!/bin/bash
set -euo pipefail

cd /work 2>/dev/null || cd /

# Run preprocess script directly from the image
exec bash /usr/local/bin/preprocess_sample_docker.sh "$@"
