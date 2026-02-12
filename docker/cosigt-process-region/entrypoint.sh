#!/bin/bash
set -euo pipefail

cd /work 2>/dev/null || cd /

# Run process script directly from the image
exec bash /usr/local/bin/process_region_docker.sh "$@"
