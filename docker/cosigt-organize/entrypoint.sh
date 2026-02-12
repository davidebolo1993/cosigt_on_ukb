#!/bin/bash
set -euo pipefail

cd /work 2>/dev/null || cd /

# Run organize.py directly from the image
exec python /usr/local/bin/organize.py "$@"
