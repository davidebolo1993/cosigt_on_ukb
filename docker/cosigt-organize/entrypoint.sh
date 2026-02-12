#!/bin/bash
set -euo pipefail

# Clone latest scripts from GitHub
REPO_URL="${COSIGT_REPO_URL:-https://github.com/davidebolo1993/cosigt_on_ukbb.git}"
BRANCH="${COSIGT_BRANCH:-main}"

echo "Cloning COSIGT scripts from $REPO_URL (branch: $BRANCH)..."
git clone -b "$BRANCH" --depth 1 "$REPO_URL" /cosigt-scripts

# Run organize.py with passed arguments
python /cosigt-scripts/src/organize.py "$@"
