#!/bin/bash
set -euo pipefail

REPO_URL="${COSIGT_REPO_URL:-https://github.com/davidebolo1993/cosigt_on_ukbb.git}"
BRANCH="${COSIGT_BRANCH:-main}"

echo "Cloning COSIGT scripts from $REPO_URL (branch: $BRANCH)..."
git clone -b "$BRANCH" --depth 1 "$REPO_URL" /cosigt-scripts

#Run preprocess_sample.sh with passed arguments
bash /cosigt-scripts/src/preprocess_sample.sh "$@"

