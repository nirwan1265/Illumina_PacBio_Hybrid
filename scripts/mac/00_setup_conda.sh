#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-hybridseq}"
PREFIX_DIR="${2:-.conda_envs}"
PREFIX_PATH="$PWD/$PREFIX_DIR/$ENV_NAME"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH. Install Miniconda/Anaconda first." >&2
  exit 1
fi

# Create env if it doesn't exist (prefix-based for mac portability)
if [[ ! -d "$PREFIX_PATH/conda-meta" ]]; then
  conda create -y -p "$PREFIX_PATH" -c conda-forge -c bioconda \
    python=3.10 \
    art \
    pbsim \
    pbsim3 \
    unicycler \
    quast \
    samtools \
    seqkit \
    pigz \
    ncbi-datasets-cli
else
  echo "Conda env prefix '$PREFIX_PATH' already exists. Skipping create."
fi

# Show versions
conda run -p "$PREFIX_PATH" art_illumina --help | head -n 2 || true
conda run -p "$PREFIX_PATH" pbsim --help | head -n 2 || true
conda run -p "$PREFIX_PATH" pbsim3 --help | head -n 2 || true
conda run -p "$PREFIX_PATH" unicycler --help | head -n 2 || true
conda run -p "$PREFIX_PATH" quast.py --help | head -n 2 || true
conda run -p "$PREFIX_PATH" seqkit --help | head -n 2 || true



### conda activate "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/Illumina_PacBio_Hybrid/.conda_envs/hybridseq"