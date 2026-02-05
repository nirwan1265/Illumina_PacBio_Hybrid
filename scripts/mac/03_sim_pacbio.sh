#!/usr/bin/env bash
set -euo pipefail

FA="${1:?usage: sim_pacbio.sh <ref.fa> <outdir> <cov> <CLR|HIFI> [seed]}"
OUTDIR="${2:?}"
COV="${3:?}"
TYPE="${4:?}"
SEED="${5:-1}"

# Resolve reference to absolute path before changing dirs
FA="$(cd "$(dirname "$FA")" && pwd)/$(basename "$FA")"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

if command -v pbsim3 >/dev/null 2>&1; then
  if [[ "$TYPE" == "CLR" ]]; then
    echo "[pbsim3] Simulating CLR at ${COV}x"
    pbsim3 --depth "$COV" --seed "$SEED" --prefix "pbCLR_cov${COV}" "$FA"
  elif [[ "$TYPE" == "HIFI" ]]; then
    echo "[pbsim3] Simulating HiFi/CCS at ${COV}x"
    pbsim3 --depth "$COV" --seed "$SEED" --prefix "pbHIFI_cov${COV}" "$FA"
  else
    echo "TYPE must be CLR or HIFI" >&2
    exit 1
  fi

elif command -v pbsim >/dev/null 2>&1; then
  # Check if this pbsim supports --data-type
  if pbsim 2>/dev/null | grep -q -- "--data-type"; then
    if [[ "$TYPE" == "CLR" ]]; then
      echo "[pbsim] Simulating CLR at ${COV}x"
      pbsim --depth "$COV" --seed "$SEED" --data-type CLR --prefix "pbCLR_cov${COV}" "$FA"
    elif [[ "$TYPE" == "HIFI" ]]; then
      echo "[pbsim] Simulating CCS/HiFi at ${COV}x"
      pbsim --depth "$COV" --seed "$SEED" --data-type CCS --prefix "pbHIFI_cov${COV}" "$FA"
    else
      echo "TYPE must be CLR or HIFI" >&2
      exit 1
    fi
  else
    # Older pbsim: use bundled model files with pass-num to approximate CLR vs HiFi
    PB_SIM_ROOT="$(cd "$(dirname "$(command -v pbsim)")/.." && pwd)"
    MODEL_DIR="${PB_SIM_ROOT}/data"
    ERR_MODEL="${MODEL_DIR}/ERRHMM-RSII.model"
    QSH_MODEL="${MODEL_DIR}/QSHMM-RSII.model"

    if [[ ! -f "$ERR_MODEL" || ! -f "$QSH_MODEL" ]]; then
      echo "ERROR: pbsim model files not found in ${MODEL_DIR}" >&2
      exit 1
    fi

    if [[ "$TYPE" == "CLR" ]]; then
      echo "[pbsim] Simulating CLR (errhmm, pass-num 1) at ${COV}x"
      pbsim --strategy wgs --method errhmm --errhmm "$ERR_MODEL" --depth "$COV" --seed "$SEED" \
        --pass-num 1 --accuracy-mean 0.85 \
        --prefix "pbCLR_cov${COV}" --genome "$FA"
    elif [[ "$TYPE" == "HIFI" ]]; then
      echo "[pbsim] Simulating HiFi (qshmm, pass-num 10) at ${COV}x"
      pbsim --strategy wgs --method qshmm --qshmm "$QSH_MODEL" --depth "$COV" --seed "$SEED" \
        --pass-num 10 --accuracy-mean 0.99 \
        --prefix "pbHIFI_cov${COV}" --genome "$FA"
    else
      echo "TYPE must be CLR or HIFI" >&2
      exit 1
    fi
  fi
else
  echo "ERROR: need pbsim3 or pbsim installed" >&2
  exit 1
fi

# Convert BAM outputs (if any) to FASTQ, then gzip
if command -v samtools >/dev/null 2>&1; then
  for b in *.bam; do
    [[ -f "$b" ]] || break
    fq="${b%.bam}.fastq"
    samtools fastq "$b" > "$fq"
  done
fi

ls -1 *.fastq *.fq 2>/dev/null | xargs -I{} pigz -f "{}" || true
