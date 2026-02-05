#!/usr/bin/env bash
set -euo pipefail

FA="${1:?usage: sim_pacbio.sh <ref.fa> <outdir> <cov> <CLR|HIFI> [seed]}"
OUTDIR="${2:?}"
COV="${3:?}"
TYPE="${4:?}"
SEED="${5:-1}"

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
  echo "ERROR: need pbsim3 or pbsim installed" >&2
  exit 1
fi

ls -1 *.fastq 2>/dev/null | xargs -I{} pigz -f "{}" || true
