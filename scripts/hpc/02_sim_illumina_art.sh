#!/usr/bin/env bash
set -euo pipefail

FA="${1:?usage: sim_illumina_art.sh <ref.fa> <outprefix> <cov> [readlen] [ins] [sd]}"
OUTPREFIX="${2:?}"
COV="${3:?}"
READLEN="${4:-150}"
INS="${5:-350}"
SD="${6:-50}"

if ! command -v art_illumina >/dev/null 2>&1; then
  echo "ERROR: art_illumina not found in PATH" >&2
  exit 1
fi

art_illumina \
  -ss HS25 \
  -i "$FA" \
  -p \
  -l "$READLEN" \
  -f "$COV" \
  -m "$INS" \
  -s "$SD" \
  -o "$OUTPREFIX"

# ART outputs: ${OUTPREFIX}1.fq and ${OUTPREFIX}2.fq
pigz -f "${OUTPREFIX}1.fq" "${OUTPREFIX}2.fq"
