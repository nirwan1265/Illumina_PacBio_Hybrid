#!/usr/bin/env bash
set -euo pipefail

TSV="${1:-data/markers_curated.tsv}"
OUT_FA="${2:-data/markers_curated.fa}"

if [[ ! -f "$TSV" ]]; then
  echo "ERROR: marker TSV not found: $TSV" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_FA")"
: > "$OUT_FA"

# Fetch one record from NCBI E-utilities. If start/end set, use seq_start/seq_stop.
fetch_one() {
  local name="$1" acc="$2" start="$3" end="$4" strand="$5"
  local base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  local url="${base_url}?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
  if [[ -n "$start" && -n "$end" ]]; then
    url="${url}&seq_start=${start}&seq_stop=${end}"
  fi

  # Fetch and normalize header
  local tmp
  tmp="$(mktemp)"
  curl -fsSL "$url" > "$tmp"

  # Replace FASTA header with a stable name
  local hdr
  hdr=">${name}|${acc}"
  if [[ -n "$start" && -n "$end" ]]; then
    hdr+="|${start}-${end}"
  fi
  if [[ -n "$strand" ]]; then
    hdr+="|strand=${strand}"
  fi

  # Replace header line
  { echo "$hdr"; grep -v '^>' "$tmp"; } > "${tmp}.norm"

  # Reverse-complement if needed
  if [[ "$strand" == "-" ]]; then
    if command -v seqkit >/dev/null 2>&1; then
      seqkit seq -r -p "${tmp}.norm" > "${tmp}.rc"
      mv "${tmp}.rc" "${tmp}.norm"
    else
      echo "WARN: seqkit not found; skipping reverse-complement for $name" >&2
    fi
  fi

  cat "${tmp}.norm" >> "$OUT_FA"
  echo >> "$OUT_FA"
  rm -f "$tmp" "${tmp}.norm"
}

# Read TSV
# name accession start end strand notes
while IFS=$'\t' read -r name acc start end strand notes; do
  # skip header or empty/comment lines
  [[ -z "$name" ]] && continue
  [[ "$name" == "name" ]] && continue
  [[ "$name" =~ ^# ]] && continue

  fetch_one "$name" "$acc" "$start" "$end" "$strand"
  # polite delay
  sleep 0.2
done < "$TSV"

echo "Wrote curated marker FASTA: $OUT_FA"
