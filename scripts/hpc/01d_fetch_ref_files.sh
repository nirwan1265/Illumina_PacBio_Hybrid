#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-inst/extdata/ref_files}"
mkdir -p "$OUTDIR"

fetch_ncbi() {
  local acc="$1"
  local out="$2"
  local base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  local url="${base_url}?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
  echo "[fetch] ${acc} -> ${out}"
  curl -fsSL "$url" > "$out"
}

# E. coli K-12 MG1655 chromosome
fetch_ncbi "NC_000913.3" "${OUTDIR}/ecoli_k12_mg1655_chr.fa"

# Human GRCh38 chromosome 1
fetch_ncbi "NC_000001.11" "${OUTDIR}/human_grch38_chr1.fa"

# Maize (Zea mays) chromosome 1 (RefSeq)
fetch_ncbi "NC_050096.1" "${OUTDIR}/maize_b73_refseq_chr1.fa"

echo "Done. Files in: ${OUTDIR}"
