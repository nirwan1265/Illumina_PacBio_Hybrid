#!/usr/bin/env bash
set -euo pipefail

echo "RNA-seq simulation stub"
echo "----------------------"
echo "This is a placeholder for transcriptome simulation."
echo "Inputs you will likely need:"
echo "  - genome FASTA"
echo "  - GFF3 or transcript annotations"
echo "  - expression profile (TPM/FPKM or counts)"
echo
if [[ $# -lt 2 ]]; then
  echo "Usage: 09_rnaseq_stub.sh <genome.fa> <gff3> [expr.tsv]" >&2
  exit 1
fi

echo "genome: $1"
echo "gff3:   $2"
if [[ $# -ge 3 ]]; then
  echo "expr:   $3"
fi

echo
 echo "Next steps (manual):"
 echo "- Choose simulator (e.g., Polyester, Flux Simulator, ART for RNA-seq)."
 echo "- Generate transcripts and simulate reads."
 echo "- Map and quantify (e.g., salmon/kallisto)."
 echo
 echo "If you want, I can implement a full RNA-seq pipeline here."
