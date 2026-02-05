#!/usr/bin/env bash
set -euo pipefail

TAG="${1:?usage: run_quast_both.sh <tag>}"
REF="01_simref/${TAG}.fa"

ILL_COVS=(10 20 30 40 60 80)
PB_COVS=(5 10 15 20 30 40)

for MODE in CLR HIFI; do
  OUTROOT="04_eval/${TAG}/${MODE}"
  mkdir -p "$OUTROOT"

  for ic in "${ILL_COVS[@]}"; do
    for pc in "${PB_COVS[@]}"; do
      ASM="03_assemblies/${TAG}/${MODE}/ill${ic}_pb${pc}/assembly.fasta"
      OUT="${OUTROOT}/ill${ic}_pb${pc}"

      if [[ ! -f "$ASM" ]]; then
        continue
      fi
      if [[ -f "${OUT}/report.tsv" ]]; then
        continue
      fi

      quast.py -r "$REF" -o "$OUT" --threads 4 "$ASM"
    done
  done
done
