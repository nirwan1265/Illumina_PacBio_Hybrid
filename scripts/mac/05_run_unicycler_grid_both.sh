#!/usr/bin/env bash
set -euo pipefail

TAG="${1:?usage: run_unicycler_grid_both.sh <tag> [threads]}"
THREADS="${2:-4}"

ILL_COVS=(10 20 30 40 60 80)
PB_COVS=(5 10 15 20 30 40)

for ic in "${ILL_COVS[@]}"; do
  R1="02_reads/${TAG}/illumina/${TAG}.ill_cov${ic}_1.fq.gz"
  R2="02_reads/${TAG}/illumina/${TAG}.ill_cov${ic}_2.fq.gz"

  for pc in "${PB_COVS[@]}"; do
    PB_DIR_CLR="02_reads/${TAG}/pacbio_CLR/cov${pc}"
    PB_CLR="$(ls -1 ${PB_DIR_CLR}/pbCLR_cov${pc}*.fastq.gz ${PB_DIR_CLR}/pbCLR_cov${pc}*.fq.gz 2>/dev/null | head -n 1 || true)"

    OUT_CLR="03_assemblies/${TAG}/CLR/ill${ic}_pb${pc}"
    mkdir -p "$OUT_CLR"
    if [[ -n "$PB_CLR" && ! -f "${OUT_CLR}/assembly.fasta" ]]; then
      echo "Assembling (CLR) ${TAG}: Ill ${ic}x + PB_CLR ${pc}x"
      unicycler -1 "$R1" -2 "$R2" -l "$PB_CLR" -o "$OUT_CLR" -t "$THREADS" --mode normal
    fi

    PB_DIR_HIFI="02_reads/${TAG}/pacbio_HIFI/cov${pc}"
    PB_HIFI="$(ls -1 ${PB_DIR_HIFI}/pbHIFI_cov${pc}*.fastq.gz ${PB_DIR_HIFI}/pbHIFI_cov${pc}*.fq.gz 2>/dev/null | head -n 1 || true)"

    OUT_HIFI="03_assemblies/${TAG}/HIFI/ill${ic}_pb${pc}"
    mkdir -p "$OUT_HIFI"
    if [[ -n "$PB_HIFI" && ! -f "${OUT_HIFI}/assembly.fasta" ]]; then
      echo "Assembling (HIFI) ${TAG}: Ill ${ic}x + PB_HIFI ${pc}x"
      unicycler -1 "$R1" -2 "$R2" -l "$PB_HIFI" -o "$OUT_HIFI" -t "$THREADS" --mode normal
    fi
  done
done
