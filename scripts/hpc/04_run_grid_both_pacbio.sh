#!/usr/bin/env bash
set -euo pipefail

FA="${1:?usage: run_grid_both_pacbio.sh <simref.fa> <tag>}"
TAG="${2:?}"

ILL_COVS=(10 20 30 40 60 80)
PB_CLR_COVS=(5 10 15 20 30 40)
PB_HIFI_COVS=(5 10 15 20 30 40)

mkdir -p "02_reads/${TAG}/illumina" \
         "02_reads/${TAG}/pacbio_CLR" \
         "02_reads/${TAG}/pacbio_HIFI"

for ic in "${ILL_COVS[@]}"; do
  outp="02_reads/${TAG}/illumina/${TAG}.ill_cov${ic}_"
  if [[ ! -f "${outp}1.fq.gz" ]]; then
    scripts/sim_illumina_art.sh "$FA" "$outp" "$ic" 150 350 50
  fi
done

for pc in "${PB_CLR_COVS[@]}"; do
  outd="02_reads/${TAG}/pacbio_CLR/cov${pc}"
  if ! ls "$outd"/pbCLR_cov${pc}*.fastq.gz >/dev/null 2>&1; then
    scripts/sim_pacbio.sh "$FA" "$outd" "$pc" CLR 1
  fi
done

for pc in "${PB_HIFI_COVS[@]}"; do
  outd="02_reads/${TAG}/pacbio_HIFI/cov${pc}"
  if ! ls "$outd"/pbHIFI_cov${pc}*.fastq.gz >/dev/null 2>&1; then
    scripts/sim_pacbio.sh "$FA" "$outd" "$pc" HIFI 1
  fi
done
