#!/usr/bin/env bash
set -euo pipefail

FA="${1:?usage: run_grid_both_pacbio.sh <simref.fa> <tag>}"
TAG="${2:?}"

# Mac-friendly: use smaller defaults if desired by editing these arrays
ILL_COVS=(10 20 30 40 60 80)
PB_CLR_COVS=(5 10 15 20 30 40)
PB_HIFI_COVS=(5 10 15 20 30 40)

mkdir -p "02_reads/${TAG}/illumina" \
         "02_reads/${TAG}/pacbio_CLR" \
         "02_reads/${TAG}/pacbio_HIFI"

# Illumina (ART)
for ic in "${ILL_COVS[@]}"; do
  outp="02_reads/${TAG}/illumina/${TAG}.ill_cov${ic}_"
  if [[ ! -f "${outp}1.fq.gz" ]]; then
    scripts/mac/02_sim_illumina_art.sh "$FA" "$outp" "$ic" 150 350 50
  fi
done

# PacBio CLR
for pc in "${PB_CLR_COVS[@]}"; do
  outd="02_reads/${TAG}/pacbio_CLR/cov${pc}"
  if ! ls "$outd"/pbCLR_cov${pc}*.fastq.gz "$outd"/pbCLR_cov${pc}*.fq.gz >/dev/null 2>&1; then
    scripts/mac/03_sim_pacbio.sh "$FA" "$outd" "$pc" CLR 1
  fi
done

# PacBio HiFi/CCS
for pc in "${PB_HIFI_COVS[@]}"; do
  outd="02_reads/${TAG}/pacbio_HIFI/cov${pc}"
  if ! ls "$outd"/pbHIFI_cov${pc}*.fastq.gz "$outd"/pbHIFI_cov${pc}*.fq.gz >/dev/null 2>&1; then
    scripts/mac/03_sim_pacbio.sh "$FA" "$outd" "$pc" HIFI 1
  fi
done
