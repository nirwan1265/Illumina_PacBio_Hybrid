#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-hybridseq}"
TAG="${2:-smoke}"
THREADS="${3:-2}"
PREFIX_DIR="${4:-.conda_envs}"
PREFIX_PATH="$PWD/$PREFIX_DIR/$ENV_NAME"

if [[ -d "$PREFIX_PATH/conda-meta" ]]; then
  CONDA_RUN=(conda run -p "$PREFIX_PATH")
  PY="$PREFIX_PATH/bin/python"
  QUAST="$PREFIX_PATH/bin/quast.py"
else
  CONDA_RUN=(conda run -n "$ENV_NAME")
  PY="python"
  QUAST="$(command -v quast.py || echo quast.py)"
fi

# If QUAST install path includes spaces, run from a temp copy without spaces.
SITE_PACKAGES="$("$PY" - <<'PY'
import site
print(site.getsitepackages()[0])
PY
)"
SITE_PACKAGES_REAL="$("$PY" - <<'PY'
import site, os
print(os.path.realpath(site.getsitepackages()[0]))
PY
)"
if [[ "$QUAST" == *" "* || "$SITE_PACKAGES" == *" "* || "$SITE_PACKAGES_REAL" == *" "* ]]; then
  QUAST_TMP="${TMPDIR:-/tmp}/quast_pkg"
  mkdir -p "$QUAST_TMP"
  cp "$QUAST" "$QUAST_TMP/quast.py"
  if [[ -d "$SITE_PACKAGES/quast_libs" ]]; then
    rm -rf "$QUAST_TMP/quast_libs"
    cp -R "$SITE_PACKAGES/quast_libs" "$QUAST_TMP/"
  fi
  QUAST="$QUAST_TMP/quast.py"
  QUAST_PYTHONPATH="$QUAST_TMP"
fi

# Small test coverage
ILL_COVS=(10)
PB_COVS=(5)

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH." >&2
  exit 1
fi

# Minimal directory setup
mkdir -p 00_ref 01_simref 02_reads 03_assemblies 04_eval 05_summary
mkdir -p "02_reads/${TAG}/illumina" "02_reads/${TAG}/pacbio_CLR" "02_reads/${TAG}/pacbio_HIFI"

# Expect a reference in 00_ref/ecoli.fa
if [[ ! -f 00_ref/ecoli.fa ]]; then
  echo "ERROR: 00_ref/ecoli.fa not found. Put a FASTA there first." >&2
  exit 1
fi

# Make repeat genome (light)
"${CONDA_RUN[@]}" scripts/mac/01_make_tandem_repeats.py \
  --in_fa 00_ref/ecoli.fa \
  --out_fa 01_simref/${TAG}.fa \
  --n_events 2 --seg_len 500 --copies 3 --seed 1

# Simulate reads (small)
for ic in "${ILL_COVS[@]}"; do
  outp="02_reads/${TAG}/illumina/${TAG}.ill_cov${ic}_"
  if [[ ! -f "${outp}1.fq.gz" ]]; then
    "${CONDA_RUN[@]}" scripts/mac/02_sim_illumina_art.sh "01_simref/${TAG}.fa" "$outp" "$ic" 150 350 50
  fi
done

for pc in "${PB_COVS[@]}"; do
  outd="02_reads/${TAG}/pacbio_CLR/cov${pc}"
  if ! ls "$outd"/pbCLR_cov${pc}*.fastq.gz "$outd"/pbCLR_cov${pc}*.fq.gz >/dev/null 2>&1; then
    "${CONDA_RUN[@]}" scripts/mac/03_sim_pacbio.sh "01_simref/${TAG}.fa" "$outd" "$pc" CLR 1
  fi
  outd2="02_reads/${TAG}/pacbio_HIFI/cov${pc}"
  if ! ls "$outd2"/pbHIFI_cov${pc}*.fastq.gz "$outd2"/pbHIFI_cov${pc}*.fq.gz >/dev/null 2>&1; then
    "${CONDA_RUN[@]}" scripts/mac/03_sim_pacbio.sh "01_simref/${TAG}.fa" "$outd2" "$pc" HIFI 1
  fi
done

# Assemble (single pair)
for ic in "${ILL_COVS[@]}"; do
  R1="02_reads/${TAG}/illumina/${TAG}.ill_cov${ic}_1.fq.gz"
  R2="02_reads/${TAG}/illumina/${TAG}.ill_cov${ic}_2.fq.gz"

  for pc in "${PB_COVS[@]}"; do
    PB_CLR="$(ls -1 02_reads/${TAG}/pacbio_CLR/cov${pc}/pbCLR_cov${pc}*.fastq.gz 02_reads/${TAG}/pacbio_CLR/cov${pc}/pbCLR_cov${pc}*.fq.gz 2>/dev/null | head -n 1 || true)"
    OUT_CLR="03_assemblies/${TAG}/CLR/ill${ic}_pb${pc}"
    mkdir -p "$OUT_CLR"
    if [[ -n "$PB_CLR" && ! -f "${OUT_CLR}/assembly.fasta" ]]; then
      "${CONDA_RUN[@]}" unicycler -1 "$R1" -2 "$R2" -l "$PB_CLR" -o "$OUT_CLR" -t "$THREADS" --mode conservative --kmers 21,33,55 --no_miniasm --no_long_read_alignment --no_simple_bridges
    fi

    PB_HIFI="$(ls -1 02_reads/${TAG}/pacbio_HIFI/cov${pc}/pbHIFI_cov${pc}*.fastq.gz 02_reads/${TAG}/pacbio_HIFI/cov${pc}/pbHIFI_cov${pc}*.fq.gz 2>/dev/null | head -n 1 || true)"
    OUT_HIFI="03_assemblies/${TAG}/HIFI/ill${ic}_pb${pc}"
    mkdir -p "$OUT_HIFI"
    if [[ -n "$PB_HIFI" && ! -f "${OUT_HIFI}/assembly.fasta" ]]; then
      "${CONDA_RUN[@]}" unicycler -1 "$R1" -2 "$R2" -l "$PB_HIFI" -o "$OUT_HIFI" -t "$THREADS" --mode conservative --kmers 21,33,55 --no_miniasm --no_long_read_alignment --no_simple_bridges
    fi
  done
done

# QUAST on both
for MODE in CLR HIFI; do
  ASM="03_assemblies/${TAG}/${MODE}/ill10_pb5/assembly.fasta"
  OUT="04_eval/${TAG}/${MODE}/ill10_pb5"
  if [[ -f "$ASM" && ! -f "${OUT}/report.tsv" ]]; then
    if [[ -n "${QUAST_PYTHONPATH:-}" ]]; then
      PYTHONPATH="$QUAST_PYTHONPATH" "$PY" "$QUAST" -r "01_simref/${TAG}.fa" -o "$OUT" --threads 1 "$ASM"
    else
      "$PY" "$QUAST" -r "01_simref/${TAG}.fa" -o "$OUT" --threads 1 "$ASM"
    fi
  fi
done

# Summarize
"${CONDA_RUN[@]}" scripts/mac/07_summarize_quast.R "04_eval/${TAG}/CLR"  "05_summary/${TAG}.CLR.quast_summary.csv" || true
"${CONDA_RUN[@]}" scripts/mac/07_summarize_quast.R "04_eval/${TAG}/HIFI" "05_summary/${TAG}.HIFI.quast_summary.csv" || true

echo "Smoke test complete. Check 05_summary/${TAG}.*.quast_summary.csv"
