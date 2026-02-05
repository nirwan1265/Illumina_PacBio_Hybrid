#!/usr/bin/env bash
set -euo pipefail

TAG="${1:?usage: run_quast_both.sh <tag>}"
REF="01_simref/${TAG}.fa"
PY="$(command -v python)"
QUAST="$(command -v quast.py)"

if [[ -z "$PY" || -z "$QUAST" ]]; then
  echo "ERROR: python or quast.py not found in PATH. Activate your env first." >&2
  exit 1
fi

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

      if [[ -n "${QUAST_PYTHONPATH:-}" ]]; then
        PYTHONPATH="$QUAST_PYTHONPATH" "$PY" "$QUAST" -r "$REF" -o "$OUT" --threads 2 "$ASM"
      else
        "$PY" "$QUAST" -r "$REF" -o "$OUT" --threads 2 "$ASM"
      fi
    done
  done
done
