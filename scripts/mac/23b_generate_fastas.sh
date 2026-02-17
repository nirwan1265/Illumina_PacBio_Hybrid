#!/bin/bash
# Generate sample FASTAs using bcftools consensus
# Simpler version that works

set -euo pipefail

VCF_GZ="05_summary/bc2s3_controlled.vcf.gz"
REF="05_summary/bc2s3_fastas/reference.fa"
OUT_DIR="05_summary/bc2s3_fastas"
SAMPLES="${1:-sample1,sample2,sample3}"

echo "=== Generate Sample FASTAs ==="
echo "VCF: $VCF_GZ"
echo "Ref: $REF"
echo "Samples: $SAMPLES"

IFS=',' read -ra SAMPLE_ARR <<< "$SAMPLES"
N=${#SAMPLE_ARR[@]}

for i in "${!SAMPLE_ARR[@]}"; do
    S="${SAMPLE_ARR[$i]}"
    OUT="${OUT_DIR}/${S}.fa"

    echo "[$((i+1))/$N] $S"

    # Create per-sample VCF with phased genotypes
    TMP_VCF=$(mktemp).vcf.gz

    # Extract sample, keep only variant sites
    bcftools view -s "$S" -c 1 "$VCF_GZ" 2>/dev/null | bgzip > "$TMP_VCF"
    bcftools index -f "$TMP_VCF"

    # Generate haplotype 1 (first allele)
    HAP1=$(bcftools consensus -f "$REF" -H 1 "$TMP_VCF" 2>/dev/null | grep -v "^>")

    # Generate haplotype 2 (second allele)
    HAP2=$(bcftools consensus -f "$REF" -H 2 "$TMP_VCF" 2>/dev/null | grep -v "^>")

    # Write combined FASTA
    echo ">${S}_hap1" > "$OUT"
    echo "$HAP1" | tr -d '\n' | fold -w 80 >> "$OUT"
    echo "" >> "$OUT"
    echo ">${S}_hap2" >> "$OUT"
    echo "$HAP2" | tr -d '\n' | fold -w 80 >> "$OUT"
    echo "" >> "$OUT"

    rm -f "$TMP_VCF" "${TMP_VCF}.csi"
done

echo "Done!"
