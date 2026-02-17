#!/bin/bash
# Generate FASTAs for all 100 samples

set -euo pipefail

VCF_GZ="05_summary/bc2s3_controlled.vcf.gz"
REF="05_summary/bc2s3_fastas/reference.fa"
OUT_DIR="05_summary/bc2s3_fastas"

echo "=== Generating FASTAs for all samples ==="

# Get all sample names
SAMPLES=($(bcftools query -l "$VCF_GZ"))
N=${#SAMPLES[@]}

echo "Total samples: $N"

for i in "${!SAMPLES[@]}"; do
    S="${SAMPLES[$i]}"
    OUT="${OUT_DIR}/${S}.fa"

    if [[ -f "$OUT" ]]; then
        # Check if file is complete (should be ~100Mb)
        SIZE=$(stat -f%z "$OUT" 2>/dev/null || stat -c%s "$OUT" 2>/dev/null || echo 0)
        if [[ $SIZE -gt 100000000 ]]; then
            echo "[$((i+1))/$N] $S - exists (${SIZE} bytes), skipping"
            continue
        fi
    fi

    echo "[$((i+1))/$N] $S"

    # Extract sample VCF
    TMP="/tmp/${S}_tmp.vcf.gz"
    bcftools view -s "$S" "$VCF_GZ" -Oz -o "$TMP" 2>/dev/null
    bcftools index -f "$TMP" 2>/dev/null

    # Generate haplotypes
    bcftools consensus -f "$REF" -H 1 "$TMP" 2>/dev/null > "/tmp/${S}_h1.fa"
    bcftools consensus -f "$REF" -H 2 "$TMP" 2>/dev/null > "/tmp/${S}_h2.fa"

    # Combine
    {
        echo ">${S}_hap1"
        grep -v "^>" "/tmp/${S}_h1.fa"
        echo ">${S}_hap2"
        grep -v "^>" "/tmp/${S}_h2.fa"
    } > "$OUT"

    # Cleanup
    rm -f "$TMP" "${TMP}.csi" "/tmp/${S}_h1.fa" "/tmp/${S}_h2.fa"
done

echo ""
echo "=== Done ==="
echo "Generated $N sample FASTAs in $OUT_DIR"
