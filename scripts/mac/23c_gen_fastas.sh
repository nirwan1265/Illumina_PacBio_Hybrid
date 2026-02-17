#!/bin/bash
# Generate sample FASTAs - simplified version

VCF_GZ="05_summary/bc2s3_controlled.vcf.gz"
REF="05_summary/bc2s3_fastas/reference.fa"
OUT_DIR="05_summary/bc2s3_fastas"

# Get samples from argument or default
SAMPLES="${1:-sample1,sample2,sample3}"

echo "Generating FASTAs for: $SAMPLES"

# Convert comma-separated to array
IFS=',' read -ra SARR <<< "$SAMPLES"

for S in "${SARR[@]}"; do
    echo "Processing $S..."

    TMP="/tmp/${S}_tmp.vcf.gz"
    OUT="${OUT_DIR}/${S}.fa"

    # Extract sample (keep all sites including ref-only for proper consensus)
    bcftools view -s "$S" "$VCF_GZ" 2>/dev/null | bgzip > "$TMP"
    bcftools index -f "$TMP" 2>/dev/null

    # Generate haplotypes
    bcftools consensus -f "$REF" -H 1 "$TMP" 2>/dev/null > "/tmp/${S}_h1.fa"
    bcftools consensus -f "$REF" -H 2 "$TMP" 2>/dev/null > "/tmp/${S}_h2.fa"

    # Combine
    echo ">${S}_hap1" > "$OUT"
    grep -v "^>" "/tmp/${S}_h1.fa" >> "$OUT"
    echo ">${S}_hap2" >> "$OUT"
    grep -v "^>" "/tmp/${S}_h2.fa" >> "$OUT"

    # Cleanup
    rm -f "$TMP" "${TMP}.csi" "/tmp/${S}_h1.fa" "/tmp/${S}_h2.fa"

    echo "  Created: $OUT"
done

echo "Done!"
