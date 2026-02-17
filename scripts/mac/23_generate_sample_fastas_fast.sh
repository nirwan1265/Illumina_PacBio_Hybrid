#!/bin/bash
# Fast generation of reference and sample FASTAs using bcftools consensus
#
# Usage: ./23_generate_sample_fastas_fast.sh [options]

set -euo pipefail

VCF="05_summary/bc2s3_controlled.vcf"
OUT_DIR="05_summary/bc2s3_fastas"
GENOME_SIZE=100000000
SEED=42
SAMPLES=""  # Empty = all samples

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --vcf <path>          VCF file [default: $VCF]"
    echo "  --out-dir <path>      Output directory [default: $OUT_DIR]"
    echo "  --genome-size <int>   Genome size [default: $GENOME_SIZE]"
    echo "  --samples <csv>       Specific samples (default: all)"
    echo "  --seed <int>          Random seed [default: $SEED]"
    echo "  --help                Show help"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf) VCF="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --genome-size) GENOME_SIZE="$2"; shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --help) usage ;;
        *) echo "Unknown: $1"; usage ;;
    esac
done

echo "=== Fast Sample FASTA Generation ==="
echo "VCF:         $VCF"
echo "Output:      $OUT_DIR"
echo "Genome size: $GENOME_SIZE"
echo ""

mkdir -p "$OUT_DIR"

REF="${OUT_DIR}/reference.fa"

# Step 1: Generate random reference (fast with awk)
if [[ ! -f "$REF" ]]; then
    echo "Generating reference sequence..."
    awk -v size="$GENOME_SIZE" -v seed="$SEED" '
    BEGIN {
        srand(seed)
        bases = "ACGT"
        print ">chr1"
        for (i = 1; i <= size; i++) {
            printf "%s", substr(bases, int(rand()*4)+1, 1)
            if (i % 80 == 0) printf "\n"
        }
        if (size % 80 != 0) printf "\n"
    }' > "$REF"
    echo "  Created reference: $REF"
else
    echo "  Reference exists: $REF"
fi

# Index reference
[[ ! -f "${REF}.fai" ]] && samtools faidx "$REF"

# Step 2: Compress and index VCF if needed
VCF_GZ="${VCF}.gz"
if [[ ! -f "$VCF_GZ" ]]; then
    echo "Compressing VCF..."
    bgzip -c "$VCF" > "$VCF_GZ"
    bcftools index "$VCF_GZ"
fi

# Step 3: Get sample list
if [[ -n "$SAMPLES" ]]; then
    IFS=',' read -ra SAMPLE_ARR <<< "$SAMPLES"
else
    mapfile -t SAMPLE_ARR < <(bcftools query -l "$VCF_GZ")
fi

N_SAMPLES=${#SAMPLE_ARR[@]}
echo "Processing $N_SAMPLES samples..."

# Step 4: Generate per-sample FASTAs using bcftools consensus
for i in "${!SAMPLE_ARR[@]}"; do
    SAMPLE="${SAMPLE_ARR[$i]}"
    SAMPLE_FA="${OUT_DIR}/${SAMPLE}.fa"

    if [[ -f "$SAMPLE_FA" ]]; then
        echo "  [$((i+1))/$N_SAMPLES] $SAMPLE - exists, skipping"
        continue
    fi

    echo "  [$((i+1))/$N_SAMPLES] $SAMPLE"

    # Create temporary per-sample VCF for haplotype 1 (first allele)
    TMP_HAP1=$(mktemp).vcf.gz
    TMP_HAP2=$(mktemp).vcf.gz

    # Extract sample and set to haplotype 1 (0|0 -> 0, 0|1 -> 0, 1|1 -> 1)
    bcftools view -s "$SAMPLE" "$VCF_GZ" 2>/dev/null | \
        awk -F'\t' -v OFS='\t' '
        /^#/ { print; next }
        {
            gt = $10
            split(gt, parts, ":")
            split(parts[1], alleles, /[|\/]/)
            if (alleles[1] == "1") parts[1] = "1"
            else parts[1] = "0"
            $10 = parts[1]
            for (i=2; i<=length(parts); i++) $10 = $10 ":" parts[i]
            print
        }' | bgzip > "$TMP_HAP1"
    bcftools index "$TMP_HAP1"

    # Extract sample and set to haplotype 2 (0|0 -> 0, 0|1 -> 1, 1|1 -> 1)
    bcftools view -s "$SAMPLE" "$VCF_GZ" 2>/dev/null | \
        awk -F'\t' -v OFS='\t' '
        /^#/ { print; next }
        {
            gt = $10
            split(gt, parts, ":")
            split(parts[1], alleles, /[|\/]/)
            if (alleles[2] == "1" || (length(alleles) == 1 && alleles[1] == "1")) parts[1] = "1"
            else if (length(alleles) > 1 && alleles[2] == "1") parts[1] = "1"
            else parts[1] = "0"
            $10 = parts[1]
            for (i=2; i<=length(parts); i++) $10 = $10 ":" parts[i]
            print
        }' | bgzip > "$TMP_HAP2"
    bcftools index "$TMP_HAP2"

    # Generate haplotype sequences
    HAP1_SEQ=$(bcftools consensus -f "$REF" "$TMP_HAP1" 2>/dev/null)
    HAP2_SEQ=$(bcftools consensus -f "$REF" "$TMP_HAP2" 2>/dev/null)

    # Write combined FASTA
    {
        echo ">${SAMPLE}_hap1"
        echo "$HAP1_SEQ" | grep -v "^>" | tr -d '\n' | fold -w 80
        echo ""
        echo ">${SAMPLE}_hap2"
        echo "$HAP2_SEQ" | grep -v "^>" | tr -d '\n' | fold -w 80
        echo ""
    } > "$SAMPLE_FA"

    # Cleanup
    rm -f "$TMP_HAP1" "$TMP_HAP2" "${TMP_HAP1}.csi" "${TMP_HAP2}.csi"
done

echo ""
echo "=== Done ==="
echo "Reference: $REF"
echo "Samples:   ${OUT_DIR}/sample*.fa"
