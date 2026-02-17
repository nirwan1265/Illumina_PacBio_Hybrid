#!/bin/bash
# Simulate reads, align, and call GLs with bcftools
#
# Usage: ./24_simulate_reads_bcftools.sh [options]
#
# Requires: art_illumina, bwa, samtools, bcftools (all in simitall conda)

set -euo pipefail

# Defaults
FASTA_DIR="05_summary/bc2s3_fastas"
OUT_DIR="05_summary/bc2s3_reads"
COVERAGES="10,8,4,1,0.5"
THREADS=4
SAMPLE_LIST=""  # Empty = all samples

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --fasta-dir <path>    Directory with sample FASTAs [default: $FASTA_DIR]"
    echo "  --out-dir <path>      Output directory [default: $OUT_DIR]"
    echo "  --coverages <csv>     Coverage levels [default: $COVERAGES]"
    echo "  --threads <int>       Threads for bwa/samtools [default: $THREADS]"
    echo "  --samples <csv>       Specific samples to process (default: all)"
    echo "  --help                Show this help"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta-dir) FASTA_DIR="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --coverages) COVERAGES="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --samples) SAMPLE_LIST="$2"; shift 2 ;;
        --help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

echo "=== Read Simulation + GL Calling Pipeline ==="
echo "FASTA dir:  $FASTA_DIR"
echo "Output dir: $OUT_DIR"
echo "Coverages:  $COVERAGES"
echo "Threads:    $THREADS"
echo ""

# Check reference exists
REF="${FASTA_DIR}/reference.fa"
if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference not found: $REF"
    echo "Run 23_generate_sample_fastas.R first"
    exit 1
fi

# Create output directories
mkdir -p "${OUT_DIR}/bams" "${OUT_DIR}/gl"

# Index reference if needed
if [[ ! -f "${REF}.bwt" ]]; then
    echo "Indexing reference with bwa..."
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "Indexing reference with samtools..."
    samtools faidx "$REF"
fi

# Get sample list
if [[ -n "$SAMPLE_LIST" ]]; then
    IFS=',' read -ra SAMPLES <<< "$SAMPLE_LIST"
else
    SAMPLES=($(ls "${FASTA_DIR}"/sample*.fa 2>/dev/null | xargs -n1 basename | sed 's/.fa$//' | sort -V))
fi

echo "Processing ${#SAMPLES[@]} samples"

# Parse coverages
IFS=',' read -ra COV_ARRAY <<< "$COVERAGES"

# Process each sample at each coverage
TOTAL=$((${#SAMPLES[@]} * ${#COV_ARRAY[@]}))
CURRENT=0

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_FA="${FASTA_DIR}/${SAMPLE}.fa"

    if [[ ! -f "$SAMPLE_FA" ]]; then
        echo "WARNING: Sample FASTA not found: $SAMPLE_FA, skipping"
        continue
    fi

    for COV in "${COV_ARRAY[@]}"; do
        CURRENT=$((CURRENT + 1))
        echo ""
        echo "[$CURRENT/$TOTAL] ${SAMPLE} @ ${COV}x"

        # Output paths
        READS_PREFIX="${OUT_DIR}/${SAMPLE}_cov${COV}"
        BAM="${OUT_DIR}/bams/${SAMPLE}_cov${COV}.bam"
        GL_FILE="${OUT_DIR}/gl/${SAMPLE}_cov${COV}.gl.vcf.gz"

        # Skip if GL file already exists
        if [[ -f "$GL_FILE" ]]; then
            echo "  GL file exists, skipping"
            continue
        fi

        # Step 1: Simulate reads with ART
        # Simulate from both haplotypes at half coverage each
        HALF_COV=$(echo "$COV / 2" | bc -l)

        echo "  Simulating reads (${COV}x = ${HALF_COV}x per haplotype)..."

        # Extract haplotypes to temp files
        TMP_HAP1=$(mktemp).fa
        TMP_HAP2=$(mktemp).fa

        # Use seqkit to extract each haplotype
        seqkit grep -p "_hap1" "$SAMPLE_FA" > "$TMP_HAP1"
        seqkit grep -p "_hap2" "$SAMPLE_FA" > "$TMP_HAP2"

        # Simulate from each haplotype
        art_illumina -ss HS25 -i "$TMP_HAP1" -p -l 150 -f "$HALF_COV" -m 350 -s 50 \
            -o "${READS_PREFIX}_hap1_" -q 2>/dev/null || true

        art_illumina -ss HS25 -i "$TMP_HAP2" -p -l 150 -f "$HALF_COV" -m 350 -s 50 \
            -o "${READS_PREFIX}_hap2_" -q 2>/dev/null || true

        # Combine reads from both haplotypes
        cat "${READS_PREFIX}_hap1_1.fq" "${READS_PREFIX}_hap2_1.fq" > "${READS_PREFIX}_R1.fq"
        cat "${READS_PREFIX}_hap1_2.fq" "${READS_PREFIX}_hap2_2.fq" > "${READS_PREFIX}_R2.fq"

        # Clean up temp files
        rm -f "$TMP_HAP1" "$TMP_HAP2"
        rm -f "${READS_PREFIX}_hap1_"*.fq "${READS_PREFIX}_hap2_"*.fq
        rm -f "${READS_PREFIX}_hap1_"*.aln "${READS_PREFIX}_hap2_"*.aln 2>/dev/null || true

        # Step 2: Align with bwa
        echo "  Aligning reads..."
        bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" "$REF" \
            "${READS_PREFIX}_R1.fq" "${READS_PREFIX}_R2.fq" 2>/dev/null | \
            samtools sort -@ "$THREADS" -o "$BAM" -

        samtools index "$BAM"

        # Clean up FASTQ files
        rm -f "${READS_PREFIX}_R1.fq" "${READS_PREFIX}_R2.fq"

        # Step 3: Call GLs with bcftools mpileup
        echo "  Calling genotype likelihoods..."
        bcftools mpileup -f "$REF" -a FORMAT/AD,FORMAT/DP "$BAM" 2>/dev/null | \
            bcftools call -m -v -Oz -o "$GL_FILE" 2>/dev/null

        bcftools index "$GL_FILE"

        echo "  Done: $GL_FILE"
    done
done

echo ""
echo "=== Complete ==="
echo "BAMs:      ${OUT_DIR}/bams/"
echo "GL VCFs:   ${OUT_DIR}/gl/"
echo ""
echo "GL files contain PL (Phred-scaled likelihoods)."
echo "To convert PL to GL (log10): GL = -PL/10"
