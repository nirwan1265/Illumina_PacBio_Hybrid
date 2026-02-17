#!/bin/bash
# Simulate reads at multiple coverages and call GLs with bcftools
#
# Usage: ./25_simulate_reads_call_gl.sh [options]
#
# Requires: art_illumina, bwa/minimap2, samtools, bcftools

set -euo pipefail

# Ensure conda tools are in PATH
export PATH="/opt/anaconda3/envs/simitall/bin:$PATH"

FASTA_DIR="05_summary/bc2s3_fastas"
OUT_DIR="05_summary/bc2s3_reads"
COVERAGES="10,8,4,1,0.5"
THREADS=4
SAMPLES=""
N_SAMPLES=""   # Limit to first N samples

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --fasta-dir <path>    Sample FASTAs directory [default: $FASTA_DIR]"
    echo "  --out-dir <path>      Output directory [default: $OUT_DIR]"
    echo "  --coverages <csv>     Coverage levels [default: $COVERAGES]"
    echo "  --threads <int>       Threads [default: $THREADS]"
    echo "  --n-samples <int>     Limit to first N samples [default: all]"
    echo "  --samples <csv>       Specific samples to process (e.g., sample1,sample5)"
    echo "  --help                Show help"
    echo ""
    echo "Examples:"
    echo "  $0 --n-samples 10 --coverages 4,1 --threads 4"
    echo "  $0 --samples sample1,sample2,sample3 --coverages 8"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta-dir) FASTA_DIR="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --coverages) COVERAGES="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --n-samples) N_SAMPLES="$2"; shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        --help) usage ;;
        *) echo "Unknown: $1"; usage ;;
    esac
done

REF="${FASTA_DIR}/reference.fa"
# Use the original VCF to restrict variant calling to known positions (MUCH faster)
SITES_VCF="05_summary/bc2s3_controlled.vcf.gz"

echo "=== Read Simulation + GL Calling ==="
echo "Reference: $REF"
echo "Sites:     $SITES_VCF"
echo "Output:    $OUT_DIR"
echo "Coverages: $COVERAGES"
echo "Threads:   $THREADS"
echo ""

# Check reference
if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference not found: $REF"
    exit 1
fi

# Check sites VCF
if [[ ! -f "$SITES_VCF" ]]; then
    echo "WARNING: Sites VCF not found: $SITES_VCF"
    echo "Will call at ALL positions (very slow for large genomes)"
    SITES_VCF=""
fi

mkdir -p "${OUT_DIR}/bams" "${OUT_DIR}/vcf"

# Index reference if needed
if [[ ! -f "${REF}.bwt" ]]; then
    echo "Indexing reference with bwa..."
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# Get samples
if [[ -n "$SAMPLES" ]]; then
    IFS=',' read -ra SAMPLE_ARR <<< "$SAMPLES"
elif [[ -n "$N_SAMPLES" ]]; then
    SAMPLE_ARR=($(ls "${FASTA_DIR}"/sample*.fa 2>/dev/null | xargs -n1 basename | sed 's/.fa$//' | sort -V | head -n "$N_SAMPLES"))
else
    SAMPLE_ARR=($(ls "${FASTA_DIR}"/sample*.fa 2>/dev/null | xargs -n1 basename | sed 's/.fa$//' | sort -V))
fi

echo "Samples:   ${#SAMPLE_ARR[@]}"

IFS=',' read -ra COV_ARR <<< "$COVERAGES"

TOTAL=$((${#SAMPLE_ARR[@]} * ${#COV_ARR[@]}))
CURRENT=0

for S in "${SAMPLE_ARR[@]}"; do
    SAMPLE_FA="${FASTA_DIR}/${S}.fa"

    if [[ ! -f "$SAMPLE_FA" ]]; then
        echo "WARNING: $SAMPLE_FA not found, skipping"
        continue
    fi

    for COV in "${COV_ARR[@]}"; do
        CURRENT=$((CURRENT + 1))

        BAM="${OUT_DIR}/bams/${S}_cov${COV}.bam"
        VCF="${OUT_DIR}/vcf/${S}_cov${COV}.vcf.gz"

        if [[ -f "$VCF" ]]; then
            echo "[$CURRENT/$TOTAL] ${S} @ ${COV}x - exists, skipping"
            continue
        fi

        echo "[$CURRENT/$TOTAL] ${S} @ ${COV}x"

        # Half coverage per haplotype
        HALF_COV=$(echo "scale=2; $COV / 2" | bc)

        # Extract haplotypes
        echo "  Extracting haplotypes..."
        seqkit grep -r -p "hap1" "$SAMPLE_FA" > "/tmp/${S}_h1.fa" 2>/dev/null
        seqkit grep -r -p "hap2" "$SAMPLE_FA" > "/tmp/${S}_h2.fa" 2>/dev/null

        # Simulate reads from each haplotype
        echo "  Simulating reads (${HALF_COV}x per haplotype)..."
        art_illumina -ss HS25 -i "/tmp/${S}_h1.fa" -p -l 150 -f "$HALF_COV" -m 350 -s 50 \
            -o "/tmp/${S}_h1_" -q 2>&1 | grep -v "^$" | head -5 || true
        art_illumina -ss HS25 -i "/tmp/${S}_h2.fa" -p -l 150 -f "$HALF_COV" -m 350 -s 50 \
            -o "/tmp/${S}_h2_" -q 2>&1 | grep -v "^$" | head -5 || true

        # Combine reads
        cat "/tmp/${S}_h1_1.fq" "/tmp/${S}_h2_1.fq" > "/tmp/${S}_R1.fq" 2>/dev/null
        cat "/tmp/${S}_h1_2.fq" "/tmp/${S}_h2_2.fq" > "/tmp/${S}_R2.fq" 2>/dev/null

        # Align
        echo "  Aligning reads..."
        bwa mem -t "$THREADS" -R "@RG\\tID:${S}\\tSM:${S}" "$REF" \
            "/tmp/${S}_R1.fq" "/tmp/${S}_R2.fq" 2>/dev/null | \
            samtools sort -@ "$THREADS" -o "$BAM" -

        samtools index "$BAM"

        # Call with bcftools (outputs PL which are Phred-scaled likelihoods)
        # Use -T to restrict to known SNP positions (MUCH faster than calling everywhere)
        # PL is output by default with -m caller
        echo "  Calling variants..."
        if [[ -n "$SITES_VCF" ]]; then
            # Fast mode: only call at known sites
            bcftools mpileup -f "$REF" -T "$SITES_VCF" -a FORMAT/AD,FORMAT/DP "$BAM" 2>/dev/null | \
                bcftools call -m -a GQ -Oz -o "$VCF"
        else
            # Slow mode: call everywhere
            bcftools mpileup -f "$REF" -a FORMAT/AD,FORMAT/DP "$BAM" 2>/dev/null | \
                bcftools call -m -a GQ -Oz -o "$VCF"
        fi

        bcftools index -f "$VCF"

        # Cleanup temp files
        rm -f /tmp/${S}_h*.fa /tmp/${S}_h*_*.fq /tmp/${S}_h*_*.aln /tmp/${S}_R*.fq 2>/dev/null

        echo "  -> $VCF"
    done
done

echo ""
echo "=== Complete ==="
echo "BAMs: ${OUT_DIR}/bams/"
echo "VCFs: ${OUT_DIR}/vcf/"
echo ""
echo "VCF files contain PL (Phred-scaled likelihoods)."
echo "To convert PL to GL (log10): GL = -PL/10"
echo ""
echo "Example to extract GLs:"
echo "  bcftools query -f '%CHROM\\t%POS\\t[%PL\\t]\\n' file.vcf.gz"
