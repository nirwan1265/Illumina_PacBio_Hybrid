#!/bin/bash
# =============================================================================
# BC2S3 Introgression Simulation Pipeline
# =============================================================================
#
# Complete pipeline to simulate BC2S3-like populations with:
#   1. Controlled introgressions (~12.5% donor genome)
#   2. Known truth files for benchmarking
#   3. Simulated Illumina reads at multiple coverages
#   4. Variant calls with genotype likelihoods (PL/GL)
#
# Usage: ./30_bc2s3_pipeline.sh [options]
#
# Example:
#   ./30_bc2s3_pipeline.sh --n-samples 20 --coverages 4,1 --threads 4
#
# =============================================================================

set -euo pipefail

# Ensure conda tools are in PATH
export PATH="/opt/anaconda3/envs/simitall/bin:$PATH"

# =============================================================================
# Default Parameters
# =============================================================================
N_SAMPLES=100
GENOME_SIZE=100000000      # 100 Mb
SNP_DENSITY=0.01           # 1 SNP per 100 bp = 1M SNPs
COVERAGES="4,1,0.5"
THREADS=4
OUT_PREFIX="05_summary/bc2s3"
SEED=42
SKIP_FASTA=false
SKIP_READS=false
ONLY_SAMPLES=""            # Process only specific samples (csv)

# =============================================================================
# Usage
# =============================================================================
usage() {
    cat << EOF
Usage: $0 [options]

BC2S3 Introgression Simulation Pipeline

Options:
  --n-samples <int>      Number of samples to simulate [default: $N_SAMPLES]
  --genome-size <int>    Genome size in bp [default: $GENOME_SIZE]
  --snp-density <float>  SNPs per bp [default: $SNP_DENSITY]
  --coverages <csv>      Coverage levels for read simulation [default: $COVERAGES]
  --threads <int>        Number of threads [default: $THREADS]
  --out-prefix <path>    Output prefix [default: $OUT_PREFIX]
  --seed <int>           Random seed [default: $SEED]
  --skip-fasta           Skip FASTA generation (use existing)
  --skip-reads           Skip read simulation (only generate VCF + FASTAs)
  --only-samples <csv>   Only process specific samples (e.g., "sample1,sample2")
  --help                 Show this help

Pipeline Steps:
  1. Generate VCF with controlled introgressions (~12.5% donor)
  2. Create diploid FASTAs for each sample
  3. Simulate Illumina reads at each coverage
  4. Call variants with bcftools (outputs PL/GL)

Output Structure:
  <out_prefix>_controlled.vcf.gz           - Multi-sample VCF with true genotypes
  <out_prefix>_controlled.introgressions.tsv - Truth file with introgression locations
  <out_prefix>_fastas/                     - Diploid FASTAs per sample
  <out_prefix>_reads/bams/                 - Aligned BAM files
  <out_prefix>_reads/vcf/                  - VCF files with GL/PL values

Examples:
  # Quick test (5 samples, 2 coverages)
  $0 --n-samples 5 --coverages 4,1 --threads 4

  # Full simulation (100 samples, 4 coverages)
  $0 --n-samples 100 --coverages 8,4,1,0.5 --threads 8

  # Resume from existing FASTAs
  $0 --skip-fasta --coverages 4,1

EOF
    exit 1
}

# =============================================================================
# Parse Arguments
# =============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --n-samples) N_SAMPLES="$2"; shift 2 ;;
        --genome-size) GENOME_SIZE="$2"; shift 2 ;;
        --snp-density) SNP_DENSITY="$2"; shift 2 ;;
        --coverages) COVERAGES="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --skip-fasta) SKIP_FASTA=true; shift ;;
        --skip-reads) SKIP_READS=true; shift ;;
        --only-samples) ONLY_SAMPLES="$2"; shift 2 ;;
        --help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# =============================================================================
# Derived paths
# =============================================================================
VCF_OUT="${OUT_PREFIX}_controlled.vcf.gz"
FASTA_DIR="${OUT_PREFIX}_fastas"
READS_DIR="${OUT_PREFIX}_reads"

echo "=============================================="
echo "  BC2S3 Introgression Simulation Pipeline"
echo "=============================================="
echo ""
echo "Parameters:"
echo "  Samples:     $N_SAMPLES"
echo "  Genome:      $GENOME_SIZE bp"
echo "  SNP density: $SNP_DENSITY"
echo "  Coverages:   $COVERAGES"
echo "  Threads:     $THREADS"
echo "  Output:      $OUT_PREFIX"
echo "  Seed:        $SEED"
echo ""

# =============================================================================
# Step 1: Generate VCF with Controlled Introgressions
# =============================================================================
if [[ -f "$VCF_OUT" ]] && [[ "$SKIP_FASTA" == "true" ]]; then
    echo "=== Step 1: VCF exists, skipping ==="
    echo "  $VCF_OUT"
else
    echo "=== Step 1: Generating VCF with Controlled Introgressions ==="
    echo ""

    Rscript scripts/mac/22_simulate_controlled_introgressions.R \
        --genome_size "$GENOME_SIZE" \
        --n_samples "$N_SAMPLES" \
        --snp_density "$SNP_DENSITY" \
        --out_prefix "${OUT_PREFIX}_controlled" \
        --seed "$SEED"

    echo ""
    echo "  VCF:    $VCF_OUT"
    echo "  Truth:  ${OUT_PREFIX}_controlled.introgressions.tsv"
    echo ""
fi

# =============================================================================
# Step 2: Generate Diploid FASTAs
# =============================================================================
if [[ "$SKIP_FASTA" == "true" ]]; then
    echo "=== Step 2: FASTA generation skipped ==="
else
    echo "=== Step 2: Generating Diploid FASTAs ==="
    echo ""

    mkdir -p "$FASTA_DIR"

    REF="${FASTA_DIR}/reference.fa"

    # Create reference if it doesn't exist
    if [[ ! -f "$REF" ]]; then
        echo "Creating reference FASTA..."
        # Extract reference alleles from VCF and create FASTA
        python3 << PYEOF
import gzip

# Read VCF and build reference sequence
ref_seq = ['N'] * $GENOME_SIZE

with gzip.open("$VCF_OUT", 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        pos = int(parts[1])
        ref = parts[3]
        if pos <= $GENOME_SIZE:
            ref_seq[pos-1] = ref[0]

# Fill gaps with random bases
import random
random.seed($SEED)
bases = 'ACGT'
for i in range($GENOME_SIZE):
    if ref_seq[i] == 'N':
        ref_seq[i] = random.choice(bases)

# Write FASTA
with open("$REF", 'w') as f:
    f.write(">chr1\\n")
    seq = ''.join(ref_seq)
    for i in range(0, len(seq), 80):
        f.write(seq[i:i+80] + "\\n")

print(f"Reference created: {len(seq)} bp")
PYEOF

        samtools faidx "$REF"
        bgzip -c "$VCF_OUT" > "${VCF_OUT}.tmp" 2>/dev/null || true
        bcftools index -f "$VCF_OUT"
    fi

    # Get sample list
    SAMPLES=($(bcftools query -l "$VCF_OUT" | head -n "$N_SAMPLES"))

    echo "Generating FASTAs for ${#SAMPLES[@]} samples..."
    echo ""

    CURRENT=0
    for S in "${SAMPLES[@]}"; do
        CURRENT=$((CURRENT + 1))
        SAMPLE_FA="${FASTA_DIR}/${S}.fa"

        if [[ -f "$SAMPLE_FA" ]]; then
            SIZE=$(stat -f%z "$SAMPLE_FA" 2>/dev/null || stat -c%s "$SAMPLE_FA" 2>/dev/null)
            echo "[$CURRENT/${#SAMPLES[@]}] $S - exists (${SIZE} bytes), skipping"
            continue
        fi

        echo "[$CURRENT/${#SAMPLES[@]}] $S"

        # Extract sample VCF
        TMP_VCF="/tmp/${S}.vcf.gz"
        bcftools view -s "$S" "$VCF_OUT" -Oz -o "$TMP_VCF" 2>/dev/null
        bcftools index -f "$TMP_VCF" 2>/dev/null

        # Generate haplotype 1
        bcftools consensus -f "$REF" -H 1 "$TMP_VCF" 2>/dev/null | \
            sed "s/^>.*/>chr1_${S}_hap1/" > "/tmp/${S}_h1.fa"

        # Generate haplotype 2
        bcftools consensus -f "$REF" -H 2 "$TMP_VCF" 2>/dev/null | \
            sed "s/^>.*/>chr1_${S}_hap2/" > "/tmp/${S}_h2.fa"

        # Combine
        cat "/tmp/${S}_h1.fa" "/tmp/${S}_h2.fa" > "$SAMPLE_FA"

        # Cleanup
        rm -f "$TMP_VCF" "${TMP_VCF}.csi" "/tmp/${S}_h1.fa" "/tmp/${S}_h2.fa"
    done

    echo ""
    echo "  FASTAs: $FASTA_DIR/"
    echo ""
fi

# =============================================================================
# Step 3: Simulate Reads and Call Variants
# =============================================================================
if [[ "$SKIP_READS" == "true" ]]; then
    echo "=== Step 3: Read simulation skipped ==="
else
    echo "=== Step 3: Simulating Reads and Calling Variants ==="
    echo ""

    mkdir -p "${READS_DIR}/bams" "${READS_DIR}/vcf"

    REF="${FASTA_DIR}/reference.fa"

    # Index reference if needed
    if [[ ! -f "${REF}.bwt" ]]; then
        echo "Indexing reference with bwa..."
        bwa index "$REF"
    fi

    if [[ ! -f "${REF}.fai" ]]; then
        samtools faidx "$REF"
    fi

    # Get samples
    if [[ -n "$ONLY_SAMPLES" ]]; then
        IFS=',' read -ra SAMPLE_ARR <<< "$ONLY_SAMPLES"
    else
        SAMPLE_ARR=($(ls "${FASTA_DIR}"/sample*.fa 2>/dev/null | xargs -n1 basename | sed 's/.fa$//' | sort -V | head -n "$N_SAMPLES"))
    fi

    IFS=',' read -ra COV_ARR <<< "$COVERAGES"

    TOTAL=$((${#SAMPLE_ARR[@]} * ${#COV_ARR[@]}))
    CURRENT=0

    echo "Processing ${#SAMPLE_ARR[@]} samples x ${#COV_ARR[@]} coverages = $TOTAL jobs"
    echo ""

    for S in "${SAMPLE_ARR[@]}"; do
        SAMPLE_FA="${FASTA_DIR}/${S}.fa"

        if [[ ! -f "$SAMPLE_FA" ]]; then
            echo "WARNING: $SAMPLE_FA not found, skipping"
            continue
        fi

        for COV in "${COV_ARR[@]}"; do
            CURRENT=$((CURRENT + 1))

            BAM="${READS_DIR}/bams/${S}_cov${COV}.bam"
            VCF="${READS_DIR}/vcf/${S}_cov${COV}.vcf.gz"

            if [[ -f "$VCF" ]]; then
                echo "[$CURRENT/$TOTAL] ${S} @ ${COV}x - exists, skipping"
                continue
            fi

            echo "[$CURRENT/$TOTAL] ${S} @ ${COV}x"

            # Half coverage per haplotype
            HALF_COV=$(echo "scale=2; $COV / 2" | bc)

            # Extract haplotypes
            seqkit grep -r -p "hap1" "$SAMPLE_FA" > "/tmp/${S}_h1.fa" 2>/dev/null
            seqkit grep -r -p "hap2" "$SAMPLE_FA" > "/tmp/${S}_h2.fa" 2>/dev/null

            # Simulate reads from each haplotype
            art_illumina -ss HS25 -i "/tmp/${S}_h1.fa" -p -l 150 -f "$HALF_COV" -m 350 -s 50 \
                -o "/tmp/${S}_h1_" -q 2>/dev/null || true
            art_illumina -ss HS25 -i "/tmp/${S}_h2.fa" -p -l 150 -f "$HALF_COV" -m 350 -s 50 \
                -o "/tmp/${S}_h2_" -q 2>/dev/null || true

            # Combine reads
            cat "/tmp/${S}_h1_1.fq" "/tmp/${S}_h2_1.fq" > "/tmp/${S}_R1.fq" 2>/dev/null
            cat "/tmp/${S}_h1_2.fq" "/tmp/${S}_h2_2.fq" > "/tmp/${S}_R2.fq" 2>/dev/null

            # Align
            bwa mem -t "$THREADS" -R "@RG\\tID:${S}\\tSM:${S}" "$REF" \
                "/tmp/${S}_R1.fq" "/tmp/${S}_R2.fq" 2>/dev/null | \
                samtools sort -@ "$THREADS" -o "$BAM" -

            samtools index "$BAM"

            # Call variants at known sites only (fast)
            # PL is output by default with -m, use -a GQ for genotype quality
            if [[ -f "$VCF_OUT" ]]; then
                bcftools mpileup -f "$REF" -T "$VCF_OUT" -a FORMAT/AD,FORMAT/DP "$BAM" 2>/dev/null | \
                    bcftools call -m -a GQ -Oz -o "$VCF"
            else
                bcftools mpileup -f "$REF" -a FORMAT/AD,FORMAT/DP "$BAM" 2>/dev/null | \
                    bcftools call -m -a GQ -Oz -o "$VCF"
            fi

            bcftools index -f "$VCF"

            # Cleanup temp files
            rm -f /tmp/${S}_h*.fa /tmp/${S}_h*_*.fq /tmp/${S}_h*_*.aln /tmp/${S}_R*.fq 2>/dev/null

            echo "  -> $VCF"
        done
    done
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "  Pipeline Complete"
echo "=============================================="
echo ""
echo "Outputs:"
echo "  VCF (truth):      $VCF_OUT"
echo "  Introgressions:   ${OUT_PREFIX}_controlled.introgressions.tsv"
echo "  Sample FASTAs:    ${FASTA_DIR}/"
echo "  BAMs:             ${READS_DIR}/bams/"
echo "  VCFs (with GL):   ${READS_DIR}/vcf/"
echo ""
echo "VCF files contain PL (Phred-scaled likelihoods)."
echo "To convert to GL (log10): GL = -PL/10"
echo ""
echo "Extract GLs:"
echo "  bcftools query -f '%CHROM\\t%POS\\t[%PL\\t]\\n' <file>.vcf.gz"
echo ""
