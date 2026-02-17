#!/bin/bash
# Simulate BC2S3 population with controlled introgression and generate reads
# at multiple coverage levels (10x, 8x, 4x, 1x, 0.5x)
#
# Usage: ./21_simulate_bc2s3_reads.sh [options]
#
# Requires: R, ART (art_illumina), pigz (all in simitall conda env)

set -euo pipefail

# Defaults
PANEL="inst/extdata/ref_files/aligned/b73_zdiplo_chr10_panel.fa"
OUT_PREFIX="05_summary/bc2s3_maize"
N_OFFSPRING=50
SEED=42
COVERAGES="10,8,4,1,0.5"
READ_TYPE="illumina"  # illumina, pacbio_clr, pacbio_hifi, or all

# Introgression settings
FIX_LOCUS=""  # e.g., "1000000:2000000" to fix a 1Mb region from Zdiplo
BACKGROUND_SELECTION="false"

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --panel <path>         Haplotype panel FASTA, default: $PANEL"
    echo "  --out-prefix <path>    Output prefix, default: $OUT_PREFIX"
    echo "  --n-offspring <int>    Number of offspring, default: $N_OFFSPRING"
    echo "  --coverages <csv>      Coverage levels (comma-sep), default: $COVERAGES"
    echo "  --read-type <type>     illumina, pacbio_clr, pacbio_hifi, or all"
    echo "  --fix-locus <s:e>      Fix introgression region (start:end bp)"
    echo "  --background-selection Enable background selection for minimal donor"
    echo "  --seed <int>           Random seed, default: $SEED"
    echo "  --help                 Show this help"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --panel) PANEL="$2"; shift 2 ;;
        --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
        --n-offspring) N_OFFSPRING="$2"; shift 2 ;;
        --coverages) COVERAGES="$2"; shift 2 ;;
        --read-type) READ_TYPE="$2"; shift 2 ;;
        --fix-locus) FIX_LOCUS="$2"; shift 2 ;;
        --background-selection) BACKGROUND_SELECTION="true"; shift ;;
        --seed) SEED="$2"; shift 2 ;;
        --help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

echo "=== BC2S3 Breeding + Read Simulation Pipeline ==="
echo "Panel:       $PANEL"
echo "Output:      $OUT_PREFIX"
echo "Offspring:   $N_OFFSPRING"
echo "Coverages:   $COVERAGES"
echo "Read type:   $READ_TYPE"
echo "Fix locus:   ${FIX_LOCUS:-none}"
echo "BG select:   $BACKGROUND_SELECTION"
echo "Seed:        $SEED"
echo ""

# Check panel exists
if [[ ! -f "$PANEL" ]]; then
    echo "ERROR: Panel not found: $PANEL"
    echo "Run 20_prepare_aligned_haplotypes.sh first"
    exit 1
fi

# Create output directories
BREED_DIR="$(dirname "$OUT_PREFIX")"
READS_DIR="02_reads/bc2s3"
mkdir -p "$BREED_DIR" "$READS_DIR"

# Step 1: Run breeding simulation
echo "[1/3] Running BC2S3 breeding simulation..."

BREED_ARGS=(
    --haplotype_fa "$PANEL"
    --out_prefix "$OUT_PREFIX"
    --parents "B73,Zdiplo"
    --sequence "F1,BC:P1:1,BC:P1:1,SELF:1,SELF:1,SELF:1"
    --n_offspring "$N_OFFSPRING"
    --vcf_out "${OUT_PREFIX}.vcf"
    --graph_out "${OUT_PREFIX}.mmd"
    --qc_out "${OUT_PREFIX}.qc.tsv"
    --seed "$SEED"
)

# Add introgression options
if [[ -n "$FIX_LOCUS" ]]; then
    BREED_ARGS+=(--fix_locus "$FIX_LOCUS" --fix_allele donor)
fi

if [[ "$BACKGROUND_SELECTION" == "true" ]]; then
    BREED_ARGS+=(--background_selection --selection_pool 100 --marker_step 500)
fi

Rscript scripts/mac/11_simulate_breeding.R "${BREED_ARGS[@]}"

echo "  Breeding complete: ${OUT_PREFIX}.fa"

# Step 2: Split bred FASTA into individual samples
echo ""
echo "[2/3] Splitting bred FASTA into individual samples..."

SAMPLES_DIR="${OUT_PREFIX}_samples"
mkdir -p "$SAMPLES_DIR"

# Use seqkit to split by sample (each sample has 2 haplotypes)
# First, get unique sample IDs
seqkit seq -n "${OUT_PREFIX}.fa" | sed 's/_hap[12]$//' | sort -u > "${SAMPLES_DIR}/sample_ids.txt"
N_SAMPLES=$(wc -l < "${SAMPLES_DIR}/sample_ids.txt" | tr -d ' ')
echo "  Found $N_SAMPLES samples"

# Create per-sample FASTA files (diploid: both haplotypes combined)
while read -r SAMPLE; do
    SAMPLE_FA="${SAMPLES_DIR}/${SAMPLE}.fa"
    seqkit grep -r -p "^${SAMPLE}_hap" "${OUT_PREFIX}.fa" > "$SAMPLE_FA"
done < "${SAMPLES_DIR}/sample_ids.txt"

echo "  Split into ${SAMPLES_DIR}/"

# Step 3: Simulate reads at each coverage
echo ""
echo "[3/3] Simulating reads at coverages: $COVERAGES"

IFS=',' read -ra COV_ARRAY <<< "$COVERAGES"

simulate_illumina() {
    local REF="$1"
    local OUTPFX="$2"
    local COV="$3"

    art_illumina -ss HS25 -i "$REF" -p -l 150 -f "$COV" -m 350 -s 50 -o "$OUTPFX" -q 2>/dev/null

    # Compress
    [[ -f "${OUTPFX}1.fq" ]] && pigz -f "${OUTPFX}1.fq"
    [[ -f "${OUTPFX}2.fq" ]] && pigz -f "${OUTPFX}2.fq"
}

simulate_pacbio_clr() {
    local REF="$1"
    local OUTDIR="$2"
    local COV="$3"

    mkdir -p "$OUTDIR"
    pbsim --strategy wgs --method qshmm --qshmm QSHMM-RSII.model \
        --depth "$COV" --genome "$REF" --prefix "${OUTDIR}/clr" 2>/dev/null || true

    # Convert sam to fastq if needed
    if [[ -f "${OUTDIR}/clr_0001.sam" ]]; then
        samtools fastq "${OUTDIR}/clr_0001.sam" 2>/dev/null | pigz > "${OUTDIR}/clr.fq.gz"
        rm -f "${OUTDIR}/clr_0001.sam" "${OUTDIR}/clr_0001.ref" "${OUTDIR}/clr_0001.maf" 2>/dev/null || true
    fi
}

simulate_pacbio_hifi() {
    local REF="$1"
    local OUTDIR="$2"
    local COV="$3"

    mkdir -p "$OUTDIR"
    pbsim --strategy wgs --method qshmm --qshmm QSHMM-SEQUEL.model \
        --depth "$COV" --accuracy-mean 0.999 --genome "$REF" --prefix "${OUTDIR}/hifi" 2>/dev/null || true

    if [[ -f "${OUTDIR}/hifi_0001.sam" ]]; then
        samtools fastq "${OUTDIR}/hifi_0001.sam" 2>/dev/null | pigz > "${OUTDIR}/hifi.fq.gz"
        rm -f "${OUTDIR}/hifi_0001.sam" "${OUTDIR}/hifi_0001.ref" "${OUTDIR}/hifi_0001.maf" 2>/dev/null || true
    fi
}

# Progress tracking
TOTAL_JOBS=$((N_SAMPLES * ${#COV_ARRAY[@]}))
CURRENT=0

while read -r SAMPLE; do
    SAMPLE_FA="${SAMPLES_DIR}/${SAMPLE}.fa"

    for COV in "${COV_ARRAY[@]}"; do
        CURRENT=$((CURRENT + 1))
        echo "  [$CURRENT/$TOTAL_JOBS] ${SAMPLE} @ ${COV}x"

        # Illumina
        if [[ "$READ_TYPE" == "illumina" || "$READ_TYPE" == "all" ]]; then
            ILL_OUT="${READS_DIR}/${SAMPLE}_ill${COV}x_"
            if [[ ! -f "${ILL_OUT}1.fq.gz" ]]; then
                simulate_illumina "$SAMPLE_FA" "$ILL_OUT" "$COV"
            fi
        fi

        # PacBio CLR
        if [[ "$READ_TYPE" == "pacbio_clr" || "$READ_TYPE" == "all" ]]; then
            CLR_OUT="${READS_DIR}/${SAMPLE}_clr${COV}x"
            if [[ ! -f "${CLR_OUT}/clr.fq.gz" ]]; then
                simulate_pacbio_clr "$SAMPLE_FA" "$CLR_OUT" "$COV"
            fi
        fi

        # PacBio HiFi
        if [[ "$READ_TYPE" == "pacbio_hifi" || "$READ_TYPE" == "all" ]]; then
            HIFI_OUT="${READS_DIR}/${SAMPLE}_hifi${COV}x"
            if [[ ! -f "${HIFI_OUT}/hifi.fq.gz" ]]; then
                simulate_pacbio_hifi "$SAMPLE_FA" "$HIFI_OUT" "$COV"
            fi
        fi
    done
done < "${SAMPLES_DIR}/sample_ids.txt"

echo ""
echo "=== Complete ==="
echo ""
echo "Outputs:"
echo "  Bred population:  ${OUT_PREFIX}.fa"
echo "  Metadata:         ${OUT_PREFIX}.meta.tsv"
echo "  VCF:              ${OUT_PREFIX}.vcf"
echo "  QC report:        ${OUT_PREFIX}.qc.tsv"
echo "  Crossing graph:   ${OUT_PREFIX}.mmd"
echo "  Per-sample FASTA: ${SAMPLES_DIR}/"
echo "  Simulated reads:  ${READS_DIR}/"
echo ""
echo "Read files pattern:"
echo "  Illumina: ${READS_DIR}/sampleN_illXx_1.fq.gz, ${READS_DIR}/sampleN_illXx_2.fq.gz"
if [[ "$READ_TYPE" == "pacbio_clr" || "$READ_TYPE" == "all" ]]; then
    echo "  PacBio CLR: ${READS_DIR}/sampleN_clrXx/clr.fq.gz"
fi
if [[ "$READ_TYPE" == "pacbio_hifi" || "$READ_TYPE" == "all" ]]; then
    echo "  PacBio HiFi: ${READS_DIR}/sampleN_hifiXx/hifi.fq.gz"
fi
