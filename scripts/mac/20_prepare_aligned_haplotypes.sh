#!/bin/bash
# Prepare aligned haplotype panel from two reference genomes
# Uses minimap2 for whole-genome alignment and extracts syntenic regions
#
# Usage: ./20_prepare_aligned_haplotypes.sh [options]
#
# Requires: minimap2, samtools, seqkit (all in simitall conda env)

set -euo pipefail

# Defaults
B73_GZ="ref/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
ZDIPLO_GZ="ref/Zd-Momo-REFERENCE-PanAnd-1.0.fa.gz"
CHR="chr10"
OUTDIR="inst/extdata/ref_files/aligned"
MIN_BLOCK=100000  # Minimum syntenic block size (100kb)
SEED=42

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --b73 <path>       B73 reference FASTA (gzipped), default: $B73_GZ"
    echo "  --zdiplo <path>    Zdiploperennis reference FASTA (gzipped), default: $ZDIPLO_GZ"
    echo "  --chr <name>       Chromosome to extract, default: $CHR"
    echo "  --outdir <path>    Output directory, default: $OUTDIR"
    echo "  --min-block <int>  Minimum syntenic block size (bp), default: $MIN_BLOCK"
    echo "  --help             Show this help"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --b73) B73_GZ="$2"; shift 2 ;;
        --zdiplo) ZDIPLO_GZ="$2"; shift 2 ;;
        --chr) CHR="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --min-block) MIN_BLOCK="$2"; shift 2 ;;
        --help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

echo "=== Aligned Haplotype Panel Preparation ==="
echo "B73:        $B73_GZ"
echo "Zdiplo:     $ZDIPLO_GZ"
echo "Chromosome: $CHR"
echo "Output:     $OUTDIR"
echo "Min block:  $MIN_BLOCK bp"
echo ""

# Create output directory
mkdir -p "$OUTDIR"

# Step 1: Extract target chromosome from both genomes
echo "[1/5] Extracting $CHR from both genomes..."

B73_CHR="$OUTDIR/b73_${CHR}.fa"
ZDIPLO_CHR="$OUTDIR/zdiplo_${CHR}.fa"

# Remove old files if they exist but are empty
[[ -f "$B73_CHR" && ! -s "$B73_CHR" ]] && rm "$B73_CHR"
[[ -f "$ZDIPLO_CHR" && ! -s "$ZDIPLO_CHR" ]] && rm "$ZDIPLO_CHR"

if [[ ! -f "$B73_CHR" ]]; then
    seqkit grep -p "${CHR}" "$B73_GZ" | seqkit replace -p '.+' -r 'B73' > "$B73_CHR"
    echo "  Extracted B73 $CHR"
else
    echo "  B73 $CHR already exists, skipping"
fi

if [[ ! -f "$ZDIPLO_CHR" ]]; then
    seqkit grep -p "${CHR}" "$ZDIPLO_GZ" | seqkit replace -p '.+' -r 'Zdiplo' > "$ZDIPLO_CHR"
    echo "  Extracted Zdiplo $CHR"
else
    echo "  Zdiplo $CHR already exists, skipping"
fi

# Get chromosome lengths
B73_LEN=$(seqkit stats -T "$B73_CHR" | tail -1 | cut -f5)
ZDIPLO_LEN=$(seqkit stats -T "$ZDIPLO_CHR" | tail -1 | cut -f5)
echo "  B73 $CHR length:    $B73_LEN bp"
echo "  Zdiplo $CHR length: $ZDIPLO_LEN bp"

# Step 2: Align Zdiplo to B73 using minimap2
echo ""
echo "[2/5] Aligning Zdiplo to B73 with minimap2..."

PAF="$OUTDIR/zdiplo_to_b73_${CHR}.paf"
if [[ ! -f "$PAF" ]]; then
    minimap2 -x asm20 -c --cs "$B73_CHR" "$ZDIPLO_CHR" > "$PAF"
    echo "  Alignment complete"
else
    echo "  Alignment file exists, skipping"
fi

# Step 3: Extract largest syntenic block
echo ""
echo "[3/5] Finding syntenic blocks (min ${MIN_BLOCK} bp)..."

# Parse PAF to find co-linear blocks
# PAF columns: qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, matches, alnlen, mapq
BLOCKS_TSV="$OUTDIR/syntenic_blocks_${CHR}.tsv"

# Extract blocks (header separate from data to avoid sorting issues)
echo -e "b73_start\tb73_end\tzdiplo_start\tzdiplo_end\tstrand\tlength\tmatches" > "$BLOCKS_TSV"

awk -v min="$MIN_BLOCK" '
$5 == "+" && ($4 - $3) >= min {
    len = $4 - $3
    print $8, $9, $3, $4, $5, len, $10
}
' OFS="\t" "$PAF" | sort -k6,6nr >> "$BLOCKS_TSV"

N_BLOCKS=$(tail -n +2 "$BLOCKS_TSV" | wc -l | tr -d ' ')
echo "  Found $N_BLOCKS syntenic blocks >= ${MIN_BLOCK} bp"

if [[ $N_BLOCKS -eq 0 ]]; then
    echo "ERROR: No syntenic blocks found. Try reducing --min-block"
    exit 1
fi

# Get the largest block (skip header, take first data line)
BEST_BLOCK=$(tail -n +2 "$BLOCKS_TSV" | head -1)
B73_START=$(echo "$BEST_BLOCK" | cut -f1)
B73_END=$(echo "$BEST_BLOCK" | cut -f2)
ZD_START=$(echo "$BEST_BLOCK" | cut -f3)
ZD_END=$(echo "$BEST_BLOCK" | cut -f4)
BLOCK_LEN=$(echo "$BEST_BLOCK" | cut -f6)

echo "  Largest block: B73 ${B73_START}-${B73_END}, Zdiplo ${ZD_START}-${ZD_END} (${BLOCK_LEN} bp)"

# Step 4: Extract aligned regions
echo ""
echo "[4/5] Extracting aligned regions..."

B73_REGION="$OUTDIR/b73_${CHR}_aligned.fa"
ZDIPLO_REGION="$OUTDIR/zdiplo_${CHR}_aligned.fa"

# Extract regions (1-based for seqkit)
B73_START_1=$((B73_START + 1))
ZD_START_1=$((ZD_START + 1))

seqkit subseq -r ${B73_START_1}:${B73_END} "$B73_CHR" | seqkit replace -p '.+' -r 'B73' > "$B73_REGION"
seqkit subseq -r ${ZD_START_1}:${ZD_END} "$ZDIPLO_CHR" | seqkit replace -p '.+' -r 'Zdiplo' > "$ZDIPLO_REGION"

# Get actual extracted lengths
B73_EXTR_LEN=$(seqkit stats -T "$B73_REGION" | tail -1 | cut -f5)
ZD_EXTR_LEN=$(seqkit stats -T "$ZDIPLO_REGION" | tail -1 | cut -f5)

echo "  B73 extracted:    $B73_EXTR_LEN bp"
echo "  Zdiplo extracted: $ZD_EXTR_LEN bp"

# Step 5: Create final haplotype panel (same length required)
echo ""
echo "[5/5] Creating haplotype panel..."

PANEL="$OUTDIR/b73_zdiplo_${CHR}_panel.fa"

# Use the shorter length for both
if [[ $B73_EXTR_LEN -le $ZD_EXTR_LEN ]]; then
    FINAL_LEN=$B73_EXTR_LEN
else
    FINAL_LEN=$ZD_EXTR_LEN
fi

echo "  Trimming both to $FINAL_LEN bp for equal-length panel"

# Trim to same length and combine
{
    seqkit subseq -r 1:${FINAL_LEN} "$B73_REGION"
    seqkit subseq -r 1:${FINAL_LEN} "$ZDIPLO_REGION"
} > "$PANEL"

# Verify panel
echo ""
echo "=== Final Panel ==="
seqkit stats -T "$PANEL" | column -t

echo ""
echo "Panel created: $PANEL"
echo ""
echo "Next steps:"
echo "  1. Run breeding simulation:"
echo "     Rscript scripts/mac/11_simulate_breeding.R \\"
echo "       --haplotype_fa $PANEL \\"
echo "       --out_prefix 05_summary/bc2s3_maize \\"
echo "       --parents B73,Zdiplo \\"
echo "       --sequence 'F1,BC:P1:1,BC:P1:1,SELF:1,SELF:1,SELF:1' \\"
echo "       --n_offspring 50 \\"
echo "       --vcf_out 05_summary/bc2s3_maize.vcf \\"
echo "       --seed 42"
echo ""
echo "  2. Simulate reads at multiple coverages (see 21_simulate_bc2s3_reads.sh)"
