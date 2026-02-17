#!/bin/bash
# Fix reference FASTA to match VCF REF alleles

set -euo pipefail

VCF_GZ="05_summary/bc2s3_controlled.vcf.gz"
REF="05_summary/bc2s3_fastas/reference.fa"
REF_FIXED="05_summary/bc2s3_fastas/reference_fixed.fa"

echo "Fixing reference to match VCF REF alleles..."

# Load reference into memory
REF_SEQ=$(grep -v "^>" "$REF" | tr -d '\n')
echo "Reference length: ${#REF_SEQ}"

# Create a temp file with variant positions and ref alleles
bcftools query -f '%POS\t%REF\n' "$VCF_GZ" > /tmp/snp_refs.tsv

# Apply REF alleles to reference using Python (faster than bash)
python3 << 'PYEOF'
import sys

# Read reference
with open("05_summary/bc2s3_fastas/reference.fa") as f:
    lines = f.readlines()
    header = lines[0]
    seq = ''.join(line.strip() for line in lines[1:])

seq_list = list(seq)
print(f"Original sequence length: {len(seq_list)}")

# Read and apply SNP ref alleles
count = 0
with open("/tmp/snp_refs.tsv") as f:
    for line in f:
        pos, ref = line.strip().split('\t')
        pos = int(pos) - 1  # 0-based
        if pos < len(seq_list):
            seq_list[pos] = ref
            count += 1

print(f"Applied {count} REF alleles")

# Write fixed reference
with open("05_summary/bc2s3_fastas/reference_fixed.fa", 'w') as f:
    f.write(">chr1\n")
    seq = ''.join(seq_list)
    for i in range(0, len(seq), 80):
        f.write(seq[i:i+80] + '\n')

print("Done!")
PYEOF

# Replace old reference
mv "$REF_FIXED" "$REF"
samtools faidx "$REF"

echo "Reference fixed!"
