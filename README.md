# Illumina_PacBio_Hybrid

## Introduction

This project provides a comprehensive toolkit for **genomic data simulation**, encompassing everything from bacterial genome benchmarking to complex plant breeding schemes. Whether you need to simulate Illumina short reads, PacBio long reads, cross-species introgression, or complete breeding populations with genotype likelihoods, this toolkit has you covered.

### What You Can Simulate

| Category | Capabilities |
|----------|-------------|
| **Genome Assembly** | Hybrid assembly benchmarking with Illumina + PacBio (CLR/HiFi) |
| **Read Simulation** | Illumina PE (ART), PacBio CLR/HiFi (PBSIM), coverage grids |
| **GWAS Cohorts** | Population structure, LD blocks, causal variants, phenotypes |
| **Breeding Populations** | F1, BC, RIL, NIL, DH, MAGIC, NAM schemes |
| **Cross-Species Breeding** | Align divergent haplotypes, simulate controlled introgressions |
| **Genotype Likelihoods** | Full pipeline: FASTA → reads → BAM → VCF with PL/GL values |

### Key Features

- **Controlled Introgression Simulation**: Generate BC2S3-like populations with known donor segments (~12.5% introgression)
- **Cross-Species Alignment**: Align chromosomes from different species (e.g., B73 maize × teosinte) for interspecific breeding simulation
- **Complete GL Pipeline**: From simulated diploid FASTAs through Illumina read simulation to bcftools variant calling
- **Flexible Breeding Schemes**: F1, selfing, backcrossing, doubled haploids, MAGIC, NAM with Mermaid graph outputs
- **Multi-Coverage Grids**: Simulate reads at multiple coverage levels (10x, 8x, 4x, 1x, 0.5x) for benchmarking

---

## Quick Start

### Installation

Create a conda environment with all required tools:

```bash
conda create -n simitall -y -c conda-forge -c bioconda \
  r-base r-jsonlite r-reticulate r-optparse r-data.table \
  art pbsim unicycler quast samtools bcftools bwa seqkit pigz simupop minimap2
```

Activate:
```bash
conda activate simitall
```

### Workspace Note (Mac)

Use a path **without spaces** for reliable SPAdes/QUAST runs:
```
/Users/yourname/Documents/Github/Illumina_PacBio_Hybrid
```

---

## Simulation Pipelines

### 1. Hybrid Assembly Benchmarking

Benchmark hybrid assembly performance across Illumina short reads and PacBio long reads (CLR and HiFi).

**Pipeline:**
```
Reference → Repeat Stress → Read Simulation → Hybrid Assembly → QUAST Evaluation
```

**Scripts (run in order):**
```bash
# 1. Create repeat-stress reference
Rscript scripts/mac/01_make_tandem_repeats.R \
  --in_fa 00_ref/ecoli.fa \
  --out_fa 01_simref/ecoli_repMed.fa \
  --n_events 15 --seg_len 1000 --copies 5

# 2. Generate annotations (optional)
Rscript scripts/mac/01b_generate_annotations.R \
  --genome_fa 01_simref/ecoli_repMed.fa \
  --out_gff3 01_simref/ecoli_repMed.gff3

# 3. Simulate reads (Illumina + PacBio grids)
Rscript scripts/mac/04_run_grid_both_pacbio.R

# 4. Run hybrid assembly
Rscript scripts/mac/05_run_unicycler_grid_both.R

# 5. Evaluate with QUAST
Rscript scripts/mac/06_run_quast_both.R

# 6. Summarize results
Rscript scripts/mac/07_summarize_quast.R \
  04_eval/ecoli_repMed/CLR \
  05_summary/ecoli_repMed.CLR.quast_summary.csv
```

---

### 2. GWAS Cohort Simulation

Simulate a GWAS-style cohort with population structure, LD blocks, and phenotypes.

**Outputs:**
- VCF with genotypes
- Genotype matrix (`.geno.tsv`)
- Phenotype table (`.pheno.tsv`)
- Recombination map (optional)

**Example (diploid, quantitative trait):**
```bash
Rscript scripts/mac/10_simulate_gwas_cohort.R \
  --genome_fa inst/extdata/ref_files/human_grch38_chr1.fa \
  --out_prefix 05_summary/gwas_chr1 \
  --n_samples 200 \
  --ploidy 2 \
  --snp_rate 0.001 \
  --indel_rate 0.0001 \
  --phenotype quantitative \
  --n_causal 50 \
  --effect_sd 0.6 \
  --n_pops 3 \
  --fst 0.05 \
  --ld_block_size 50000
```

**Options:**
| Parameter | Description |
|-----------|-------------|
| `--n_pops`, `--pop_sizes`, `--fst` | Population structure |
| `--ld_block_size`, `--ld_haplotypes` | LD blocks and strength |
| `--recomb_map_in/out` | Recombination map support |
| `--phenotype` | `quantitative` or `binary` |
| `--n_causal`, `--effect_sd` | Causal variants and effect sizes |

---

### 3. Breeding Population Simulation

Simulate F1, selfing, backcross, and advanced breeding schemes from haplotype panels.

**Supported Schemes:**
- **F1**: First filial generation from two parents
- **SELF**: Self-fertilization (SSD or SIB mating)
- **BC**: Backcross to P1 or P2
- **RIL**: Recombinant Inbred Lines
- **NIL**: Near-Isogenic Lines
- **DH**: Doubled Haploids
- **MAGIC**: Multi-parent Advanced Generation Inter-Cross
- **NAM**: Nested Association Mapping

**Example (F1 → Self 2 generations → Backcross):**
```bash
Rscript scripts/mac/11_simulate_breeding.R \
  --haplotype_fa inst/extdata/ref_files/random_panel.fa \
  --out_prefix 05_summary/breeding_demo \
  --parents hap1,hap2 \
  --sequence F1,SELF:2,BC:P1:1 \
  --n_offspring 100 \
  --seed 7
```

**MAGIC Example:**
```bash
Rscript scripts/mac/11_simulate_breeding.R \
  --haplotype_fa inst/extdata/ref_files/random_panel.fa \
  --out_prefix 05_summary/magic_demo \
  --scheme MAGIC \
  --n_founders 6 \
  --n_offspring 200
```

**NAM Example:**
```bash
Rscript scripts/mac/11_simulate_breeding.R \
  --haplotype_fa inst/extdata/ref_files/random_panel.fa \
  --out_prefix 05_summary/nam_demo \
  --scheme NAM \
  --founders hap1,hap2,hap3,hap4,hap5 \
  --n_offspring 200
```

**Outputs:**
- `<prefix>.fa` - Haplotype sequences (`sampleX_hap1`, `sampleX_hap2`)
- `<prefix>.meta.tsv` - Sample metadata (generation, scheme)
- `<prefix>.vcf` - VCF with genotypes
- `<prefix>.mmd` - Mermaid crossing diagram (optional)

---

### 4. Cross-Species Introgression Simulation

The most advanced feature: simulate BC2S3-like populations with controlled introgressions from a divergent donor species.

**Use Case:** Simulate introgression breeding from wild relatives (e.g., B73 maize × teosinte)

#### Step 1: Align Divergent Haplotypes

Align chromosomes from two species to find syntenic regions:

```bash
./scripts/mac/20_prepare_aligned_haplotypes.sh \
  --ref1 ref/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
  --ref2 ref/Zd-Momo-REFERENCE-PanAnd-1.0.fa.gz \
  --chr1 chr10 \
  --chr2 chr10 \
  --out-dir 05_summary/aligned_haps \
  --target-len 100000000
```

**Outputs:**
- Syntenic blocks (PAF format)
- Aligned haplotype FASTAs of equal length
- Coordinate mapping files

#### Step 2: Simulate Controlled Introgressions

Generate a BC2S3-like population with known introgression locations:

```bash
Rscript scripts/mac/22_simulate_controlled_introgressions.R \
  --genome_size 100000000 \
  --n_samples 100 \
  --snp_density 0.01 \
  --out_prefix 05_summary/bc2s3_controlled \
  --seed 42
```

**Population Composition:**
- 50 **mixed samples**: Contains REF, HET, and ALT genotypes
- 25 **ref_alt samples**: Inbred lines (homozygous only)
- 25 **ref_het samples**: Heterozygous carriers

**Outputs:**
- `bc2s3_controlled.vcf.gz` - Multi-sample VCF (1M SNPs × 100 samples)
- `bc2s3_controlled.introgressions.tsv` - Truth file with exact introgression locations
- `bc2s3_controlled.sample_summary.tsv` - Per-sample statistics
- `bc2s3_controlled.introgression_regions.bed` - BED file for visualization

#### Step 3: Generate Sample FASTAs

Convert VCF to diploid FASTA files (one per sample with both haplotypes):

```bash
./scripts/mac/23e_gen_all_fastas.sh
```

This generates:
- `05_summary/bc2s3_fastas/reference.fa` - Reference genome (100 Mb)
- `05_summary/bc2s3_fastas/sample1.fa` through `sample100.fa` - Diploid FASTAs (~200 MB each)

**Technical Details:**
- Uses `bcftools consensus` to apply variants to reference
- Generates `>sampleX_hap1` and `>sampleX_hap2` sequences per sample
- Total output: ~19 GB for 100 samples

#### Step 4: Simulate Reads and Call Genotype Likelihoods

Complete pipeline from FASTA to variant calls with GL/PL values:

```bash
./scripts/mac/25_simulate_reads_call_gl.sh \
  --fasta-dir 05_summary/bc2s3_fastas \
  --out-dir 05_summary/bc2s3_reads \
  --coverages 10,8,4,1,0.5 \
  --threads 4
```

**Pipeline Steps:**
1. Extract haplotypes from diploid FASTA (`seqkit grep`)
2. Simulate Illumina PE reads per haplotype at half coverage (`art_illumina`)
3. Combine reads and align to reference (`bwa mem`)
4. Call variants with genotype likelihoods (`bcftools mpileup + call`)

**Outputs:**
- `bc2s3_reads/bams/sampleX_covY.bam` - Aligned reads
- `bc2s3_reads/vcf/sampleX_covY.vcf.gz` - VCF with PL values

**Example Output (500 jobs = 100 samples × 5 coverages):**
```
[1/500] sample1 @ 10x
  -> 05_summary/bc2s3_reads/vcf/sample1_cov10.vcf.gz
[2/500] sample1 @ 8x
  -> 05_summary/bc2s3_reads/vcf/sample1_cov8.vcf.gz
...
```

**Working with GL/PL Values:**

VCF files contain PL (Phred-scaled likelihoods). To convert to GL (log10):
```
GL = -PL/10
```

Extract GLs:
```bash
bcftools query -f '%CHROM\t%POS\t[%PL\t]\n' file.vcf.gz
```

---

## Folder Layout

Created by scripts as needed:

```
├── 00_ref/              # Reference FASTA files
├── 01_simref/           # Simulated references with repeats
├── 02_reads/            # Simulated reads
├── 03_assemblies/       # Unicycler assemblies
├── 04_eval/             # QUAST evaluation outputs
├── 05_summary/          # Final outputs, VCFs, summaries
│   ├── bc2s3_fastas/    # Sample FASTA files
│   └── bc2s3_reads/     # BAMs and VCFs with GLs
├── inst/extdata/
│   ├── ref_files/       # Downloaded reference files
│   └── panels/          # Haplotype panels
├── ref/                 # Large reference genomes
└── scripts/
    ├── mac/             # Mac scripts
    └── hpc/             # HPC scripts
```

---

## Complete Script Reference

### Genome Simulation
| Script | Description |
|--------|-------------|
| `01_make_tandem_repeats.R` | Create repeat-stress genomes with tandem/motif repeats |
| `01b_generate_annotations.R` | Generate GFF3 annotations, operons, TSS, plasmids |
| `01c_fetch_marker_panel.R` | Download curated marker gene panel |
| `01d_fetch_ref_files.R` | Download reference chromosomes (E. coli, Human, Maize) |

### Read Simulation
| Script | Description |
|--------|-------------|
| `02_sim_illumina_art.R` | Simulate Illumina PE reads (ART) |
| `03_sim_pacbio.R` | Simulate PacBio CLR/HiFi reads (PBSIM) |
| `04_run_grid_both_pacbio.R` | Run coverage grid for Illumina + PacBio |

### Assembly & Evaluation
| Script | Description |
|--------|-------------|
| `05_run_unicycler_grid_both.R` | Hybrid assembly with Unicycler |
| `06_run_quast_both.R` | Evaluate assemblies with QUAST |
| `07_summarize_quast.R` | Aggregate QUAST results to CSV |

### Population Genetics
| Script | Description |
|--------|-------------|
| `10_simulate_gwas_cohort.R` | GWAS cohort with LD, population structure |
| `11_simulate_breeding.R` | F1, BC, RIL, NIL, DH, MAGIC, NAM breeding |
| `11a_generate_random_haplotype_panel.R` | Generate random haplotype panel |
| `12_simupop_api.R` | SimuPOP integration for complex breeding |

### Cross-Species Introgression Pipeline
| Script | Description |
|--------|-------------|
| `20_prepare_aligned_haplotypes.sh` | Align divergent species chromosomes (minimap2) |
| `22_simulate_controlled_introgressions.R` | Generate BC2S3-like populations with known introgressions |
| `23d_fix_reference.sh` | Fix reference FASTA to match VCF REF alleles |
| `23e_gen_all_fastas.sh` | Generate diploid FASTAs for all samples |
| `25_simulate_reads_call_gl.sh` | Full pipeline: FASTA → reads → BAM → VCF with GL |

---

## Detailed Parameter Reference

### 01_make_tandem_repeats.R

```bash
Rscript scripts/mac/01_make_tandem_repeats.R [options]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mode` | `tandem`, `motif`, or `both` | `tandem` |
| `--n_events` | Number of duplication events | 10 |
| `--seg_len` | Duplicated segment length (bp) | 1000 |
| `--copies` | Total copies in tandem | 3 |
| `--motif` | Repeat motif string (e.g., `ATTA`) | - |
| `--motif_repeat` | Repeats per block | 10 |
| `--ploidy` | Number of genome copies | 2 |
| `--ploidy_mode` | `identical` or `diverged` | `identical` |
| `--ploidy_snp_rate` | SNP rate for diverged copies | 0.001 |
| `--random_genome` | Generate random genome | false |
| `--random_length` | Random genome length | 1000000 |
| `--random_gc` | Random genome GC content | 0.5 |

### 22_simulate_controlled_introgressions.R

```bash
Rscript scripts/mac/22_simulate_controlled_introgressions.R [options]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--genome_size` | Simulated genome size (bp) | 100000000 |
| `--n_samples` | Number of samples | 100 |
| `--snp_density` | SNPs per bp | 0.01 |
| `--frac_mixed` | Fraction of mixed samples | 0.5 |
| `--frac_ref_alt` | Fraction of ref_alt samples | 0.25 |
| `--frac_ref_het` | Fraction of ref_het samples | 0.25 |
| `--min_introgression_len` | Minimum introgression length | 100000 |
| `--max_introgression_len` | Maximum introgression length | 5000000 |
| `--donor_fraction` | Target donor genome fraction | 0.125 |

### 25_simulate_reads_call_gl.sh

```bash
./scripts/mac/25_simulate_reads_call_gl.sh [options]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--fasta-dir` | Directory with sample FASTAs | `05_summary/bc2s3_fastas` |
| `--out-dir` | Output directory | `05_summary/bc2s3_reads` |
| `--coverages` | Coverage levels (CSV) | `10,8,4,1,0.5` |
| `--threads` | Number of threads | 4 |
| `--samples` | Specific samples to process | all |

---

## Example Data Sources

### Maize (Panzea HapMap / GBS)
- [Panzea Genotypes](https://www.panzea.org/genotypes) - HapMap, GBS, SNP arrays
- [Panzea Genotype Search](https://www.panzea.org/genotype-search) - HapMap, SNP50, NAM SNPs

### Rice (3,000 Rice Genomes)
- [3K RG on AWS](https://iric.irri.org/resources/3krg-in-aws) - Public dataset access
- [Project Overview](https://iric.irri.org/projects/3000-rice-genomes-project)

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Path with spaces breaks SPAdes/QUAST | Use a no-spaces path |
| `mapfile: command not found` on macOS | Scripts use compatible array syntax |
| bwa/samtools not found | Activate conda: `conda activate simitall` |
| Empty BAM files | Check seqkit haplotype extraction |
| Reference mismatch in bcftools | Run `23d_fix_reference.sh` first |

---

## Citation

If you use these tools in a publication, please cite:
- **ART** - Illumina read simulation
- **PBSIM / PBSIM3** - PacBio read simulation
- **Unicycler** - Hybrid assembly
- **QUAST** - Assembly evaluation
- **bcftools** - Variant calling
- **minimap2** - Sequence alignment
- **SimuPOP** - Population genetics simulation

---

## License

MIT License - See LICENSE file for details.
