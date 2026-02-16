# SIMulate IT ALL!!!

## Introduction

This project provides a comprehensive toolkit for **genomic and transcriptomic data simulation**, encompassing everything from bacterial genome benchmarking to complex plant breeding schemes. Whether you need to simulate Illumina short reads, PacBio long reads, cross-species introgression, GWAS cohorts, or complete breeding populations with genotype likelihoods, this toolkit has you covered.

### Simulation Capabilities

| Category | Capabilities |
|----------|-------------|
| **Genome Simulation** | Repeat-stress genomes, random genomes, ploidy simulation |
| **DNA-seq Simulation** | Illumina PE (ART), PacBio CLR/HiFi (PBSIM), coverage grids |
| **GWAS Cohorts** | Population structure, LD blocks, causal variants, phenotypes |
| **Breeding Populations** | F1, BC, RIL, NIL, DH, MAGIC, NAM schemes |
| **Cross-Species Breeding** | Align divergent haplotypes, simulate controlled introgressions |
| **RNA-seq Simulation** | Bulk RNA-seq, scRNA-seq (recommended external tools) |

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

## 1. Genome Simulation

Simulate reference genomes with controlled complexity: tandem repeats, motif repeats, random sequences, and polyploid genomes.

### 1.1 Create Repeat-Stress Genomes

```bash
Rscript scripts/mac/01_make_tandem_repeats.R \
  --in_fa 00_ref/ecoli.fa \
  --out_fa 01_simref/ecoli_repMed.fa \
  --mode tandem \
  --n_events 15 \
  --seg_len 1000 \
  --copies 5
```

#### Parameters: `01_make_tandem_repeats.R`

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--in_fa` | Input reference FASTA | required |
| `--out_fa` | Output FASTA with repeats | required |
| `--mode` | `tandem`, `motif`, or `both` | `tandem` |
| `--n_events` | Number of duplication events | 10 |
| `--seg_len` | Duplicated segment length (bp) | 1000 |
| `--copies` | Total copies per tandem array | 3 |
| `--motif` | Repeat motif string (e.g., `ATTA`) | - |
| `--motif_repeat` | Repeats per motif block | 10 |
| `--ploidy` | Number of genome copies | 2 |
| `--ploidy_mode` | `identical` or `diverged` | `identical` |
| `--ploidy_snp_rate` | SNP rate for diverged copies | 0.001 |
| `--random_genome` | Generate random genome instead | false |
| `--random_length` | Random genome length (bp) | 1000000 |
| `--random_gc` | Random genome GC content | 0.5 |

### 1.2 Generate Annotations

```bash
Rscript scripts/mac/01b_generate_annotations.R \
  --genome_fa 01_simref/ecoli_repMed.fa \
  --out_gff3 01_simref/ecoli_repMed.gff3
```

### 1.3 Download Reference Files

```bash
Rscript scripts/mac/01d_fetch_ref_files.R \
  --species human \
  --chromosome chr1 \
  --out_dir inst/extdata/ref_files
```

---

## 2. DNA-seq Simulation (Genomics)

Simulate sequencing reads from reference genomes using industry-standard simulators.

### 2.1 Illumina Short Reads (ART)

```bash
Rscript scripts/mac/02_sim_illumina_art.R \
  01_simref/ecoli_repMed.fa \
  02_reads/illumina/ecoli_cov30 \
  30 \
  150 \
  350 \
  50
```

#### Parameters: `02_sim_illumina_art.R`

| Position | Parameter | Description | Default |
|----------|-----------|-------------|---------|
| 1 | `ref_fa` | Reference FASTA | required |
| 2 | `outprefix` | Output prefix | required |
| 3 | `cov` | Coverage depth | required |
| 4 | `readlen` | Read length (bp) | 150 |
| 5 | `ins` | Insert size (bp) | 350 |
| 6 | `sd` | Insert size SD | 50 |

**Outputs:** `<prefix>1.fq.gz`, `<prefix>2.fq.gz`

### 2.2 PacBio Long Reads (PBSIM)

```bash
Rscript scripts/mac/03_sim_pacbio.R \
  01_simref/ecoli_repMed.fa \
  02_reads/pacbio_hifi \
  20 \
  HIFI \
  42
```

#### Parameters: `03_sim_pacbio.R`

| Position | Parameter | Description | Default |
|----------|-----------|-------------|---------|
| 1 | `ref_fa` | Reference FASTA | required |
| 2 | `outdir` | Output directory | required |
| 3 | `cov` | Coverage depth | required |
| 4 | `type` | `HIFI` or `CLR` | required |
| 5 | `seed` | Random seed | 1 |

**Outputs:** `<outdir>/pb{HIFI,CLR}_cov<cov>_*.fastq.gz`

### 2.3 Coverage Grid Simulation

Run systematic coverage combinations for benchmarking:

```bash
Rscript scripts/mac/04_run_grid_both_pacbio.R
```

**Default Coverage Levels:**
- Illumina: 10x, 20x, 30x, 40x, 60x, 80x
- PacBio: 5x, 10x, 15x, 20x, 30x, 40x

### 2.4 Hybrid Assembly Benchmarking

Complete pipeline from reference to assembly evaluation:

```bash
# 1. Create repeat-stress reference
Rscript scripts/mac/01_make_tandem_repeats.R --in_fa 00_ref/ecoli.fa --out_fa 01_simref/ecoli_repMed.fa

# 2. Simulate reads (coverage grid)
Rscript scripts/mac/04_run_grid_both_pacbio.R

# 3. Run Unicycler hybrid assembly
Rscript scripts/mac/05_run_unicycler_grid_both.R

# 4. Evaluate with QUAST
Rscript scripts/mac/06_run_quast_both.R

# 5. Summarize results
Rscript scripts/mac/07_summarize_quast.R 04_eval/ecoli_repMed/CLR 05_summary/ecoli_repMed.CLR.quast_summary.csv
```

---

## 3. GWAS Cohort Simulation

Simulate GWAS-style cohorts with population structure, LD blocks, and phenotypes.

### 3.1 Basic GWAS Cohort

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

#### Parameters: `10_simulate_gwas_cohort.R`

**Required:**

| Parameter | Description |
|-----------|-------------|
| `--genome_fa` | Reference FASTA (single-contig preferred) |
| `--out_prefix` | Output prefix (.vcf, .geno.tsv, .pheno.tsv) |

**Cohort:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--n_samples` | Number of samples | 100 |
| `--ploidy` | Ploidy level (1=haploid, 2=diploid) | 2 |

**Variants:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--snp_rate` | Per-bp SNP rate | 0.001 |
| `--indel_rate` | Per-bp indel rate | 0.0001 |
| `--indel_maxlen` | Maximum indel length | 3 |
| `--af_beta1` | Allele-freq Beta alpha | 0.8 |
| `--af_beta2` | Allele-freq Beta beta | 0.8 |

**Population Structure:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--n_pops` | Number of populations | 1 |
| `--pop_sizes` | Population sizes (CSV or fractions) | equal |
| `--fst` | Fst divergence between populations | 0.0 |
| `--pop_effect_shift` | Phenotype shift between populations | 0.0 |

**LD Structure:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--ld_block_size` | LD block size (bp, 0=no LD) | 0 |
| `--ld_haplotypes` | Haplotypes per LD block | 6 |

**Recombination:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--recomb_map_in` | Input recombination map TSV | - |
| `--recomb_map_out` | Output recombination map TSV | - |
| `--recomb_rate_mean` | Mean cM/Mb | 1.0 |
| `--recomb_rate_sd` | SD of recombination rate | 0.3 |
| `--recomb_hotspots` | Number of hotspots | 3 |
| `--recomb_hotspot_mult` | Hotspot multiplier | 5.0 |

**Outputs:**
- `<prefix>.vcf` - VCF with genotypes
- `<prefix>.geno.tsv` - Genotype matrix
- `<prefix>.pheno.tsv` - Phenotype table

### 3.2 Phenotype Simulation

Simulate quantitative or binary traits with causal variants.

#### Phenotype Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--phenotype` | `none`, `binary`, or `quantitative` | quantitative |
| `--n_causal` | Number of causal variants | 20 |
| `--effect_sd` | Effect size SD | 0.5 |
| `--case_frac` | Target case fraction (binary only) | 0.5 |

**Quantitative Trait:**
```bash
Rscript scripts/mac/10_simulate_gwas_cohort.R \
  --genome_fa ref.fa \
  --out_prefix gwas_quant \
  --phenotype quantitative \
  --n_causal 100 \
  --effect_sd 0.3
```

**Binary Trait (Case-Control):**
```bash
Rscript scripts/mac/10_simulate_gwas_cohort.R \
  --genome_fa ref.fa \
  --out_prefix gwas_binary \
  --phenotype binary \
  --n_causal 50 \
  --effect_sd 0.8 \
  --case_frac 0.3
```

---

## 4. Breeding Population Simulation

Simulate F1, selfing, backcross, and advanced breeding schemes from haplotype panels.

### 4.1 Supported Breeding Schemes

| Scheme | Description |
|--------|-------------|
| **F1** | First filial generation from two parents |
| **SELF** | Self-fertilization (SSD) |
| **SIB** | Sib-mating |
| **BC** | Backcross to P1 or P2 |
| **RIL** | Recombinant Inbred Lines |
| **NIL** | Near-Isogenic Lines |
| **DH** | Doubled Haploids |
| **MAGIC** | Multi-parent Advanced Generation Inter-Cross |
| **NAM** | Nested Association Mapping |

### 4.2 Basic Breeding Simulation

**F1 → Self 2 generations → Backcross:**
```bash
Rscript scripts/mac/11_simulate_breeding.R \
  --haplotype_fa inst/extdata/ref_files/random_panel.fa \
  --out_prefix 05_summary/breeding_demo \
  --parents hap1,hap2 \
  --sequence F1,SELF:2,BC:P1:1 \
  --n_offspring 100 \
  --seed 7
```

#### Parameters: `11_simulate_breeding.R`

**Required:**

| Parameter | Description |
|-----------|-------------|
| `--haplotype_fa` | Multi-FASTA panel (same length haplotypes) |
| `--out_prefix` | Output prefix (.fa, .meta.tsv) |

**Parents/Population:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--parents` | Two haplotype IDs or indices (CSV) | random |
| `--n_offspring` | Number of offspring | 100 |
| `--founders` | Founder hap IDs for MAGIC/NAM | - |
| `--n_founders` | Number of founders | 4 |

**Breeding Sequence:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--sequence` | Breeding sequence (e.g., `F1,SELF:3,BC:P1:2`) | F1,SELF:1 |
| `--scheme` | High-level scheme: `F2`, `MAGIC`, `NAM`, `RIL`, `NIL`, `DH` | - |
| `--self_generations` | Selfing generations for RIL/NIL | 6 |
| `--backcross_generations` | BC generations for NIL | 3 |
| `--ril_mating` | `SSD` or `SIB` | SSD |

**Selection:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--background_selection` | Select NILs by minimal donor background | false |
| `--selection_pool` | Candidates per generation | 50 |
| `--marker_step` | Marker spacing for selection (bp) | 1000 |
| `--introgression_target_len` | Target donor tract length (bp) | - |
| `--selection_loci` | Selected loci positions (CSV) | - |
| `--selection_model` | `add`, `dom`, or `rec` | add |
| `--selection_strength` | Fitness penalty per risk allele | 0 |
| `--distortion_rate` | Segregation distortion rate | 0 |

**Fixed Locus:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--fix_locus` | Force locus from donor/recipient (start:end) | - |
| `--fix_allele` | `donor`, `recipient`, or haplotype ID | - |

**Recombination:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--recomb_map_in` | Recombination map TSV (pos_bp, cM) | - |
| `--recomb_rate_mean` | Mean cM/Mb | 1.0 |
| `--recomb_rate_sd` | SD of recombination rate | 0.3 |
| `--recomb_hotspots` | Number of hotspots | 3 |
| `--recomb_hotspot_mult` | Hotspot multiplier | 5.0 |
| `--interference_shape` | Gamma shape for crossover interference | 1.0 |

**Genotype Errors:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--genotype_error` | Per-genotype error rate | 0 |
| `--missing_rate` | Per-genotype missingness | 0 |
| `--sv_rate` | Structural variant rate | 0 |
| `--sv_maxlen` | Max SV length | 1000 |

**Outputs:**

| Parameter | Description |
|-----------|-------------|
| `--vcf_out` | Output VCF path |
| `--graph_out` | Mermaid crossing graph path |
| `--graph_format` | `mmd`, `svg`, or `png` |
| `--qc_out` | QC report TSV |
| `--ascertainment` | `founders` or `all` |

### 4.3 Advanced Schemes

**MAGIC:**
```bash
Rscript scripts/mac/11_simulate_breeding.R \
  --haplotype_fa inst/extdata/ref_files/random_panel.fa \
  --out_prefix 05_summary/magic_demo \
  --scheme MAGIC \
  --n_founders 6 \
  --n_offspring 200
```

**NAM:**
```bash
Rscript scripts/mac/11_simulate_breeding.R \
  --haplotype_fa inst/extdata/ref_files/random_panel.fa \
  --out_prefix 05_summary/nam_demo \
  --scheme NAM \
  --founders hap1,hap2,hap3,hap4,hap5 \
  --n_offspring 200
```

**RIL with Background Selection:**
```bash
Rscript scripts/mac/11_simulate_breeding.R \
  --haplotype_fa panel.fa \
  --out_prefix 05_summary/ril_selected \
  --scheme RIL \
  --background_selection \
  --selection_pool 100 \
  --introgression_target_len 500000
```

### 4.4 Cross-Species Introgression Simulation

Simulate BC2S3-like populations with controlled introgressions from divergent species.

**Use Case:** Introgression breeding from wild relatives (e.g., B73 maize × teosinte)

#### Step 1: Align Divergent Haplotypes

```bash
./scripts/mac/20_prepare_aligned_haplotypes.sh \
  --ref1 ref/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
  --ref2 ref/Zd-Momo-REFERENCE-PanAnd-1.0.fa.gz \
  --chr1 chr10 \
  --chr2 chr10 \
  --out-dir 05_summary/aligned_haps \
  --target-len 100000000
```

#### Step 2: Simulate Controlled Introgressions

```bash
Rscript scripts/mac/22_simulate_controlled_introgressions.R \
  --genome_size 100000000 \
  --n_samples 100 \
  --snp_density 0.01 \
  --out_prefix 05_summary/bc2s3_controlled \
  --seed 42
```

**Parameters: `22_simulate_controlled_introgressions.R`**

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

**Population Composition:**
- 50% **mixed samples**: Contains REF, HET, and ALT genotypes
- 25% **ref_alt samples**: Inbred lines (homozygous only)
- 25% **ref_het samples**: Heterozygous carriers

#### Step 3: Generate Sample FASTAs

```bash
./scripts/mac/23e_gen_all_fastas.sh
```

#### Step 4: Simulate Reads and Call Genotype Likelihoods

```bash
./scripts/mac/25_simulate_reads_call_gl.sh \
  --fasta-dir 05_summary/bc2s3_fastas \
  --out-dir 05_summary/bc2s3_reads \
  --coverages 10,8,4,1,0.5 \
  --threads 4
```

**Parameters: `25_simulate_reads_call_gl.sh`**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--fasta-dir` | Directory with sample FASTAs | 05_summary/bc2s3_fastas |
| `--out-dir` | Output directory | 05_summary/bc2s3_reads |
| `--coverages` | Coverage levels (CSV) | 10,8,4,1,0.5 |
| `--threads` | Number of threads | 4 |
| `--samples` | Specific samples to process | all |

**Outputs:**
- `bc2s3_controlled.vcf.gz` - Multi-sample VCF
- `bc2s3_controlled.introgressions.tsv` - Truth file with exact introgression locations
- `bc2s3_controlled.sample_summary.tsv` - Per-sample statistics
- `bc2s3_reads/vcf/sampleX_covY.vcf.gz` - VCF with PL values

---

## 5. RNA-seq Simulation (Transcriptomics)

For RNA-seq simulation, we recommend using established external tools that integrate well with this pipeline.

### 5.1 Bulk RNA-seq Simulation

#### Recommended Tools

| Tool | Description | Installation |
|------|-------------|--------------|
| **Polyester** | R/Bioconductor RNA-seq simulator | `BiocManager::install("polyester")` |
| **RSEM** | RNA-seq simulator with RSEM models | `conda install -c bioconda rsem` |
| **Flux Simulator** | Full transcriptome simulation | [download](http://sammeth.net/confluence/display/SIM/Home) |

#### Polyester Example

```r
# Install
BiocManager::install("polyester")

library(polyester)

# Simulate from reference transcriptome
simulate_experiment(
  fasta = "transcriptome.fa",
  outdir = "simulated_reads",
  num_reps = c(3, 3),  # 3 reps per condition
  reads_per_transcript = 300,
  fold_changes = fold_changes_matrix,
  paired = TRUE,
  readlen = 100
)
```

#### RSEM Simulation

```bash
# Simulate reads from RSEM model
rsem-simulate-reads \
  reference/rsem_ref \
  estimated_model.stat/rsem_ref.model \
  estimated_isoforms.results \
  0.05 \
  10000000 \
  simulated_reads
```

### 5.2 Single-Cell RNA-seq Simulation

#### Recommended Tools

| Tool | Description | Installation |
|------|-------------|--------------|
| **Splatter** | R/Bioconductor scRNA-seq simulator | `BiocManager::install("splatter")` |
| **scDesign2** | Statistical scRNA-seq simulator | `devtools::install_github("JSB-UCLA/scDesign2")` |
| **SymSim** | Symmetric scRNA-seq simulator | `devtools::install_github("YosefLab/SymSim")` |
| **SERGIO** | GRN-based scRNA-seq simulator | `pip install sergio` |

#### Splatter Example

```r
# Install
BiocManager::install("splatter")

library(splatter)

# Simulate single-cell data
params <- newSplatParams(
  nGenes = 10000,
  batchCells = 1000,
  group.prob = c(0.3, 0.3, 0.4),
  de.prob = 0.1,
  de.facLoc = 0.3
)

sim <- splatSimulate(params, method = "groups")

# Extract counts
counts <- counts(sim)
```

#### scDesign2 Example

```r
library(scDesign2)

# Fit model to reference data
copula_result <- fit_model_scDesign2(
  reference_counts,
  cell_type_labels,
  marginal = "nb"
)

# Simulate new cells
sim_counts <- simulate_count_scDesign2(
  copula_result,
  n_cell_new = 5000
)
```

### 5.3 Integration with This Pipeline

After simulating RNA-seq reads with external tools, you can use this pipeline for:

1. **Alignment and quantification** via standard tools (STAR, Salmon, kallisto)
2. **Downstream analysis** with simulated ground truth

```bash
# Example: Align Polyester output
STAR --genomeDir star_index \
  --readFilesIn simulated_reads/sample_01_1.fasta simulated_reads/sample_01_2.fasta \
  --outFileNamePrefix aligned/sample_01_
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

### Genome Simulation Scripts
| Script | Description |
|--------|-------------|
| `01_make_tandem_repeats.R` | Create repeat-stress genomes with tandem/motif repeats |
| `01b_generate_annotations.R` | Generate GFF3 annotations, operons, TSS, plasmids |
| `01c_fetch_marker_panel.R` | Download curated marker gene panel |
| `01d_fetch_ref_files.R` | Download reference chromosomes (E. coli, Human, Maize) |

### DNA-seq Simulation Scripts
| Script | Description |
|--------|-------------|
| `02_sim_illumina_art.R` | Simulate Illumina PE reads (ART) |
| `03_sim_pacbio.R` | Simulate PacBio CLR/HiFi reads (PBSIM) |
| `04_run_grid_both_pacbio.R` | Run coverage grid for Illumina + PacBio |

### Assembly & Evaluation Scripts
| Script | Description |
|--------|-------------|
| `05_run_unicycler_grid_both.R` | Hybrid assembly with Unicycler |
| `06_run_quast_both.R` | Evaluate assemblies with QUAST |
| `07_summarize_quast.R` | Aggregate QUAST results to CSV |

### GWAS & Population Genetics Scripts
| Script | Description |
|--------|-------------|
| `10_simulate_gwas_cohort.R` | GWAS cohort with LD, population structure, phenotypes |
| `11_simulate_breeding.R` | F1, BC, RIL, NIL, DH, MAGIC, NAM breeding |
| `11a_generate_random_haplotype_panel.R` | Generate random haplotype panel |
| `12_simupop_api.R` | SimuPOP integration for complex breeding |

### Cross-Species Introgression Scripts
| Script | Description |
|--------|-------------|
| `20_prepare_aligned_haplotypes.sh` | Align divergent species chromosomes (minimap2) |
| `22_simulate_controlled_introgressions.R` | Generate BC2S3-like populations with known introgressions |
| `23d_fix_reference.sh` | Fix reference FASTA to match VCF REF alleles |
| `23e_gen_all_fastas.sh` | Generate diploid FASTAs for all samples |
| `25_simulate_reads_call_gl.sh` | Full pipeline: FASTA → reads → BAM → VCF with GL |

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
- **Polyester** - RNA-seq simulation (if used)
- **Splatter** - scRNA-seq simulation (if used)

---

## License

MIT License - See LICENSE file for details.
