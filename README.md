# Illumina_PacBio_Hybrid

This repo contains a Mac‑friendly and HPC‑friendly **hybrid simulation → assembly → evaluation** pipeline for bacterial genomes (E. coli–like). It simulates **Illumina PE** reads, **PacBio CLR + HiFi/CCS** reads, runs **Unicycler** hybrid assembly, and evaluates assemblies with **QUAST**.

The scripts are sequentially numbered so you always know what to run first.

## Workspace note (Mac)
Use a path **without spaces** for reliable SPAdes/QUAST runs.
Example:
`/Users/nirwantandukar/Documents/Github/Illumina_PacBio_Hybrid`

## Folder layout
These are created by the scripts as needed:
- `00_ref` reference FASTA
- `01_simref` simulated reference with tandem repeats
- `02_reads` simulated reads
- `03_assemblies` assemblies
- `04_eval` QUAST outputs
- `05_summary` CSV summaries
- `scripts/mac` Mac scripts
- `scripts/hpc` HPC scripts

## Scripts (Mac)
Run in order:
1. `scripts/mac/00_setup_conda.sh`
2. `scripts/mac/01_make_tandem_repeats.py`
3. `scripts/mac/02_sim_illumina_art.sh`
4. `scripts/mac/03_sim_pacbio.sh`
5. `scripts/mac/04_run_grid_both_pacbio.sh`
6. `scripts/mac/05_run_unicycler_grid_both.sh`
7. `scripts/mac/06_run_quast_both.sh`
8. `scripts/mac/07_summarize_quast.R`
9. `scripts/mac/08_smoke_test.sh` (optional quick check)

## Scripts (HPC)
The same logic lives in `scripts/hpc` with identical numbers except the conda setup and smoke test.

## Quick start (Mac)
Create a conda env locally in the repo (no global writes):

```bash
scripts/mac/00_setup_conda.sh hybridseq .conda_envs
```

Place a reference genome FASTA at:
`00_ref/ecoli.fa`

Create a repeat‑stress reference:

```bash
scripts/mac/01_make_tandem_repeats.py \
  --in_fa 00_ref/ecoli.fa \
  --out_fa 01_simref/ecoli_repMed.fa \
  --n_events 15 --seg_len 1000 --copies 5 --seed 1
```

Simulate reads (Illumina + CLR + HiFi):

```bash
scripts/mac/04_run_grid_both_pacbio.sh 01_simref/ecoli_repMed.fa ecoli_repMed
```

Assemble hybrid grids:

```bash
scripts/mac/05_run_unicycler_grid_both.sh ecoli_repMed 4
```

Run QUAST:

```bash
scripts/mac/06_run_quast_both.sh ecoli_repMed
```

Summarize:

```bash
scripts/mac/07_summarize_quast.R 04_eval/ecoli_repMed/CLR  05_summary/ecoli_repMed.CLR.quast_summary.csv
scripts/mac/07_summarize_quast.R 04_eval/ecoli_repMed/HIFI 05_summary/ecoli_repMed.HIFI.quast_summary.csv
```

## Smoke test (Mac)
A fast, minimal end‑to‑end run:

```bash
scripts/mac/08_smoke_test.sh hybridseq smoke 2 .conda_envs
```

## Parameter explanations
### 01_make_tandem_repeats.py
- `--n_events` number of tandem repeat insertions
- `--seg_len` length of each duplicated segment (bp)
- `--copies` tandem copy count per event
- `--seed` RNG seed

### 02_sim_illumina_art.sh
Arguments:
- `<ref.fa>` reference FASTA
- `<outprefix>` output prefix
- `<cov>` coverage (fold)
- `[readlen]` read length, default 150
- `[ins]` mean insert size, default 350
- `[sd]` insert SD, default 50

### 03_sim_pacbio.sh
Arguments:
- `<ref.fa>` reference FASTA
- `<outdir>` output directory
- `<cov>` coverage (fold)
- `<CLR|HIFI>` read type
- `[seed]` RNG seed

Notes:
- Supports older `pbsim` (no `--data-type`) and uses bundled model files.
- If `pbsim` outputs `.bam`, the script converts to FASTQ using `samtools`.

### 04_run_grid_both_pacbio.sh
Runs a coverage grid for Illumina + PacBio CLR + PacBio HiFi.
Edit the arrays inside the script to change coverage grids.

### 05_run_unicycler_grid_both.sh
Runs Unicycler hybrid assemblies for each Illumina × PacBio pair.
The Mac smoke test uses conservative settings for speed.

### 06_run_quast_both.sh
Runs QUAST against the known simulated reference.
On macOS, this script copies QUAST to a temp path if the environment path contains spaces.

### 07_summarize_quast.R
Aggregates `report.tsv` files into a single CSV.
If some metrics are missing (e.g., minimap2 didn’t compile), values are `NA`.

## Common issues
- **Path with spaces** can break SPAdes/QUAST. Use a no‑spaces path.
- **macOS arm64** can fail to install old QUAST via conda. Use pip install of QUAST.
- **QUAST minimap2 build warnings** in smoke tests are OK; basic stats still work.

## Citation
If you use these tools in a publication, please cite:
- ART
- PBSIM / PBSIM3
- Unicycler
- QUAST

