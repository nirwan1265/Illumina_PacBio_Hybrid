# Helpers -----------------------------------------------------------------

script_path <- function(fname, script_dir = file.path("scripts", "mac")) {
  path <- file.path(script_dir, fname)
  if (!file.exists(path)) stop("Script not found: ", path)
  path
}

build_cli_args <- function(..., .args = list()) {
  args <- c(.args, list(...))
  out <- character(0)
  for (nm in names(args)) {
    val <- args[[nm]]
    flag <- paste0("--", nm)
    if (is.logical(val)) {
      if (isTRUE(val)) out <- c(out, flag)
    } else if (length(val) > 1) {
      for (v in val) out <- c(out, flag, as.character(v))
    } else if (!is.null(val)) {
      out <- c(out, flag, as.character(val))
    }
  }
  out
}

run_rscript <- function(script, args = character(0)) {
  status <- system2("Rscript", c(script, args), stdout = "", stderr = "")
  if (!identical(status, 0L)) stop("Rscript failed: ", script)
  invisible(TRUE)
}

#' Simulate a genome with tandem/motif repeats (script wrapper)
#'
#' Wrapper around `scripts/mac/01_make_tandem_repeats.R`. Use this to insert
#' tandem repeats, motif repeats, or both into an input FASTA, or to generate
#' a random genome with repeat insertions.
#'
#' @param out_fa Character. Output FASTA path (required).
#' @param in_fa Character. Input FASTA path (optional unless `random_genome=TRUE`).
#' @param random_genome Logical. If TRUE, generate a random genome instead of using `in_fa`.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed as `--key value` to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' simulate_genome(out_fa = "01_simref/ecoli_repMed.fa", in_fa = "00_ref/ecoli.fa", n_events = 15, seg_len = 1000, copies = 5)
#' }
#' @export
simulate_genome <- function(out_fa, in_fa = NULL, random_genome = FALSE, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("01_make_tandem_repeats.R", script_dir)
  args <- build_cli_args(
    out_fa = out_fa,
    in_fa = in_fa,
    random_genome = random_genome,
    .args = list(...)
  )
  run_rscript(script, args)
}

#' Generate annotations on a genome (script wrapper)
#'
#' Wrapper around `scripts/mac/01b_generate_annotations.R` to add genes, operons,
#' TSS/promoters, rRNA/tRNA, CRISPR arrays, riboswitches, terminators, and plasmids.
#'
#' @param in_fa Character. Input FASTA path.
#' @param out_prefix Character. Output prefix for GFF3/TSV/FASTA outputs.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' simulate_annotations(in_fa = "01_simref/ecoli_repMed.fa", out_prefix = "sim_annotations/ecoli_repMed")
#' }
#' @export
simulate_annotations <- function(in_fa, out_prefix, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("01b_generate_annotations.R", script_dir)
  args <- build_cli_args(
    in_fa = in_fa,
    out_prefix = out_prefix,
    .args = list(...)
  )
  run_rscript(script, args)
}

#' Fetch curated marker panel FASTA (script wrapper)
#'
#' Wrapper around `scripts/mac/01c_fetch_marker_panel.R` to download or update
#' marker sequences used for plasmid marker simulation.
#'
#' @param out_fa Character. Output FASTA path.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' fetch_marker_panel(out_fa = "inst/extdata/marker_panel.fasta")
#' }
#' @export
fetch_marker_panel <- function(out_fa, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("01c_fetch_marker_panel.R", script_dir)
  args <- build_cli_args(out_fa = out_fa, .args = list(...))
  run_rscript(script, args)
}

#' Fetch example reference FASTA files (script wrapper)
#'
#' Wrapper around `scripts/mac/01d_fetch_ref_files.R` to download example
#' bacterial, human, and maize reference contigs into `inst/extdata/ref_files`.
#'
#' @param out_dir Character. Output directory.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' fetch_ref_files(out_dir = "inst/extdata/ref_files")
#' }
#' @export
fetch_ref_files <- function(out_dir, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("01d_fetch_ref_files.R", script_dir)
  args <- build_cli_args(out_dir = out_dir, .args = list(...))
  run_rscript(script, args)
}

#' Summarize QUAST reports into a single CSV (script wrapper)
#'
#' Wrapper around `scripts/mac/07_summarize_quast.R`.
#'
#' @param quast_root Character. Directory containing per-assembly QUAST `report.tsv` files.
#' @param out_csv Character. Output CSV path.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' summarize_quast(quast_root = "04_eval/ecoli_repMed/HIFI", out_csv = "05_summary/ecoli_repMed.HIFI.csv")
#' }
#' @export
summarize_quast <- function(quast_root, out_csv, script_dir = file.path("scripts", "mac")) {
  script <- script_path("07_summarize_quast.R", script_dir)
  args <- c(quast_root, out_csv)
  run_rscript(script, args)
}

#' Simulate a GWAS cohort with LD blocks (script wrapper)
#'
#' Wrapper around `scripts/mac/10_simulate_gwas_cohort.R`.
#'
#' @param out_prefix Character. Output prefix for VCF/phenotype/QC files.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' simulate_gwas_cohort(out_prefix = "gwas/out", n_samples = 200, n_markers = 5000)
#' }
#' @export
simulate_gwas_cohort <- function(out_prefix, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("10_simulate_gwas_cohort.R", script_dir)
  args <- build_cli_args(out_prefix = out_prefix, .args = list(...))
  run_rscript(script, args)
}

#' Generate a random haplotype panel (script wrapper)
#'
#' Wrapper around `scripts/mac/11a_generate_random_haplotype_panel.R`.
#'
#' @param out_fa Character. Output FASTA path.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' generate_random_haplotype_panel(out_fa = "inst/extdata/panels/panel.fa", n_haplotypes = 20, length = 50000)
#' }
#' @export
generate_random_haplotype_panel <- function(out_fa, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("11a_generate_random_haplotype_panel.R", script_dir)
  args <- build_cli_args(out_fa = out_fa, .args = list(...))
  run_rscript(script, args)
}

#' Simulate breeding schemes (script wrapper)
#'
#' Wrapper around `scripts/mac/11_simulate_breeding.R` with support for
#' biparental, RIL, NIL, NAM, MAGIC, and DH designs.
#'
#' @param out_prefix Character. Output prefix for VCF and QC outputs.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' simulate_breeding(out_prefix = "breeding/out", scheme = "RIL", n_lines = 200, n_generations = 6)
#' }
#' @export
simulate_breeding <- function(out_prefix, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("11_simulate_breeding.R", script_dir)
  args <- build_cli_args(out_prefix = out_prefix, .args = list(...))
  run_rscript(script, args)
}

#' Run SimuPOP API helpers (script wrapper)
#'
#' Wrapper around `scripts/mac/12_simupop_api.R` for simplified SimuPOP presets
#' and VCF export. Requires the `simupop` Python package available to reticulate.
#'
#' @param out_prefix Character. Output prefix for VCF and metadata.
#' @param script_dir Character. Directory containing mac scripts (default `scripts/mac`).
#' @param ... Additional named options passed to the script.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' simupop_api(out_prefix = "simupop/out", scheme = "F2", n_ind = 200)
#' }
#' @export
simupop_api <- function(out_prefix, script_dir = file.path("scripts", "mac"), ...) {
  script <- script_path("12_simupop_api.R", script_dir)
  args <- build_cli_args(out_prefix = out_prefix, .args = list(...))
  run_rscript(script, args)
}
