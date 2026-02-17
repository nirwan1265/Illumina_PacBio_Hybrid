#!/usr/bin/env Rscript
# Mac wrapper for Nanopore read simulation using Badread
# Usage: Rscript scripts/mac/03b_sim_nanopore.R <ref_fa> <outdir> <coverage> [options]

args <- commandArgs(trailingOnly = TRUE)
cli_path <- file.path(dirname(dirname(dirname(sys.frame(1)$ofile))), "inst", "cli", "cli_03b_sim_nanopore.R")

# Source and run CLI
if (file.exists(cli_path)) {
  source(cli_path)
  main_03b_sim_nanopore(args)
} else {
  # Fallback: run directly
  source(file.path(system.file(package = "simitall"), "cli", "cli_03b_sim_nanopore.R"))
  main_03b_sim_nanopore(args)
}
