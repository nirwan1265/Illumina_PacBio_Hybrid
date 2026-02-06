#!/usr/bin/env Rscript

# R runner for QUAST evaluation

run_quast_grid <- function(tag) {
  ill_covs <- c(10, 20, 30, 40, 60, 80)
  pb_covs <- c(5, 10, 15, 20, 30, 40)
  ref <- file.path("01_simref", paste0(tag, ".fa"))

  for (mode in c("CLR", "HIFI")) {
    outroot <- file.path("04_eval", tag, mode)
    dir.create(outroot, showWarnings = FALSE, recursive = TRUE)

    for (ic in ill_covs) {
      for (pc in pb_covs) {
        asm <- file.path("03_assemblies", tag, mode, paste0("ill", ic, "_pb", pc), "assembly.fasta")
        if (!file.exists(asm)) next
        outdir <- file.path(outroot, paste0("ill", ic, "_pb", pc))
        if (file.exists(file.path(outdir, "report.tsv"))) next

        cmd <- sprintf("quast.py -r %s -o %s --threads 4 %s",
                       shQuote(ref), shQuote(outdir), shQuote(asm))
        system(cmd)
      }
    }
  }
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  run_quast_grid(args[1])
}
