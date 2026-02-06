#!/usr/bin/env Rscript

# R runner for read simulation grid

source("scripts/mac/02_sim_illumina_art.R")
source("scripts/mac/03_sim_pacbio.R")

run_grid_both <- function(simref_fa, tag) {
  ill_covs <- c(10, 20, 30, 40, 60, 80)
  pb_covs <- c(5, 10, 15, 20, 30, 40)

  dir.create(file.path("02_reads", tag, "illumina"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path("02_reads", tag, "pacbio_CLR"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path("02_reads", tag, "pacbio_HIFI"), showWarnings = FALSE, recursive = TRUE)

  for (ic in ill_covs) {
    outp <- file.path("02_reads", tag, "illumina", paste0(tag, ".ill_cov", ic, "_"))
    if (!file.exists(paste0(outp, "1.fq.gz"))) {
      sim_illumina_art(simref_fa, outp, ic, 150, 350, 50)
    }
  }

  for (pc in pb_covs) {
    outd <- file.path("02_reads", tag, "pacbio_CLR", paste0("cov", pc))
    if (!length(list.files(outd, pattern = paste0("pbCLR_cov", pc, ".*fastq.gz"), full.names = TRUE))) {
      sim_pacbio(simref_fa, outd, pc, "CLR", 1)
    }
  }

  for (pc in pb_covs) {
    outd <- file.path("02_reads", tag, "pacbio_HIFI", paste0("cov", pc))
    if (!length(list.files(outd, pattern = paste0("pbHIFI_cov", pc, ".*fastq.gz"), full.names = TRUE))) {
      sim_pacbio(simref_fa, outd, pc, "HIFI", 1)
    }
  }

  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  run_grid_both(args[1], args[2])
}
