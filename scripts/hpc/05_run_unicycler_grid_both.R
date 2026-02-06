#!/usr/bin/env Rscript

# R runner for Unicycler hybrid assemblies

run_unicycler_grid <- function(tag, threads = 8) {
  ill_covs <- c(10, 20, 30, 40, 60, 80)
  pb_covs <- c(5, 10, 15, 20, 30, 40)

  for (ic in ill_covs) {
    r1 <- file.path("02_reads", tag, "illumina", paste0(tag, ".ill_cov", ic, "_1.fq.gz"))
    r2 <- file.path("02_reads", tag, "illumina", paste0(tag, ".ill_cov", ic, "_2.fq.gz"))

    for (pc in pb_covs) {
      # CLR
      pb_dir <- file.path("02_reads", tag, "pacbio_CLR", paste0("cov", pc))
      pb_fastq <- list.files(pb_dir, pattern = paste0("pbCLR_cov", pc, ".*fastq.gz"), full.names = TRUE)
      if (length(pb_fastq) > 0) {
        out_dir <- file.path("03_assemblies", tag, "CLR", paste0("ill", ic, "_pb", pc))
        if (!file.exists(file.path(out_dir, "assembly.fasta"))) {
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          cmd <- sprintf("unicycler -1 %s -2 %s -l %s -o %s -t %s --mode normal",
                         shQuote(r1), shQuote(r2), shQuote(pb_fastq[1]), shQuote(out_dir), threads)
          system(cmd)
        }
      }

      # HIFI
      pb_dir <- file.path("02_reads", tag, "pacbio_HIFI", paste0("cov", pc))
      pb_fastq <- list.files(pb_dir, pattern = paste0("pbHIFI_cov", pc, ".*fastq.gz"), full.names = TRUE)
      if (length(pb_fastq) > 0) {
        out_dir <- file.path("03_assemblies", tag, "HIFI", paste0("ill", ic, "_pb", pc))
        if (!file.exists(file.path(out_dir, "assembly.fasta"))) {
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          cmd <- sprintf("unicycler -1 %s -2 %s -l %s -o %s -t %s --mode normal",
                         shQuote(r1), shQuote(r2), shQuote(pb_fastq[1]), shQuote(out_dir), threads)
          system(cmd)
        }
      }
    }
  }
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  tag <- args[1]
  threads <- ifelse(length(args) >= 2, as.integer(args[2]), 8)
  run_unicycler_grid(tag, threads)
}
