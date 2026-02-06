#!/usr/bin/env Rscript

# R wrapper for PacBio simulation (pbsim/pbsim3)

sim_pacbio <- function(ref_fa, outdir, cov, type = "HIFI", seed = 1) {
  if (!file.exists(ref_fa)) stop("Reference FASTA not found: ", ref_fa)
  type <- toupper(type)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # prefer pbsim3 if present
  has_pbsim3 <- nzchar(Sys.which("pbsim3"))
  has_pbsim <- nzchar(Sys.which("pbsim"))

  if (has_pbsim3) {
    prefix <- ifelse(type == "HIFI", "pbHIFI", "pbCLR")
    cmd <- sprintf("pbsim3 --depth %s --seed %s --prefix %s_cov%s %s", cov, seed, prefix, cov, shQuote(ref_fa))
  } else if (has_pbsim) {
    dtype <- ifelse(type == "HIFI", "CCS", "CLR")
    prefix <- ifelse(type == "HIFI", "pbHIFI", "pbCLR")
    cmd <- sprintf("pbsim --depth %s --seed %s --data-type %s --prefix %s_cov%s %s", cov, seed, dtype, prefix, cov, shQuote(ref_fa))
  } else {
    stop("pbsim3 or pbsim not found in PATH")
  }

  old <- getwd(); setwd(outdir); on.exit(setwd(old), add = TRUE)
  status <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (status != 0) stop("PacBio simulation failed")

  # gzip FASTQ outputs
  fq <- list.files(outdir, pattern = "\\.fastq$", full.names = TRUE)
  if (length(fq) > 0) {
    for (f in fq) system(sprintf("pigz -f %s", shQuote(f)))
  }
  invisible(outdir)
}

# CLI entry
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  ref_fa <- args[1]
  outdir <- args[2]
  cov <- as.numeric(args[3])
  type <- args[4]
  seed <- ifelse(length(args) >= 5, as.numeric(args[5]), 1)
  sim_pacbio(ref_fa, outdir, cov, type, seed)
}
