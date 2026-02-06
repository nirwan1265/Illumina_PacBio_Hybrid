main_02_sim_illumina_art <- function(args = commandArgs(trailingOnly = TRUE)) {

# R wrapper for ART Illumina simulation (direct CLI)

sim_illumina_art <- function(ref_fa, outprefix, cov, readlen = 150, ins = 350, sd = 50) {
  if (!file.exists(ref_fa)) stop("Reference FASTA not found: ", ref_fa)
  cmd <- sprintf(
    "art_illumina -ss HS25 -i %s -p -l %s -f %s -m %s -s %s -o %s",
    shQuote(ref_fa), readlen, cov, ins, sd, shQuote(outprefix)
  )
  status <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (status != 0) stop("ART simulation failed")

  r1 <- paste0(outprefix, "1.fq")
  r2 <- paste0(outprefix, "2.fq")
  if (file.exists(r1)) system(sprintf("pigz -f %s", shQuote(r1)))
  if (file.exists(r2)) system(sprintf("pigz -f %s", shQuote(r2)))
  invisible(list(r1 = paste0(r1, ".gz"), r2 = paste0(r2, ".gz")))
}

# CLI entry
args <- args
if (length(args) >= 3) {
  ref_fa <- args[1]
  outprefix <- args[2]
  cov <- as.numeric(args[3])
  readlen <- ifelse(length(args) >= 4, as.numeric(args[4]), 150)
  ins <- ifelse(length(args) >= 5, as.numeric(args[5]), 350)
  sd <- ifelse(length(args) >= 6, as.numeric(args[6]), 50)
  sim_illumina_art(ref_fa, outprefix, cov, readlen, ins, sd)
}
}
