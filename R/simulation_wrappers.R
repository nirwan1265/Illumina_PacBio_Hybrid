#' Simulate Illumina reads with ART
#'
#' Wrapper around the `art_illumina` command for paired-end read simulation.
#' This function runs ART and compresses output FASTQ files with pigz if available.
#'
#' @param ref_fa Character. Reference FASTA path.
#' @param outprefix Character. Output prefix (ART will write `<outprefix>1.fq` and `<outprefix>2.fq`).
#' @param cov Numeric. Fold coverage.
#' @param readlen Integer. Read length (default 150).
#' @param ins Integer. Mean insert size (default 350).
#' @param sd Integer. Insert size standard deviation (default 50).
#'
#' @return A list with `r1` and `r2` file paths to gzipped FASTQ files.
#' @examples
#' \dontrun{
#' sim_illumina_art("01_simref/ecoli_repMed.fa", "02_reads/ecoli/illumina/ecoli.ill_cov30_", 30)
#' }
#' @export
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
  list(r1 = paste0(r1, ".gz"), r2 = paste0(r2, ".gz"))
}

#' Simulate PacBio reads with PBSIM/PBSIM3
#'
#' Wrapper around `pbsim` or `pbsim3` for CLR/HiFi simulation.
#'
#' @param ref_fa Character. Reference FASTA path.
#' @param outdir Character. Output directory.
#' @param cov Numeric. Fold coverage.
#' @param type Character. `"CLR"` or `"HIFI"`.
#' @param seed Integer. Random seed (default 1).
#'
#' @return Character. Output directory path.
#' @examples
#' \dontrun{
#' sim_pacbio("01_simref/ecoli_repMed.fa", "02_reads/ecoli/pacbio_HIFI/cov20", 20, type = "HIFI")
#' }
#' @export
sim_pacbio <- function(ref_fa, outdir, cov, type = "HIFI", seed = 1) {
  if (!file.exists(ref_fa)) stop("Reference FASTA not found: ", ref_fa)
  type <- toupper(type)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

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

  fq <- list.files(outdir, pattern = "\\.fastq$", full.names = TRUE)
  if (length(fq) > 0) {
    for (f in fq) system(sprintf("pigz -f %s", shQuote(f)))
  }
  outdir
}

#' Run Illumina + PacBio grid simulation
#'
#' Simulates read sets across fixed Illumina and PacBio coverage grids.
#'
#' @param simref_fa Character. Reference FASTA.
#' @param tag Character. Tag used for output subfolders.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_grid_both("01_simref/ecoli_repMed.fa", "ecoli_repMed")
#' }
#' @export
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

  TRUE
}

#' Run Unicycler hybrid assemblies across a grid
#'
#' @param tag Character. Tag used for input/output folders.
#' @param threads Integer. Number of threads.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_unicycler_grid("ecoli_repMed", threads = 8)
#' }
#' @export
run_unicycler_grid <- function(tag, threads = 8) {
  ill_covs <- c(10, 20, 30, 40, 60, 80)
  pb_covs <- c(5, 10, 15, 20, 30, 40)

  for (ic in ill_covs) {
    r1 <- file.path("02_reads", tag, "illumina", paste0(tag, ".ill_cov", ic, "_1.fq.gz"))
    r2 <- file.path("02_reads", tag, "illumina", paste0(tag, ".ill_cov", ic, "_2.fq.gz"))

    for (pc in pb_covs) {
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
  TRUE
}

#' Run QUAST evaluation across assembly grid
#'
#' @param tag Character. Tag used for input/output folders.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_quast_grid("ecoli_repMed")
#' }
#' @export
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
  TRUE
}
