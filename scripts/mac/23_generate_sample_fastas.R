#!/usr/bin/env Rscript
# Generate reference and per-sample FASTAs from controlled introgression simulation
# These FASTAs can then be used for read simulation
#
# Usage: Rscript 23_generate_sample_fastas.R --vcf <path> --out_dir <path>

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[idx + 1]
}

vcf_file <- get_arg("--vcf", "05_summary/bc2s3_controlled.vcf")
out_dir <- get_arg("--out_dir", "05_summary/bc2s3_fastas")
genome_size <- as.integer(get_arg("--genome_size", 100000000))
seed <- as.integer(get_arg("--seed", 42))

set.seed(seed)

cat("=== Generate Sample FASTAs ===\n")
cat("VCF:", vcf_file, "\n")
cat("Output dir:", out_dir, "\n")
cat("Genome size:", genome_size, "\n\n")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Generate random reference sequence
cat("Generating reference sequence...\n")
bases <- c("A", "C", "G", "T")
ref_seq <- paste(sample(bases, genome_size, replace = TRUE, prob = c(0.3, 0.2, 0.2, 0.3)), collapse = "")

# Parse VCF to get SNP positions and genotypes
cat("Parsing VCF...\n")
vcf_con <- file(vcf_file, "r")

# Skip header lines
line <- readLines(vcf_con, n = 1)
while (grepl("^##", line)) {
  line <- readLines(vcf_con, n = 1)
}

# Parse column header
header <- strsplit(line, "\t")[[1]]
sample_names <- header[10:length(header)]
n_samples <- length(sample_names)
cat("Found", n_samples, "samples\n")

# Read all variants
cat("Reading variants...\n")
variants <- list()
chunk_size <- 100000
total_read <- 0

while (TRUE) {
  lines <- readLines(vcf_con, n = chunk_size)
  if (length(lines) == 0) break

  for (line in lines) {
    fields <- strsplit(line, "\t")[[1]]
    pos <- as.integer(fields[2])
    ref <- fields[4]
    alt <- fields[5]

    # Extract genotypes (GT is first field before :)
    gts <- sapply(fields[10:length(fields)], function(x) strsplit(x, ":")[[1]][1])

    variants[[length(variants) + 1]] <- list(pos = pos, ref = ref, alt = alt, gts = gts)
  }

  total_read <- total_read + length(lines)
  if (total_read %% 100000 == 0) cat("  Read", total_read, "variants\n")
}
close(vcf_con)

cat("Total variants:", length(variants), "\n")

# Insert variants into reference at their positions
cat("Inserting variants into reference...\n")
ref_chars <- strsplit(ref_seq, "")[[1]]

for (i in seq_along(variants)) {
  v <- variants[[i]]
  ref_chars[v$pos] <- v$ref
}
ref_seq <- paste(ref_chars, collapse = "")

# Write reference FASTA
ref_file <- file.path(out_dir, "reference.fa")
cat("Writing reference:", ref_file, "\n")
write_fasta <- function(path, id, seq, width = 80) {
  con <- file(path, "w")
  writeLines(paste0(">", id), con)
  for (i in seq(1, nchar(seq), by = width)) {
    writeLines(substr(seq, i, min(nchar(seq), i + width - 1)), con)
  }
  close(con)
}
write_fasta(ref_file, "chr1", ref_seq)

# Generate per-sample FASTAs (diploid: two haplotypes)
cat("\nGenerating sample FASTAs...\n")

for (s in 1:n_samples) {
  sample_name <- sample_names[s]

  # Start with reference
  hap1_chars <- ref_chars
  hap2_chars <- ref_chars

  # Apply variants based on genotype
  for (v in variants) {
    gt <- v$gts[s]
    pos <- v$pos

    if (gt == "0/1" || gt == "0|1") {
      # Het: one ref, one alt
      hap1_chars[pos] <- v$ref
      hap2_chars[pos] <- v$alt
    } else if (gt == "1/1" || gt == "1|1") {
      # Hom alt: both alt
      hap1_chars[pos] <- v$alt
      hap2_chars[pos] <- v$alt
    }
    # 0/0 stays as ref (already set)
  }

  hap1_seq <- paste(hap1_chars, collapse = "")
  hap2_seq <- paste(hap2_chars, collapse = "")

  # Write diploid FASTA (both haplotypes for read simulation)
  sample_file <- file.path(out_dir, paste0(sample_name, ".fa"))
  con <- file(sample_file, "w")

  # Write hap1
  writeLines(paste0(">", sample_name, "_hap1"), con)
  for (i in seq(1, nchar(hap1_seq), by = 80)) {
    writeLines(substr(hap1_seq, i, min(nchar(hap1_seq), i + 80 - 1)), con)
  }

  # Write hap2
  writeLines(paste0(">", sample_name, "_hap2"), con)
  for (i in seq(1, nchar(hap2_seq), by = 80)) {
    writeLines(substr(hap2_seq, i, min(nchar(hap2_seq), i + 80 - 1)), con)
  }

  close(con)

  if (s %% 10 == 0) cat("  Generated", s, "/", n_samples, "samples\n")
}

cat("\n=== Done ===\n")
cat("Reference:", ref_file, "\n")
cat("Sample FASTAs:", out_dir, "/sample*.fa\n")
cat("\nNext: Run read simulation with 24_simulate_reads_bcftools.sh\n")
