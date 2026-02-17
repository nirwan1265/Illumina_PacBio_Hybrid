#!/usr/bin/env Rscript
# Simulate BC2S3-like population with controlled introgressions
# Creates known introgression locations for benchmarking HMM/detection methods
#
# Usage: Rscript 22_simulate_controlled_introgressions.R [options]

# Simple argument parsing without optparse
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[idx + 1]
}

opt <- list(
  out_prefix = get_arg("--out_prefix", "05_summary/bc2s3_controlled"),
  genome_size = as.integer(get_arg("--genome_size", 100000000)),
  n_samples = as.integer(get_arg("--n_samples", 100)),
  snp_density = as.numeric(get_arg("--snp_density", 0.01)),
  donor_frac = as.numeric(get_arg("--donor_frac", 0.125)),
  n_introgressions_mean = as.numeric(get_arg("--n_introgressions_mean", 3)),
  n_introgressions_sd = as.numeric(get_arg("--n_introgressions_sd", 1.5)),
  introgression_size_mean = as.numeric(get_arg("--introgression_size_mean", 4000000)),
  introgression_size_sd = as.numeric(get_arg("--introgression_size_sd", 2000000)),
  introgression_size_min = as.numeric(get_arg("--introgression_size_min", 500000)),
  introgression_size_max = as.numeric(get_arg("--introgression_size_max", 15000000)),
  frac_mixed = as.numeric(get_arg("--frac_mixed", 0.5)),
  frac_ref_alt = as.numeric(get_arg("--frac_ref_alt", 0.25)),
  frac_ref_het = as.numeric(get_arg("--frac_ref_het", 0.25)),
  seed = as.integer(get_arg("--seed", 42))
)
set.seed(opt$seed)

cat("=== Controlled Introgression Simulation ===\n")
cat("Genome size:", opt$genome_size / 1e6, "Mb\n")
cat("Samples:", opt$n_samples, "\n")
cat("SNP density:", opt$snp_density, "\n")
cat("Donor fraction:", opt$donor_frac, "\n")
cat("Sample types: mixed=", opt$frac_mixed, ", ref_alt=", opt$frac_ref_alt,
    ", ref_het=", opt$frac_ref_het, "\n", sep = "")
cat("\n")

# Create output directory
dir.create(dirname(opt$out_prefix), showWarnings = FALSE, recursive = TRUE)

# Generate SNP positions
n_snps <- as.integer(opt$genome_size * opt$snp_density)
snp_positions <- sort(sample(1:opt$genome_size, n_snps, replace = FALSE))
cat("Generated", n_snps, "SNP positions\n")

# Generate reference and alternate alleles
bases <- c("A", "C", "G", "T")
ref_alleles <- sample(bases, n_snps, replace = TRUE)
alt_alleles <- sapply(ref_alleles, function(r) sample(setdiff(bases, r), 1))

# Assign sample types
n_mixed <- round(opt$n_samples * opt$frac_mixed)
n_ref_alt <- round(opt$n_samples * opt$frac_ref_alt)
n_ref_het <- opt$n_samples - n_mixed - n_ref_alt

sample_types <- c(rep("mixed", n_mixed), rep("ref_alt", n_ref_alt), rep("ref_het", n_ref_het))
sample_types <- sample(sample_types)  # shuffle

cat("Sample types: mixed=", n_mixed, ", ref_alt=", n_ref_alt, ", ref_het=", n_ref_het, "\n")

# Function to generate introgressions for a sample
generate_introgressions <- function(genome_size, n_mean, n_sd, size_mean, size_sd, size_min, size_max) {
  n_intro <- max(1, round(rnorm(1, n_mean, n_sd)))

  intros <- data.frame(start = integer(), end = integer(), stringsAsFactors = FALSE)
  attempts <- 0
  max_attempts <- 1000

  while (nrow(intros) < n_intro && attempts < max_attempts) {
    size <- max(size_min, min(size_max, round(rnorm(1, size_mean, size_sd))))
    start <- sample(1:(genome_size - size), 1)
    end <- start + size - 1

    # Check for overlap with existing introgressions
    overlaps <- any(sapply(1:nrow(intros), function(i) {
      !(end < intros$start[i] || start > intros$end[i])
    }))

    if (nrow(intros) == 0 || !overlaps) {
      intros <- rbind(intros, data.frame(start = start, end = end))
    }
    attempts <- attempts + 1
  }

  intros[order(intros$start), ]
}

# Generate genotypes for all samples
genotypes <- matrix("0/0", nrow = n_snps, ncol = opt$n_samples)
colnames(genotypes) <- paste0("sample", 1:opt$n_samples)

# Store introgression truth
introgression_truth <- data.frame(
  sample = character(),
  introgression_id = integer(),
  start = integer(),
  end = integer(),
  size = integer(),
  genotype = character(),  # HET or HOM_ALT
  n_snps = integer(),
  stringsAsFactors = FALSE
)

cat("\nGenerating introgressions...\n")

for (i in 1:opt$n_samples) {
  sample_name <- paste0("sample", i)
  sample_type <- sample_types[i]

  # Generate introgressions
  intros <- generate_introgressions(
    opt$genome_size,
    opt$n_introgressions_mean,
    opt$n_introgressions_sd,
    opt$introgression_size_mean,
    opt$introgression_size_sd,
    opt$introgression_size_min,
    opt$introgression_size_max
  )

  for (j in 1:nrow(intros)) {
    start <- intros$start[j]
    end <- intros$end[j]

    # Find SNPs in this introgression
    snp_idx <- which(snp_positions >= start & snp_positions <= end)

    if (length(snp_idx) == 0) next

    # Determine genotype based on sample type
    if (sample_type == "mixed") {
      # Randomly choose HET or HOM_ALT for each introgression
      geno_type <- sample(c("HET", "HOM_ALT"), 1, prob = c(0.4, 0.6))
    } else if (sample_type == "ref_alt") {
      # Only HOM_ALT (inbred)
      geno_type <- "HOM_ALT"
    } else {
      # Only HET
      geno_type <- "HET"
    }

    # Assign genotypes
    if (geno_type == "HET") {
      genotypes[snp_idx, i] <- "0/1"
    } else {
      genotypes[snp_idx, i] <- "1/1"
    }

    # Record truth
    introgression_truth <- rbind(introgression_truth, data.frame(
      sample = sample_name,
      introgression_id = j,
      start = start,
      end = end,
      size = end - start + 1,
      genotype = geno_type,
      n_snps = length(snp_idx),
      stringsAsFactors = FALSE
    ))
  }

  if (i %% 20 == 0) cat("  Processed", i, "/", opt$n_samples, "samples\n")
}

cat("Generated", nrow(introgression_truth), "introgressions total\n")

# Calculate summary stats
cat("\n=== Summary Statistics ===\n")

for (i in 1:opt$n_samples) {
  geno_counts <- table(genotypes[, i])
  ref_frac <- as.numeric(geno_counts["0/0"]) / n_snps
  het_frac <- ifelse("0/1" %in% names(geno_counts), as.numeric(geno_counts["0/1"]) / n_snps, 0)
  alt_frac <- ifelse("1/1" %in% names(geno_counts), as.numeric(geno_counts["1/1"]) / n_snps, 0)
}

# Per-sample summary
sample_summary <- data.frame(
  sample = colnames(genotypes),
  type = sample_types,
  n_ref = apply(genotypes, 2, function(x) sum(x == "0/0")),
  n_het = apply(genotypes, 2, function(x) sum(x == "0/1")),
  n_alt = apply(genotypes, 2, function(x) sum(x == "1/1")),
  stringsAsFactors = FALSE
)
sample_summary$frac_ref <- sample_summary$n_ref / n_snps
sample_summary$frac_het <- sample_summary$n_het / n_snps
sample_summary$frac_alt <- sample_summary$n_alt / n_snps
sample_summary$frac_donor <- (sample_summary$n_het * 0.5 + sample_summary$n_alt) / n_snps

cat("Mean donor fraction:", round(mean(sample_summary$frac_donor), 4), "\n")
cat("Mean REF fraction:", round(mean(sample_summary$frac_ref), 4), "\n")
cat("Mean HET fraction:", round(mean(sample_summary$frac_het), 4), "\n")
cat("Mean ALT fraction:", round(mean(sample_summary$frac_alt), 4), "\n")

# Write VCF
vcf_file <- paste0(opt$out_prefix, ".vcf")
cat("\nWriting VCF:", vcf_file, "\n")

vcf_con <- file(vcf_file, "w")
writeLines("##fileformat=VCFv4.2", vcf_con)
writeLines(paste0("##contig=<ID=chr1,length=", opt$genome_size, ">"), vcf_con)
writeLines("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", vcf_con)
writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", vcf_con)
writeLines("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">", vcf_con)

# Sample metadata
for (i in 1:opt$n_samples) {
  writeLines(paste0("##SAMPLE=<ID=sample", i, ",TYPE=", sample_types[i], ">"), vcf_con)
}

# Header
header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                  colnames(genotypes)), collapse = "\t")
writeLines(header, vcf_con)

# Data rows (write in chunks for speed)
chunk_size <- 10000
n_chunks <- ceiling(n_snps / chunk_size)

for (chunk in 1:n_chunks) {
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, n_snps)

  for (idx in start_idx:end_idx) {
    gt_str <- paste(paste0(genotypes[idx, ], ":30"), collapse = "\t")
    row <- paste(c("chr1", snp_positions[idx], paste0("snp", idx),
                   ref_alleles[idx], alt_alleles[idx], ".", "PASS", "DP=100", "GT:DP",
                   gt_str), collapse = "\t")
    writeLines(row, vcf_con)
  }

  if (chunk %% 10 == 0) cat("  Written", end_idx, "/", n_snps, "variants\n")
}
close(vcf_con)

# Write introgression truth file
truth_file <- paste0(opt$out_prefix, ".introgressions.tsv")
cat("Writing truth file:", truth_file, "\n")
write.table(introgression_truth, truth_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Write sample summary
summary_file <- paste0(opt$out_prefix, ".sample_summary.tsv")
cat("Writing sample summary:", summary_file, "\n")
write.table(sample_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Write BED file of introgressions (for visualization)
bed_file <- paste0(opt$out_prefix, ".introgressions.bed")
cat("Writing BED file:", bed_file, "\n")
bed_data <- data.frame(
  chrom = "chr1",
  start = introgression_truth$start - 1,  # BED is 0-based
  end = introgression_truth$end,
  name = paste0(introgression_truth$sample, "_intro", introgression_truth$introgression_id),
  score = ifelse(introgression_truth$genotype == "HOM_ALT", 1000, 500),
  strand = ".",
  stringsAsFactors = FALSE
)
write.table(bed_data, bed_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Generate per-sample genotype likelihoods (GL) file for HMM input
# Format: chr pos ref alt GL_RR GL_RA GL_AA (log10 scale)
cat("\nGenerating genotype likelihoods for HMM input...\n")

gl_dir <- paste0(opt$out_prefix, "_GL")
dir.create(gl_dir, showWarnings = FALSE, recursive = TRUE)

# Error rate for GL simulation
error_rate <- 0.01

for (i in 1:opt$n_samples) {
  sample_name <- paste0("sample", i)
  gl_file <- file.path(gl_dir, paste0(sample_name, ".gl.tsv"))

  gl_data <- data.frame(
    chr = "chr1",
    pos = snp_positions,
    ref = ref_alleles,
    alt = alt_alleles,
    stringsAsFactors = FALSE
  )

  # Calculate GLs based on true genotype with some noise
  gl_rr <- numeric(n_snps)
  gl_ra <- numeric(n_snps)
  gl_aa <- numeric(n_snps)

  for (j in 1:n_snps) {
    gt <- genotypes[j, i]

    # Simulate realistic GLs with noise
    if (gt == "0/0") {
      gl_rr[j] <- log10(1 - error_rate)
      gl_ra[j] <- log10(error_rate / 2)
      gl_aa[j] <- log10(error_rate / 2)
    } else if (gt == "0/1") {
      gl_rr[j] <- log10(0.25 + runif(1, -0.1, 0.1))
      gl_ra[j] <- log10(0.5 + runif(1, -0.1, 0.1))
      gl_aa[j] <- log10(0.25 + runif(1, -0.1, 0.1))
    } else {  # 1/1
      gl_rr[j] <- log10(error_rate / 2)
      gl_ra[j] <- log10(error_rate / 2)
      gl_aa[j] <- log10(1 - error_rate)
    }
  }

  # Normalize GLs (max = 0)
  for (j in 1:n_snps) {
    max_gl <- max(gl_rr[j], gl_ra[j], gl_aa[j])
    gl_rr[j] <- gl_rr[j] - max_gl
    gl_ra[j] <- gl_ra[j] - max_gl
    gl_aa[j] <- gl_aa[j] - max_gl
  }

  gl_data$GL_RR <- round(gl_rr, 4)
  gl_data$GL_RA <- round(gl_ra, 4)
  gl_data$GL_AA <- round(gl_aa, 4)

  write.table(gl_data, gl_file, sep = "\t", row.names = FALSE, quote = FALSE)

  if (i %% 20 == 0) cat("  Written GL for", i, "/", opt$n_samples, "samples\n")
}

cat("\n=== Output Files ===\n")
cat("VCF:              ", vcf_file, "\n")
cat("Truth file:       ", truth_file, "\n")
cat("Sample summary:   ", summary_file, "\n")
cat("BED file:         ", bed_file, "\n")
cat("GL directory:     ", gl_dir, "/\n")
cat("\n")

# Print introgression size distribution
cat("=== Introgression Size Distribution ===\n")
sizes_mb <- introgression_truth$size / 1e6
cat("Min:", round(min(sizes_mb), 2), "Mb\n")
cat("Mean:", round(mean(sizes_mb), 2), "Mb\n")
cat("Max:", round(max(sizes_mb), 2), "Mb\n")
cat("Median:", round(median(sizes_mb), 2), "Mb\n")

# Print by sample type
cat("\n=== By Sample Type ===\n")
for (stype in c("mixed", "ref_alt", "ref_het")) {
  samples_of_type <- sample_summary$sample[sample_summary$type == stype]
  intros_of_type <- introgression_truth[introgression_truth$sample %in% samples_of_type, ]

  cat(stype, ":\n")
  cat("  N samples:", length(samples_of_type), "\n")
  cat("  N introgressions:", nrow(intros_of_type), "\n")
  if (nrow(intros_of_type) > 0) {
    cat("  Genotypes:", paste(names(table(intros_of_type$genotype)), "=", table(intros_of_type$genotype), collapse = ", "), "\n")
  }
}

cat("\nDone!\n")
