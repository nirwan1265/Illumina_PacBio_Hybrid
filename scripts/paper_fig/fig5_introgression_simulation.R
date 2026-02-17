#!/usr/bin/env Rscript
# Figure 5: Cross-Species Introgression Simulation
# Demonstrates BC2S3-like introgression breeding simulation
#
# Panels:
#   A: Introgression landscape across samples (heatmap)
#   B: Donor fraction distribution by sample type
#   C: Introgression segment size distribution
#   D: Genotype accuracy at different coverages
#
# Usage: Rscript scripts/paper_fig/fig5_introgression_simulation.R

library(ggplot2)
library(patchwork)
library(data.table)

# Output directory
out_dir <- "scripts/paper_fig/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Figure 5: Cross-Species Introgression Simulation ===\n")

# -----------------------------------------------------------------------------
# Simulate BC2S3-like introgression population
# -----------------------------------------------------------------------------

set.seed(42)

# Parameters (matching 22_simulate_controlled_introgressions.R)
genome_size <- 1000000  # 1 Mb for speed
n_samples <- 100
n_snps <- 10000
snp_positions <- sort(sample(1:genome_size, n_snps))

# Sample types
n_mixed <- 50      # REF, HET, ALT
n_ref_alt <- 25    # Homozygous only
n_ref_het <- 25    # REF and HET only

cat("Simulating introgression population...\n")
cat("  Genome size:", genome_size, "bp\n")
cat("  Samples:", n_samples, "\n")
cat("  SNPs:", n_snps, "\n")

# Generate introgression segments for each sample
generate_introgressions <- function(genome_size, target_frac = 0.125,
                                    min_len = 10000, max_len = 100000) {
  segments <- list()
  total_introg <- 0
  target_bp <- genome_size * target_frac

  while (total_introg < target_bp * 0.8) {
    start <- sample(1:(genome_size - min_len), 1)
    len <- runif(1, min_len, max_len)
    end <- min(start + len, genome_size)

    segments[[length(segments) + 1]] <- c(start, end)
    total_introg <- total_introg + (end - start)

    if (length(segments) > 20) break  # Safety limit
  }

  do.call(rbind, segments)
}

# Generate genotypes based on introgression segments
generate_genotypes <- function(snp_positions, introgressions, sample_type) {
  n_snps <- length(snp_positions)
  genotypes <- rep(0, n_snps)  # Default: REF/REF

  if (is.null(introgressions) || nrow(introgressions) == 0) {
    return(genotypes)
  }

  # Mark SNPs in introgression regions
  for (i in 1:nrow(introgressions)) {
    in_region <- snp_positions >= introgressions[i, 1] &
                 snp_positions <= introgressions[i, 2]

    if (sample_type == "mixed") {
      # Random genotypes: 0, 1, or 2
      genotypes[in_region] <- sample(0:2, sum(in_region), replace = TRUE,
                                     prob = c(0.25, 0.5, 0.25))
    } else if (sample_type == "ref_alt") {
      # Homozygous only: 0 or 2
      genotypes[in_region] <- sample(c(0, 2), sum(in_region), replace = TRUE)
    } else if (sample_type == "ref_het") {
      # REF or HET only: 0 or 1
      genotypes[in_region] <- sample(0:1, sum(in_region), replace = TRUE)
    }
  }

  genotypes
}

# Generate all samples
samples <- list()
sample_info <- data.frame(
  sample_id = character(),
  type = character(),
  donor_fraction = numeric(),
  n_introgressions = integer(),
  stringsAsFactors = FALSE
)

for (i in 1:n_samples) {
  if (i <= n_mixed) {
    type <- "mixed"
  } else if (i <= n_mixed + n_ref_alt) {
    type <- "ref_alt"
  } else {
    type <- "ref_het"
  }

  # Generate introgressions
  introg <- generate_introgressions(genome_size,
                                    target_frac = runif(1, 0.08, 0.18))

  # Generate genotypes
  geno <- generate_genotypes(snp_positions, introg, type)

  samples[[i]] <- list(
    genotypes = geno,
    introgressions = introg,
    type = type
  )

  # Calculate donor fraction
  donor_frac <- mean(geno > 0)

  sample_info <- rbind(sample_info, data.frame(
    sample_id = paste0("sample", i),
    type = type,
    donor_fraction = donor_frac,
    n_introgressions = nrow(introg),
    stringsAsFactors = FALSE
  ))
}

cat("Sample types:\n")
print(table(sample_info$type))

# -----------------------------------------------------------------------------
# Panel A: Introgression Landscape (Heatmap)
# -----------------------------------------------------------------------------

cat("\nGenerating Panel A: Introgression Landscape...\n")

# Create genotype matrix for heatmap
geno_matrix <- do.call(rbind, lapply(samples, function(s) s$genotypes))
rownames(geno_matrix) <- sample_info$sample_id

# Bin SNPs for visualization
n_bins <- 100
bin_size <- ceiling(n_snps / n_bins)
binned_geno <- matrix(0, nrow = n_samples, ncol = n_bins)

for (b in 1:n_bins) {
  start_idx <- (b - 1) * bin_size + 1
  end_idx <- min(b * bin_size, n_snps)
  binned_geno[, b] <- rowMeans(geno_matrix[, start_idx:end_idx, drop = FALSE])
}

# Convert to long format
heatmap_data <- expand.grid(
  sample = 1:n_samples,
  bin = 1:n_bins
)
heatmap_data$value <- as.vector(binned_geno)
heatmap_data$type <- sample_info$type[heatmap_data$sample]

# Order samples by type and donor fraction
sample_order <- order(sample_info$type, sample_info$donor_fraction)
heatmap_data$sample_ordered <- match(heatmap_data$sample, sample_order)

panel_a <- ggplot(heatmap_data, aes(x = bin, y = sample_ordered, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 1, name = "Genotype\n(0=REF, 2=ALT)") +
  scale_x_continuous(expand = c(0, 0),
                     labels = function(x) paste0(round(x/n_bins * genome_size/1e6, 1), " Mb")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "A. Introgression Landscape Across Samples",
    x = "Genomic Position",
    y = "Sample (ordered by type)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_blank(),
    legend.position = "right"
  ) +
  # Add type annotations
  annotate("text", x = -3, y = 25, label = "Mixed", angle = 90, size = 3, hjust = 0.5) +
  annotate("text", x = -3, y = 65, label = "Ref/Alt", angle = 90, size = 3, hjust = 0.5) +
  annotate("text", x = -3, y = 90, label = "Ref/Het", angle = 90, size = 3, hjust = 0.5)

# -----------------------------------------------------------------------------
# Panel B: Donor Fraction Distribution
# -----------------------------------------------------------------------------

cat("Generating Panel B: Donor Fraction Distribution...\n")

sample_info$type_label <- factor(sample_info$type,
  levels = c("mixed", "ref_alt", "ref_het"),
  labels = c("Mixed\n(REF/HET/ALT)", "Ref-Alt\n(Homozygous)", "Ref-Het\n(Heterozygous)"))

panel_b <- ggplot(sample_info, aes(x = type_label, y = donor_fraction, fill = type_label)) +
  geom_violin(alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 1) +
  geom_hline(yintercept = 0.125, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  annotate("text", x = 3.4, y = 0.125, label = "Expected\n(12.5%)",
           color = "red", size = 3, hjust = 0) +
  labs(
    title = "B. Donor Fraction by Sample Type",
    x = "",
    y = "Donor Allele Fraction"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(0.5, 3.8))

# -----------------------------------------------------------------------------
# Panel C: Introgression Segment Size Distribution
# -----------------------------------------------------------------------------

cat("Generating Panel C: Segment Size Distribution...\n")

# Collect all introgression segments
all_segments <- list()
for (i in 1:n_samples) {
  introg <- samples[[i]]$introgressions
  if (!is.null(introg) && nrow(introg) > 0) {
    for (j in 1:nrow(introg)) {
      all_segments[[length(all_segments) + 1]] <- data.frame(
        sample = i,
        type = samples[[i]]$type,
        start = introg[j, 1],
        end = introg[j, 2],
        length = introg[j, 2] - introg[j, 1]
      )
    }
  }
}
segment_df <- do.call(rbind, all_segments)

panel_c <- ggplot(segment_df, aes(x = length / 1000)) +
  geom_histogram(bins = 30, fill = "#FF7F00", color = "black",
                 alpha = 0.7, linewidth = 0.3) +
  geom_vline(xintercept = mean(segment_df$length) / 1000,
             linetype = "dashed", color = "#E41A1C", linewidth = 1) +
  annotate("text", x = mean(segment_df$length) / 1000 + 5, y = Inf,
           label = sprintf("Mean: %.1f kb", mean(segment_df$length) / 1000),
           vjust = 2, hjust = 0, color = "#E41A1C", size = 3.5) +
  labs(
    title = "C. Introgression Segment Size Distribution",
    x = "Segment Length (kb)",
    y = "Count"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12))

# -----------------------------------------------------------------------------
# Panel D: Genotype Accuracy vs Coverage
# -----------------------------------------------------------------------------

cat("Generating Panel D: Genotype Accuracy vs Coverage...\n")

# Simulate genotype calling accuracy at different coverages
# Based on binomial sampling of reads
simulate_gl_accuracy <- function(true_geno, coverage, error_rate = 0.01) {
  n_snps <- length(true_geno)
  called_geno <- numeric(n_snps)

  for (i in 1:n_snps) {
    # Expected allele frequency based on true genotype
    true_af <- true_geno[i] / 2

    # Sample reads
    n_reads <- rpois(1, coverage)
    if (n_reads == 0) {
      called_geno[i] <- NA
      next
    }

    # Alt reads with error
    n_alt <- rbinom(1, n_reads, true_af * (1 - error_rate) + (1 - true_af) * error_rate)
    obs_af <- n_alt / n_reads

    # Call genotype
    if (obs_af < 0.2) {
      called_geno[i] <- 0
    } else if (obs_af > 0.8) {
      called_geno[i] <- 2
    } else {
      called_geno[i] <- 1
    }
  }

  # Calculate accuracy
  valid <- !is.na(called_geno)
  accuracy <- mean(called_geno[valid] == true_geno[valid])
  het_accuracy <- NA
  if (sum(true_geno == 1 & valid) > 0) {
    het_accuracy <- mean(called_geno[true_geno == 1 & valid] == 1)
  }

  c(accuracy = accuracy, het_accuracy = het_accuracy, missing = mean(!valid))
}

# Test at multiple coverages
coverages <- c(0.5, 1, 2, 4, 8, 10, 15, 20, 30)
accuracy_results <- list()

for (cov in coverages) {
  cat("  Testing coverage:", cov, "x\n")

  # Use first 10 mixed samples
  for (s in 1:10) {
    true_geno <- samples[[s]]$genotypes
    acc <- simulate_gl_accuracy(true_geno, cov)

    accuracy_results[[length(accuracy_results) + 1]] <- data.frame(
      coverage = cov,
      sample = s,
      accuracy = acc["accuracy"],
      het_accuracy = acc["het_accuracy"],
      missing = acc["missing"]
    )
  }
}

accuracy_df <- do.call(rbind, accuracy_results)

# Summarize
accuracy_summary <- aggregate(cbind(accuracy, het_accuracy) ~ coverage,
                              data = accuracy_df,
                              FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                  sd = sd(x, na.rm = TRUE)))
accuracy_summary <- do.call(data.frame, accuracy_summary)
names(accuracy_summary) <- c("coverage", "accuracy_mean", "accuracy_sd",
                             "het_accuracy_mean", "het_accuracy_sd")

# Long format for plotting
acc_long <- rbind(
  data.frame(coverage = accuracy_summary$coverage,
             metric = "Overall",
             mean = accuracy_summary$accuracy_mean,
             sd = accuracy_summary$accuracy_sd),
  data.frame(coverage = accuracy_summary$coverage,
             metric = "Heterozygote",
             mean = accuracy_summary$het_accuracy_mean,
             sd = accuracy_summary$het_accuracy_sd)
)

panel_d <- ggplot(acc_long, aes(x = coverage, y = mean * 100, color = metric)) +
  geom_ribbon(aes(ymin = (mean - sd) * 100, ymax = (mean + sd) * 100, fill = metric),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#377EB8", "#E41A1C"), name = "Accuracy") +
  scale_fill_manual(values = c("#377EB8", "#E41A1C"), name = "Accuracy") +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray50") +
  annotate("text", x = 25, y = 96, label = "95% threshold", color = "gray50", size = 3) +
  labs(
    title = "D. Genotype Calling Accuracy vs Coverage",
    x = "Sequencing Coverage",
    y = "Accuracy (%)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = c(0.8, 0.3),
    legend.background = element_rect(fill = "white", color = "gray80")
  ) +
  coord_cartesian(ylim = c(50, 100))

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

cat("\nCombining panels...\n")

fig5 <- (panel_a) / (panel_b | panel_c) / (panel_d) +
  plot_layout(heights = c(1.2, 1, 1)) +
  plot_annotation(
    title = "Figure 5: Cross-Species Introgression Simulation",
    subtitle = sprintf("BC2S3-like population: %d samples, %.0f kb genome, %d SNPs",
                       n_samples, genome_size/1000, n_snps),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "fig5_introgression_simulation.pdf"), fig5,
       width = 12, height = 14, dpi = 300)
ggsave(file.path(out_dir, "fig5_introgression_simulation.png"), fig5,
       width = 12, height = 14, dpi = 300)

cat("\nFigure 5 saved to:\n")
cat("  ", file.path(out_dir, "fig5_introgression_simulation.pdf"), "\n")
cat("  ", file.path(out_dir, "fig5_introgression_simulation.png"), "\n")

# -----------------------------------------------------------------------------
# Save source data
# -----------------------------------------------------------------------------

write.csv(sample_info, file.path(out_dir, "fig5_source_data_samples.csv"), row.names = FALSE)
write.csv(segment_df, file.path(out_dir, "fig5_source_data_segments.csv"), row.names = FALSE)
write.csv(accuracy_df, file.path(out_dir, "fig5_source_data_accuracy.csv"), row.names = FALSE)

cat("\nSource data saved.\n")
cat("\n=== Done ===\n")
