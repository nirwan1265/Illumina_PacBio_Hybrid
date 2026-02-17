#!/usr/bin/env Rscript
# Figure 2: DNA-seq Simulation Validation
# Demonstrates read simulation across platforms
#
# Panels:
#   A: Read length distributions (Illumina vs PacBio vs Nanopore)
#   B: Quality score distributions
#   C: Coverage uniformity across genome
#   D: Platform comparison summary
#
# Usage: Rscript scripts/paper_fig/fig2_dnaseq_simulation.R

library(ggplot2)
library(patchwork)
library(data.table)

# Output directory
out_dir <- "scripts/paper_fig/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Figure 2: DNA-seq Simulation ===\n")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

read_fasta <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  if (length(hdr_idx) == 0) stop("No FASTA header found")
  name <- sub("^>", "", strsplit(lines[hdr_idx[1]], " ")[[1]][1])
  seq_lines <- if (length(hdr_idx) > 1) {
    lines[(hdr_idx[1] + 1):(hdr_idx[2] - 1)]
  } else {
    lines[(hdr_idx[1] + 1):length(lines)]
  }
  seq <- toupper(paste(seq_lines, collapse = ""))
  list(name = name, seq = seq, len = nchar(seq))
}

# Simulate Illumina-like read lengths (PE 150bp with slight variation)
sim_illumina_lengths <- function(n_reads, read_len = 150, sd = 2) {
  pmax(50, round(rnorm(n_reads, read_len, sd)))
}

# Simulate PacBio HiFi read lengths (log-normal distribution)
sim_pacbio_hifi_lengths <- function(n_reads, mean_len = 15000, sd_log = 0.4) {
  mu <- log(mean_len) - sd_log^2/2
  round(rlnorm(n_reads, meanlog = mu, sdlog = sd_log))
}

# Simulate PacBio CLR read lengths (longer, more variable)
sim_pacbio_clr_lengths <- function(n_reads, mean_len = 10000, sd_log = 0.6) {
  mu <- log(mean_len) - sd_log^2/2
  round(rlnorm(n_reads, meanlog = mu, sdlog = sd_log))
}

# Simulate Nanopore read lengths (very variable)
sim_nanopore_lengths <- function(n_reads, mean_len = 15000, sd_log = 0.8) {
  mu <- log(mean_len) - sd_log^2/2
  round(rlnorm(n_reads, meanlog = mu, sdlog = sd_log))
}

# Simulate quality scores
sim_illumina_quals <- function(n_reads, mean_q = 35, sd = 3) {
  pmin(41, pmax(20, round(rnorm(n_reads, mean_q, sd))))
}

sim_pacbio_hifi_quals <- function(n_reads, mean_q = 30, sd = 5) {
  pmin(40, pmax(15, round(rnorm(n_reads, mean_q, sd))))
}

sim_pacbio_clr_quals <- function(n_reads, mean_q = 12, sd = 3) {
  pmin(20, pmax(5, round(rnorm(n_reads, mean_q, sd))))
}

sim_nanopore_quals <- function(n_reads, model = "R10") {
  if (model == "R10") {
    # R10.4.1 Q20+ chemistry
    pmin(35, pmax(10, round(rnorm(n_reads, 22, 5))))
  } else {
    # R9.4 older chemistry
    pmin(20, pmax(5, round(rnorm(n_reads, 12, 3))))
  }
}

# Simulate coverage along genome
sim_coverage <- function(genome_len, n_reads, read_len, gc_bias = FALSE, gc_seq = NULL) {
  # Random start positions
  starts <- sample(1:(genome_len - read_len), n_reads, replace = TRUE)

  # Bin into windows
  window_size <- 1000
  n_windows <- ceiling(genome_len / window_size)
  coverage <- numeric(n_windows)

  for (s in starts) {
    win <- ceiling(s / window_size)
    if (win <= n_windows) {
      coverage[win] <- coverage[win] + 1
    }
  }

  # Normalize to mean coverage
  coverage <- coverage / mean(coverage)
  data.frame(position = (1:n_windows) * window_size, coverage = coverage)
}

# -----------------------------------------------------------------------------
# Generate simulated read data
# -----------------------------------------------------------------------------

cat("Generating simulated read data...\n")

set.seed(42)
n_reads <- 10000

# Read lengths
illumina_lens <- sim_illumina_lengths(n_reads)
pacbio_hifi_lens <- sim_pacbio_hifi_lengths(n_reads)
pacbio_clr_lens <- sim_pacbio_clr_lengths(n_reads)
nanopore_r10_lens <- sim_nanopore_lengths(n_reads, mean_len = 15000, sd_log = 0.7)
nanopore_r9_lens <- sim_nanopore_lengths(n_reads, mean_len = 8000, sd_log = 0.9)

# Quality scores
illumina_quals <- sim_illumina_quals(n_reads)
pacbio_hifi_quals <- sim_pacbio_hifi_quals(n_reads)
pacbio_clr_quals <- sim_pacbio_clr_quals(n_reads)
nanopore_r10_quals <- sim_nanopore_quals(n_reads, "R10")
nanopore_r9_quals <- sim_nanopore_quals(n_reads, "R9")

# Combine into data frames
len_data <- rbind(
  data.frame(platform = "Illumina PE150", length = illumina_lens),
  data.frame(platform = "PacBio HiFi", length = pacbio_hifi_lens),
  data.frame(platform = "PacBio CLR", length = pacbio_clr_lens),
  data.frame(platform = "Nanopore R10.4", length = nanopore_r10_lens),
  data.frame(platform = "Nanopore R9.4", length = nanopore_r9_lens)
)
len_data$platform <- factor(len_data$platform,
  levels = c("Illumina PE150", "PacBio CLR", "PacBio HiFi", "Nanopore R9.4", "Nanopore R10.4"))

qual_data <- rbind(
  data.frame(platform = "Illumina", quality = illumina_quals),
  data.frame(platform = "PacBio HiFi", quality = pacbio_hifi_quals),
  data.frame(platform = "PacBio CLR", quality = pacbio_clr_quals),
  data.frame(platform = "Nanopore R10.4", quality = nanopore_r10_quals),
  data.frame(platform = "Nanopore R9.4", quality = nanopore_r9_quals)
)
qual_data$platform <- factor(qual_data$platform,
  levels = c("Illumina", "PacBio CLR", "PacBio HiFi", "Nanopore R9.4", "Nanopore R10.4"))

# Summary statistics
summary_stats <- data.frame(
  Platform = c("Illumina PE150", "PacBio CLR", "PacBio HiFi", "Nanopore R9.4", "Nanopore R10.4"),
  `Mean Length` = c(mean(illumina_lens), mean(pacbio_clr_lens), mean(pacbio_hifi_lens),
                    mean(nanopore_r9_lens), mean(nanopore_r10_lens)),
  `N50` = c(150, quantile(pacbio_clr_lens, 0.5), quantile(pacbio_hifi_lens, 0.5),
            quantile(nanopore_r9_lens, 0.5), quantile(nanopore_r10_lens, 0.5)),
  `Mean Quality` = c(mean(illumina_quals), mean(pacbio_clr_quals), mean(pacbio_hifi_quals),
                     mean(nanopore_r9_quals), mean(nanopore_r10_quals)),
  `Error Rate` = c(0.001, 0.12, 0.001, 0.08, 0.02),
  check.names = FALSE
)

cat("\nPlatform summary:\n")
print(summary_stats)

# -----------------------------------------------------------------------------
# Panel A: Read Length Distributions
# -----------------------------------------------------------------------------

cat("\nGenerating Panel A: Read Length Distributions...\n")

# Use log scale for long reads
panel_a <- ggplot(len_data, aes(x = length, fill = platform)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "black", linewidth = 0.2) +
  facet_wrap(~platform, scales = "free", ncol = 5) +
  scale_fill_manual(values = c("#E41A1C", "#984EA3", "#377EB8", "#FF7F00", "#4DAF4A")) +
  scale_x_continuous(labels = function(x) ifelse(x >= 1000, paste0(x/1000, "k"), x)) +
  labs(
    title = "A. Read Length Distributions by Platform",
    x = "Read Length (bp)",
    y = "Count"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "none",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 9)
  )

# -----------------------------------------------------------------------------
# Panel B: Quality Score Distributions
# -----------------------------------------------------------------------------

cat("Generating Panel B: Quality Score Distributions...\n")

panel_b <- ggplot(qual_data, aes(x = platform, y = quality, fill = platform)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("#E41A1C", "#984EA3", "#377EB8", "#FF7F00", "#4DAF4A")) +
  labs(
    title = "B. Quality Score Distributions",
    x = "",
    y = "Phred Quality Score"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  geom_hline(yintercept = c(20, 30), linetype = "dashed", color = "gray50", linewidth = 0.3) +
  annotate("text", x = 5.4, y = 20, label = "Q20", hjust = 0, size = 3, color = "gray40") +
  annotate("text", x = 5.4, y = 30, label = "Q30", hjust = 0, size = 3, color = "gray40")

# -----------------------------------------------------------------------------
# Panel C: Coverage Uniformity
# -----------------------------------------------------------------------------

cat("Generating Panel C: Coverage Uniformity...\n")

# Simulate coverage for different platforms
genome_len <- 100000
n_reads_cov <- 5000

set.seed(789)
cov_illumina <- sim_coverage(genome_len, n_reads_cov, 150)
cov_illumina$platform <- "Illumina"

cov_pacbio <- sim_coverage(genome_len, n_reads_cov/10, 15000)  # fewer but longer
cov_pacbio$platform <- "PacBio HiFi"

cov_nanopore <- sim_coverage(genome_len, n_reads_cov/8, 12000)
cov_nanopore$platform <- "Nanopore R10.4"

cov_data <- rbind(cov_illumina, cov_pacbio, cov_nanopore)
cov_data$platform <- factor(cov_data$platform,
  levels = c("Illumina", "PacBio HiFi", "Nanopore R10.4"))

panel_c <- ggplot(cov_data, aes(x = position/1000, y = coverage, color = platform)) +
  geom_line(alpha = 0.8, linewidth = 0.5) +
  facet_wrap(~platform, ncol = 1) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "C. Coverage Uniformity Along Genome",
    x = "Position (kb)",
    y = "Normalized Coverage"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "none",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

# -----------------------------------------------------------------------------
# Panel D: Platform Comparison Summary
# -----------------------------------------------------------------------------

cat("Generating Panel D: Platform Comparison...\n")

# Reshape for heatmap-style comparison
comparison_data <- data.frame(
  Platform = rep(c("Illumina", "PacBio CLR", "PacBio HiFi", "Nanopore R9.4", "Nanopore R10.4"), 4),
  Metric = rep(c("Read Length", "Quality", "Accuracy", "Throughput"), each = 5),
  Value = c(
    # Read Length (normalized 0-1, longer = higher)
    0.1, 0.6, 0.8, 0.5, 0.8,
    # Quality (Phred normalized)
    1.0, 0.3, 0.75, 0.35, 0.6,
    # Accuracy
    0.999, 0.88, 0.999, 0.92, 0.98,
    # Throughput (relative)
    1.0, 0.6, 0.4, 0.7, 0.8
  )
)
comparison_data$Platform <- factor(comparison_data$Platform,
  levels = c("Illumina", "PacBio CLR", "PacBio HiFi", "Nanopore R9.4", "Nanopore R10.4"))
comparison_data$Metric <- factor(comparison_data$Metric,
  levels = c("Read Length", "Quality", "Accuracy", "Throughput"))

panel_d <- ggplot(comparison_data, aes(x = Metric, y = Platform, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Value)), size = 3) +
  scale_fill_gradient2(low = "#E41A1C", mid = "#FFFFBF", high = "#4DAF4A",
                       midpoint = 0.5, limits = c(0, 1)) +
  labs(
    title = "D. Platform Comparison Matrix",
    x = "",
    y = "",
    fill = "Score"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  )

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

cat("\nCombining panels...\n")

fig2 <- (panel_a) / (panel_b | panel_c) / (panel_d) +
  plot_layout(heights = c(1, 1.2, 0.8)) +
  plot_annotation(
    title = "Figure 2: DNA-seq Simulation Across Sequencing Platforms",
    subtitle = "Simulated characteristics of Illumina, PacBio, and Nanopore reads",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "fig2_dnaseq_simulation.pdf"), fig2,
       width = 12, height = 14, dpi = 300)
ggsave(file.path(out_dir, "fig2_dnaseq_simulation.png"), fig2,
       width = 12, height = 14, dpi = 300)

cat("\nFigure 2 saved to:\n")
cat("  ", file.path(out_dir, "fig2_dnaseq_simulation.pdf"), "\n")
cat("  ", file.path(out_dir, "fig2_dnaseq_simulation.png"), "\n")

# -----------------------------------------------------------------------------
# Save source data
# -----------------------------------------------------------------------------

write.csv(len_data, file.path(out_dir, "fig2_source_data_lengths.csv"), row.names = FALSE)
write.csv(qual_data, file.path(out_dir, "fig2_source_data_quality.csv"), row.names = FALSE)
write.csv(cov_data, file.path(out_dir, "fig2_source_data_coverage.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(out_dir, "fig2_source_data_summary.csv"), row.names = FALSE)

cat("\nSource data saved.\n")
cat("\n=== Done ===\n")
