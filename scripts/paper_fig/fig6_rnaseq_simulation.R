#!/usr/bin/env Rscript
# Figure S1: RNA-seq Simulation Tool Comparison
# Supplementary figure comparing recommended RNA-seq simulators
#
# Panels:
#   A: Bulk RNA-seq tool comparison
#   B: scRNA-seq tool comparison
#   C: Feature support matrix
#   D: Example count distributions
#
# Usage: Rscript scripts/paper_fig/figS1_rnaseq_recommendations.R

library(ggplot2)
library(patchwork)
library(data.table)

# Output directory
out_dir <- "scripts/paper_fig/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Figure S1: RNA-seq Simulation Recommendations ===\n")

# -----------------------------------------------------------------------------
# Panel A: Bulk RNA-seq Tools
# -----------------------------------------------------------------------------

cat("Generating Panel A: Bulk RNA-seq Tools...\n")

bulk_tools <- data.frame(
  Tool = c("Polyester", "RSEM-sim", "Flux Simulator", "BEERS"),
  Language = c("R", "C++/Perl", "Java", "Python"),
  `Ease of Use` = c(5, 3, 2, 4),
  `Realism` = c(4, 5, 5, 4),
  `Speed` = c(4, 4, 3, 4),
  `Documentation` = c(5, 4, 3, 3),
  check.names = FALSE
)

# Reshape for radar-like bar chart
bulk_long <- reshape(bulk_tools, direction = "long",
                     varying = list(c("Ease of Use", "Realism", "Speed", "Documentation")),
                     v.names = "Score",
                     timevar = "Metric",
                     times = c("Ease of Use", "Realism", "Speed", "Documentation"))
bulk_long$Metric <- factor(bulk_long$Metric)

panel_a <- ggplot(bulk_long, aes(x = Tool, y = Score, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "A. Bulk RNA-seq Simulators",
    x = "",
    y = "Score (1-5)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  coord_cartesian(ylim = c(0, 5.5))

# -----------------------------------------------------------------------------
# Panel B: scRNA-seq Tools
# -----------------------------------------------------------------------------

cat("Generating Panel B: scRNA-seq Tools...\n")

scrna_tools <- data.frame(
  Tool = c("Splatter", "scDesign2", "SymSim", "SERGIO"),
  Language = c("R", "R", "R", "Python"),
  `Ease of Use` = c(5, 4, 3, 3),
  `Cell Types` = c(5, 5, 4, 5),
  `Batch Effects` = c(5, 4, 3, 2),
  `GRN Support` = c(2, 2, 3, 5),
  check.names = FALSE
)

scrna_long <- reshape(scrna_tools, direction = "long",
                      varying = list(c("Ease of Use", "Cell Types", "Batch Effects", "GRN Support")),
                      v.names = "Score",
                      timevar = "Metric",
                      times = c("Ease of Use", "Cell Types", "Batch Effects", "GRN Support"))
scrna_long$Metric <- factor(scrna_long$Metric)

panel_b <- ggplot(scrna_long, aes(x = Tool, y = Score, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "B. Single-Cell RNA-seq Simulators",
    x = "",
    y = "Score (1-5)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  coord_cartesian(ylim = c(0, 5.5))

# -----------------------------------------------------------------------------
# Panel C: Feature Support Matrix
# -----------------------------------------------------------------------------

cat("Generating Panel C: Feature Matrix...\n")

features <- expand.grid(
  Tool = c("Polyester", "RSEM", "Splatter", "scDesign2", "SERGIO"),
  Feature = c("Differential\nExpression", "Isoforms", "Batch\nEffects",
              "Cell\nTypes", "Dropouts", "UMIs", "Spatial")
)

# Define support
features$Supported <- c(
  # Polyester
  TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
  # RSEM
  TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
  # Splatter
  TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
  # scDesign2
  TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
  # SERGIO
  TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE
)

features$Tool <- factor(features$Tool,
  levels = c("Polyester", "RSEM", "Splatter", "scDesign2", "SERGIO"))

panel_c <- ggplot(features, aes(x = Feature, y = Tool, fill = Supported)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(Supported, "\u2713", "\u2717")),
            size = 4, color = ifelse(features$Supported, "darkgreen", "gray70")) +
  scale_fill_manual(values = c("gray95", "#90EE90"), guide = "none") +
  labs(
    title = "C. Feature Support Matrix",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

# -----------------------------------------------------------------------------
# Panel D: Simulated Count Distributions
# -----------------------------------------------------------------------------

cat("Generating Panel D: Count Distributions...\n")

set.seed(42)

# Simulate typical RNA-seq count distributions
n_genes <- 5000

# Bulk RNA-seq (negative binomial)
bulk_counts <- rnbinom(n_genes, mu = 100, size = 1)

# scRNA-seq with dropouts
sc_lambda <- rgamma(n_genes, shape = 0.5, rate = 0.01)
sc_counts <- rpois(n_genes, sc_lambda)
# Add dropout
dropout_prob <- exp(-0.1 * sc_lambda)
sc_counts[runif(n_genes) < dropout_prob] <- 0

count_data <- rbind(
  data.frame(type = "Bulk RNA-seq", counts = log10(bulk_counts + 1)),
  data.frame(type = "scRNA-seq", counts = log10(sc_counts + 1))
)

panel_d <- ggplot(count_data, aes(x = counts, fill = type)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "black",
                 linewidth = 0.2, position = "identity") +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  labs(
    title = "D. Simulated Count Distributions",
    x = expression(log[10](count + 1)),
    y = "Frequency"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "none",
    strip.background = element_rect(fill = "gray90")
  )

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

cat("\nCombining panels...\n")

figS1 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title = "Figure S1: Recommended RNA-seq Simulation Tools",
    subtitle = "External tools for bulk and single-cell RNA-seq simulation",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "figS1_rnaseq_recommendations.pdf"), figS1,
       width = 14, height = 12, dpi = 300)
ggsave(file.path(out_dir, "figS1_rnaseq_recommendations.png"), figS1,
       width = 14, height = 12, dpi = 300)

cat("\nFigure S1 saved to:\n")
cat("  ", file.path(out_dir, "figS1_rnaseq_recommendations.pdf"), "\n")
cat("  ", file.path(out_dir, "figS1_rnaseq_recommendations.png"), "\n")

# -----------------------------------------------------------------------------
# Save source data
# -----------------------------------------------------------------------------

write.csv(bulk_tools, file.path(out_dir, "figS1_source_data_bulk_tools.csv"), row.names = FALSE)
write.csv(scrna_tools, file.path(out_dir, "figS1_source_data_scrna_tools.csv"), row.names = FALSE)

cat("\nSource data saved.\n")
cat("\n=== Done ===\n")
