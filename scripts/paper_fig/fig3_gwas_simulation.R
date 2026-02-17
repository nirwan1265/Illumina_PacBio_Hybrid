#!/usr/bin/env Rscript
# Figure 3: GWAS Cohort Simulation Validation
# Demonstrates population genetics simulation capabilities
#
# Panels:
#   A: Allele frequency spectrum
#   B: LD decay curve
#   C: Population structure (PCA)
#   D: Phenotype vs genotype correlation
#
# Usage: Rscript scripts/paper_fig/fig3_gwas_simulation.R

library(ggplot2)
library(patchwork)
library(data.table)

# Output directory
out_dir <- "scripts/paper_fig/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Figure 3: GWAS Cohort Simulation ===\n")

# -----------------------------------------------------------------------------
# Simulate GWAS cohort
# -----------------------------------------------------------------------------

set.seed(42)

# Parameters
n_samples <- 300
n_variants <- 5000
n_pops <- 3
fst <- 0.05
n_causal <- 50
effect_sd <- 0.5

cat("Simulating GWAS cohort...\n")
cat("  Samples:", n_samples, "\n")
cat("  Variants:", n_variants, "\n")
cat("  Populations:", n_pops, "\n")
cat("  Fst:", fst, "\n")

# Population assignments
pop_sizes <- c(100, 100, 100)
pop_ids <- rep(paste0("Pop", 1:n_pops), times = pop_sizes)
sample_ids <- paste0("sample_", 1:n_samples)

# Generate allele frequencies per population (with Fst drift)
# Ancestral frequencies from Beta distribution
p_ancestral <- rbeta(n_variants, 0.8, 0.8)

# Population-specific frequencies with drift
p_pop <- matrix(0, nrow = n_variants, ncol = n_pops)
for (k in 1:n_pops) {
  alpha <- p_ancestral * (1 - fst) / fst
  beta <- (1 - p_ancestral) * (1 - fst) / fst
  # Bound to avoid 0/Inf
  alpha <- pmax(0.01, alpha)
  beta <- pmax(0.01, beta)
  p_pop[, k] <- rbeta(n_variants, alpha, beta)
  p_pop[, k] <- pmin(0.99, pmax(0.01, p_pop[, k]))
}

# Generate genotypes (0, 1, 2)
pop_index <- match(pop_ids, paste0("Pop", 1:n_pops))
G <- matrix(0, nrow = n_variants, ncol = n_samples)
for (i in 1:n_variants) {
  for (s in 1:n_samples) {
    G[i, s] <- rbinom(1, 2, p_pop[i, pop_index[s]])
  }
}

# Variant positions (for LD calculation)
positions <- sort(sample(1:10000000, n_variants))

# Generate phenotype
causal_idx <- sample(1:n_variants, n_causal)
effects <- rnorm(n_causal, 0, effect_sd)
genetic_score <- as.numeric(t(effects) %*% G[causal_idx, ])
phenotype <- genetic_score + rnorm(n_samples, 0, 1)

# Add population effect
pop_effect <- (pop_index - 1) * 0.3
phenotype <- phenotype + pop_effect

cat("Genotype matrix:", nrow(G), "x", ncol(G), "\n")
cat("Phenotype range:", round(min(phenotype), 2), "to", round(max(phenotype), 2), "\n")

# -----------------------------------------------------------------------------
# Panel A: Allele Frequency Spectrum
# -----------------------------------------------------------------------------

cat("\nGenerating Panel A: Allele Frequency Spectrum...\n")

# Calculate minor allele frequency
maf <- apply(G, 1, function(x) {
  af <- mean(x) / 2
  min(af, 1 - af)
})

af_data <- data.frame(maf = maf)

# Expected SFS under neutral model (Watterson)
expected_sfs <- function(x, n) {
  # Folded SFS approximation
  1 / (x * (1-x))
}

panel_a <- ggplot(af_data, aes(x = maf)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#377EB8", color = "black", alpha = 0.7, linewidth = 0.3) +
  stat_function(fun = function(x) dbeta(x, 0.8, 0.8) * 2,
                color = "#E41A1C", linewidth = 1, linetype = "dashed") +
  labs(
    title = "A. Minor Allele Frequency Spectrum",
    x = "Minor Allele Frequency",
    y = "Density"
  ) +
  annotate("text", x = 0.4, y = 3, label = "Simulated", color = "#377EB8", size = 4) +
  annotate("text", x = 0.4, y = 2.5, label = "Expected (Beta)", color = "#E41A1C", size = 4) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  xlim(0, 0.5)

# -----------------------------------------------------------------------------
# Panel B: LD Decay
# -----------------------------------------------------------------------------

cat("Generating Panel B: LD Decay...\n")

# Calculate r^2 for pairs of SNPs at different distances
calc_r2 <- function(g1, g2) {
  cor(g1, g2)^2
}

# Sample pairs at different distances
distances <- c(100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000)
ld_results <- list()

for (d in distances) {
  # Find pairs approximately at this distance
  pairs_found <- 0
  r2_vals <- c()

  for (attempt in 1:500) {
    i <- sample(1:(n_variants-1), 1)
    j_candidates <- which(abs(positions - positions[i] - d) < d * 0.2)
    j_candidates <- j_candidates[j_candidates > i]

    if (length(j_candidates) > 0) {
      j <- j_candidates[1]
      r2 <- calc_r2(G[i, ], G[j, ])
      if (!is.na(r2)) {
        r2_vals <- c(r2_vals, r2)
        pairs_found <- pairs_found + 1
      }
    }
    if (pairs_found >= 50) break
  }

  if (length(r2_vals) > 0) {
    ld_results[[length(ld_results) + 1]] <- data.frame(
      distance = d,
      r2_mean = mean(r2_vals, na.rm = TRUE),
      r2_sd = sd(r2_vals, na.rm = TRUE),
      n_pairs = length(r2_vals)
    )
  }
}

ld_data <- do.call(rbind, ld_results)
ld_data$distance_kb <- ld_data$distance / 1000

panel_b <- ggplot(ld_data, aes(x = distance_kb, y = r2_mean)) +
  geom_ribbon(aes(ymin = pmax(0, r2_mean - r2_sd),
                  ymax = pmin(1, r2_mean + r2_sd)),
              fill = "#377EB8", alpha = 0.3) +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_point(color = "#377EB8", size = 3) +
  scale_x_log10(labels = function(x) paste0(x, " kb")) +
  labs(
    title = "B. Linkage Disequilibrium Decay",
    x = "Distance (log scale)",
    y = expression(Mean~r^2)
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  ylim(0, 1)

# -----------------------------------------------------------------------------
# Panel C: Population Structure (PCA)
# -----------------------------------------------------------------------------

cat("Generating Panel C: Population Structure...\n")

# Simple PCA on genotype matrix
G_centered <- t(G) - rowMeans(t(G))
pca <- prcomp(G_centered, scale. = FALSE)

pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Population = pop_ids
)

# Variance explained
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

panel_c <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(
    title = "C. Population Structure (PCA)",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

# -----------------------------------------------------------------------------
# Panel D: Genotype-Phenotype Association
# -----------------------------------------------------------------------------

cat("Generating Panel D: Genotype-Phenotype Association...\n")

# Calculate association for all variants
p_values <- sapply(1:n_variants, function(i) {
  fit <- lm(phenotype ~ G[i, ])
  summary(fit)$coefficients[2, 4]
})

# Manhattan-style data
assoc_data <- data.frame(
  position = positions,
  pvalue = p_values,
  neglog10p = -log10(p_values),
  is_causal = 1:n_variants %in% causal_idx
)

# Significance threshold
bonf_threshold <- -log10(0.05 / n_variants)

panel_d <- ggplot(assoc_data, aes(x = position/1e6, y = neglog10p)) +
  geom_point(aes(color = is_causal), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = bonf_threshold, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_color_manual(values = c("gray50", "#E41A1C"),
                     labels = c("Non-causal", "Causal"),
                     name = "") +
  labs(
    title = "D. Genotype-Phenotype Association",
    x = "Position (Mb)",
    y = expression(-log[10](p))
  ) +
  annotate("text", x = 9.5, y = bonf_threshold + 0.5,
           label = "Bonferroni", color = "red", size = 3) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = c(0.9, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

cat("\nCombining panels...\n")

fig3 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title = "Figure 3: GWAS Cohort Simulation",
    subtitle = sprintf("%d samples, %d variants, %d populations (Fst=%.2f), %d causal variants",
                       n_samples, n_variants, n_pops, fst, n_causal),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "fig3_gwas_simulation.pdf"), fig3,
       width = 12, height = 10, dpi = 300)
ggsave(file.path(out_dir, "fig3_gwas_simulation.png"), fig3,
       width = 12, height = 10, dpi = 300)

cat("\nFigure 3 saved to:\n")
cat("  ", file.path(out_dir, "fig3_gwas_simulation.pdf"), "\n")
cat("  ", file.path(out_dir, "fig3_gwas_simulation.png"), "\n")

# -----------------------------------------------------------------------------
# Save source data
# -----------------------------------------------------------------------------

write.csv(af_data, file.path(out_dir, "fig3_source_data_maf.csv"), row.names = FALSE)
write.csv(ld_data, file.path(out_dir, "fig3_source_data_ld.csv"), row.names = FALSE)
write.csv(pca_data, file.path(out_dir, "fig3_source_data_pca.csv"), row.names = FALSE)
write.csv(assoc_data, file.path(out_dir, "fig3_source_data_assoc.csv"), row.names = FALSE)

cat("\nSource data saved.\n")
cat("\n=== Done ===\n")
