#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Figure S2: Advanced Phenotype Simulation Validation
# =============================================================================
#
# Validates phenotype simulation capabilities including:
#   A: Heritability accuracy (simulated vs target h²)
#   B: Genetic architecture comparison (A, AD, ADE models)
#   C: Multi-trait genetic correlations
#   D: Binary trait liability threshold model
#   E: Gene-environment interaction variance decomposition
#   F: Covariate effects validation
#
# Usage: Rscript scripts/paper_fig/figS2_phenotype_simulation.R
# =============================================================================

library(ggplot2)
library(patchwork)
library(data.table)

# Output directory
out_dir <- "scripts/paper_fig/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Supplementary Figure S2: Phenotype Simulation Validation ===\n")

set.seed(42)

# =============================================================================
# Simulation Parameters
# =============================================================================

n_samples <- 500
n_snps <- 1000
n_qtn <- 20  # Causal variants

cat("Simulating genotypes...\n")
cat("  Samples:", n_samples, "\n")
cat("  SNPs:", n_snps, "\n")
cat("  QTNs:", n_qtn, "\n\n")

# Generate random genotype matrix (0, 1, 2)
# MAF drawn from Beta(0.8, 0.8) - U-shaped distribution
maf <- rbeta(n_snps, 0.8, 0.8)
maf <- pmax(0.05, pmin(0.5, maf))  # Constrain MAF

G <- matrix(0, n_samples, n_snps)
for (j in 1:n_snps) {
  G[, j] <- rbinom(n_samples, 2, maf[j])
}

# Select QTNs
qtn_idx <- sample(1:n_snps, n_qtn)

# =============================================================================
# Panel A: Heritability Validation
# =============================================================================

cat("Panel A: Heritability validation...\n")

simulate_phenotype_h2 <- function(G, qtn_idx, target_h2, effect_sd = 0.5) {
  n <- nrow(G)
  n_qtn <- length(qtn_idx)

  # Generate effects
  effects <- rnorm(n_qtn, 0, effect_sd)

  # Genetic values
  G_qtn <- G[, qtn_idx]
  g <- as.numeric(G_qtn %*% effects)

  # Scale to achieve target h²
  var_g <- var(g)
  if (var_g == 0) var_g <- 1

  # Residual variance to achieve h²
  var_e <- var_g * (1 - target_h2) / target_h2
  e <- rnorm(n, 0, sqrt(var_e))

  # Phenotype
  y <- g + e

  # Estimate h² (regression-based approximation)
  grm_diag <- diag(tcrossprod(scale(G_qtn))) / n_qtn
  # Simple estimate: correlation between genetic value and phenotype
  h2_est <- cor(g, y)^2

  list(y = y, g = g, h2_target = target_h2, h2_estimated = h2_est)
}

# Test multiple heritability levels
h2_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
h2_results <- data.frame()

for (h2 in h2_levels) {
  for (rep in 1:10) {
    result <- simulate_phenotype_h2(G, qtn_idx, h2)
    h2_results <- rbind(h2_results, data.frame(
      target = h2,
      estimated = result$h2_estimated,
      rep = rep
    ))
  }
}

# Summarize
h2_summary <- aggregate(estimated ~ target, h2_results,
                        function(x) c(mean = mean(x), sd = sd(x)))
h2_summary <- do.call(data.frame, h2_summary)
names(h2_summary) <- c("target", "mean", "sd")

panel_a <- ggplot(h2_summary, aes(x = target, y = mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.02, color = "#377EB8") +
  geom_point(size = 3, color = "#377EB8") +
  geom_line(color = "#377EB8") +
  labs(
    title = "A. Heritability Accuracy",
    x = expression("Target h"^2),
    y = expression("Estimated h"^2)
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1))

cat("  Mean absolute error:", round(mean(abs(h2_summary$mean - h2_summary$target)), 3), "\n")

# =============================================================================
# Panel B: Genetic Architecture Comparison
# =============================================================================

cat("\nPanel B: Genetic architecture comparison...\n")

# Additive model
simulate_additive <- function(G, qtn_idx, h2 = 0.5) {
  effects <- rnorm(length(qtn_idx), 0, 0.5)
  g <- as.numeric(G[, qtn_idx] %*% effects)
  var_e <- var(g) * (1 - h2) / h2
  y <- g + rnorm(nrow(G), 0, sqrt(var_e))
  list(y = y, var_a = var(g), var_d = 0, var_e = var_e, model = "Additive")
}

# Additive + Dominance model
simulate_add_dom <- function(G, qtn_idx, h2 = 0.5, dom_ratio = 0.3) {
  n_qtn <- length(qtn_idx)
  add_effects <- rnorm(n_qtn, 0, 0.5)
  dom_effects <- rnorm(n_qtn, 0, 0.3)

  G_qtn <- G[, qtn_idx]

  # Additive component
  g_add <- as.numeric(G_qtn %*% add_effects)

  # Dominance component (1 for heterozygotes)
  D <- (G_qtn == 1) * 1
  g_dom <- as.numeric(D %*% dom_effects)

  g <- g_add + g_dom
  var_e <- var(g) * (1 - h2) / h2
  y <- g + rnorm(nrow(G), 0, sqrt(var_e))

  list(y = y, var_a = var(g_add), var_d = var(g_dom), var_e = var_e, model = "Add + Dom")
}

# Additive + Dominance + Epistasis model
simulate_full <- function(G, qtn_idx, h2 = 0.5) {
  n_qtn <- length(qtn_idx)
  add_effects <- rnorm(n_qtn, 0, 0.4)
  dom_effects <- rnorm(n_qtn, 0, 0.2)

  G_qtn <- G[, qtn_idx]

  # Additive
  g_add <- as.numeric(G_qtn %*% add_effects)

  # Dominance
  D <- (G_qtn == 1) * 1
  g_dom <- as.numeric(D %*% dom_effects)

  # Epistasis (pairwise interactions, first 5 pairs)
  n_epi <- min(5, floor(n_qtn / 2))
  g_epi <- numeric(nrow(G))
  epi_effects <- rnorm(n_epi, 0, 0.2)
  for (i in 1:n_epi) {
    j1 <- 2 * i - 1
    j2 <- 2 * i
    g_epi <- g_epi + G_qtn[, j1] * G_qtn[, j2] * epi_effects[i]
  }

  g <- g_add + g_dom + g_epi
  var_e <- var(g) * (1 - h2) / h2
  y <- g + rnorm(nrow(G), 0, sqrt(var_e))

  list(y = y, var_a = var(g_add), var_d = var(g_dom), var_i = var(g_epi),
       var_e = var_e, model = "Add + Dom + Epi")
}

# Run simulations
arch_results <- list()
for (rep in 1:20) {
  arch_results[[length(arch_results) + 1]] <- simulate_additive(G, qtn_idx)
  arch_results[[length(arch_results) + 1]] <- simulate_add_dom(G, qtn_idx)
  arch_results[[length(arch_results) + 1]] <- simulate_full(G, qtn_idx)
}

# Extract variance components
var_data <- data.frame()
for (res in arch_results) {
  var_total <- res$var_a + ifelse(is.null(res$var_d), 0, res$var_d) +
               ifelse(is.null(res$var_i), 0, res$var_i) + res$var_e
  var_data <- rbind(var_data, data.frame(
    model = res$model,
    component = "Additive",
    variance = res$var_a / var_total
  ))
  if (!is.null(res$var_d) && res$var_d > 0) {
    var_data <- rbind(var_data, data.frame(
      model = res$model,
      component = "Dominance",
      variance = res$var_d / var_total
    ))
  }
  if (!is.null(res$var_i) && res$var_i > 0) {
    var_data <- rbind(var_data, data.frame(
      model = res$model,
      component = "Epistasis",
      variance = res$var_i / var_total
    ))
  }
  var_data <- rbind(var_data, data.frame(
    model = res$model,
    component = "Residual",
    variance = res$var_e / var_total
  ))
}

# Summarize
var_summary <- aggregate(variance ~ model + component, var_data, mean)
var_summary$component <- factor(var_summary$component,
                                 levels = c("Residual", "Epistasis", "Dominance", "Additive"))
var_summary$model <- factor(var_summary$model,
                            levels = c("Additive", "Add + Dom", "Add + Dom + Epi"))

panel_b <- ggplot(var_summary, aes(x = model, y = variance, fill = component)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Additive" = "#E41A1C", "Dominance" = "#377EB8",
                                "Epistasis" = "#4DAF4A", "Residual" = "#CCCCCC"),
                    name = "Variance\nComponent") +
  labs(
    title = "B. Genetic Architecture Models",
    x = "",
    y = "Proportion of Variance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(0, 1))

# =============================================================================
# Panel C: Multi-Trait Genetic Correlations
# =============================================================================

cat("Panel C: Multi-trait correlations...\n")

simulate_multi_trait <- function(G, qtn_idx, n_traits = 3, target_rg = 0.5, h2 = 0.5) {
  n <- nrow(G)
  n_qtn <- length(qtn_idx)
  G_qtn <- G[, qtn_idx]

  # Shared and trait-specific effects
  # rg controls the correlation of genetic effects
  shared_effects <- rnorm(n_qtn, 0, 0.5)

  Y <- matrix(0, n, n_traits)
  G_vals <- matrix(0, n, n_traits)

  for (t in 1:n_traits) {
    # Mix of shared and trait-specific
    trait_effects <- target_rg * shared_effects + sqrt(1 - target_rg^2) * rnorm(n_qtn, 0, 0.5)
    g <- as.numeric(G_qtn %*% trait_effects)
    G_vals[, t] <- g

    var_e <- var(g) * (1 - h2) / h2
    Y[, t] <- g + rnorm(n, 0, sqrt(var_e))
  }

  list(Y = Y, G = G_vals, target_rg = target_rg)
}

# Test different genetic correlations
rg_levels <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
rg_results <- data.frame()

for (rg in rg_levels) {
  for (rep in 1:10) {
    result <- simulate_multi_trait(G, qtn_idx, n_traits = 3, target_rg = rg)

    # Calculate observed genetic correlation (between genetic values)
    obs_rg_12 <- cor(result$G[, 1], result$G[, 2])
    obs_rg_13 <- cor(result$G[, 1], result$G[, 3])
    obs_rg_23 <- cor(result$G[, 2], result$G[, 3])

    rg_results <- rbind(rg_results, data.frame(
      target_rg = rg,
      observed_rg = c(obs_rg_12, obs_rg_13, obs_rg_23),
      pair = c("Trait 1-2", "Trait 1-3", "Trait 2-3"),
      rep = rep
    ))
  }
}

rg_summary <- aggregate(observed_rg ~ target_rg, rg_results,
                        function(x) c(mean = mean(x), sd = sd(x)))
rg_summary <- do.call(data.frame, rg_summary)
names(rg_summary) <- c("target", "mean", "sd")

panel_c <- ggplot(rg_summary, aes(x = target, y = mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.03, color = "#4DAF4A") +
  geom_point(size = 3, color = "#4DAF4A") +
  geom_line(color = "#4DAF4A") +
  labs(
    title = "C. Multi-Trait Genetic Correlation",
    x = expression("Target r"[G]),
    y = expression("Observed r"[G])
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  coord_fixed(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))

# =============================================================================
# Panel D: Binary Trait - Liability Threshold Model
# =============================================================================

cat("Panel D: Binary trait validation...\n")

simulate_binary_liability <- function(G, qtn_idx, prevalence = 0.1, h2_liability = 0.5) {
  n <- nrow(G)
  n_qtn <- length(qtn_idx)

  effects <- rnorm(n_qtn, 0, 0.5)
  G_qtn <- G[, qtn_idx]
  g <- as.numeric(G_qtn %*% effects)

  # Scale genetic component
  g <- scale(g)[, 1] * sqrt(h2_liability)

  # Add environmental component
  e <- rnorm(n, 0, sqrt(1 - h2_liability))

  # Liability
  liability <- g + e

  # Threshold
  threshold <- qnorm(1 - prevalence)
  case_status <- as.integer(liability > threshold)

  # Observed prevalence
  obs_prev <- mean(case_status)

  list(case = case_status, liability = liability, threshold = threshold,
       target_prev = prevalence, observed_prev = obs_prev)
}

# Test different prevalences
prev_levels <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.50)
prev_results <- data.frame()

for (prev in prev_levels) {
  for (rep in 1:20) {
    result <- simulate_binary_liability(G, qtn_idx, prevalence = prev)
    prev_results <- rbind(prev_results, data.frame(
      target = prev,
      observed = result$observed_prev,
      rep = rep
    ))
  }
}

prev_summary <- aggregate(observed ~ target, prev_results,
                          function(x) c(mean = mean(x), sd = sd(x)))
prev_summary <- do.call(data.frame, prev_summary)
names(prev_summary) <- c("target", "mean", "sd")

panel_d <- ggplot(prev_summary, aes(x = target, y = mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.02, color = "#984EA3") +
  geom_point(size = 3, color = "#984EA3") +
  geom_line(color = "#984EA3") +
  labs(
    title = "D. Binary Trait (Liability Threshold)",
    x = "Target Prevalence",
    y = "Observed Prevalence"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  coord_fixed(xlim = c(0, 0.55), ylim = c(0, 0.55))

# =============================================================================
# Panel E: Gene-Environment Interaction
# =============================================================================

cat("Panel E: Gene-environment interaction...\n")

simulate_gxe <- function(G, qtn_idx, h2 = 0.3, gxe_var = 0.1, env_type = "continuous") {
  n <- nrow(G)
  n_qtn <- length(qtn_idx)

  # Genetic effects
  effects <- rnorm(n_qtn, 0, 0.5)
  G_qtn <- G[, qtn_idx]
  g <- as.numeric(G_qtn %*% effects)
  g <- scale(g)[, 1]

  # Environment
  if (env_type == "continuous") {
    E <- rnorm(n)
  } else {
    E <- rbinom(n, 1, 0.5)
    E <- scale(E)[, 1]
  }

  # GxE interaction
  gxe <- g * E * sqrt(gxe_var)

  # Residual (what's left after G and GxE)
  total_genetic_var <- h2 + gxe_var
  var_e <- (1 - total_genetic_var) / total_genetic_var * (var(g) + var(gxe))
  e <- rnorm(n, 0, sqrt(max(0.1, var_e)))

  # Phenotype
  y <- g * sqrt(h2) + gxe + e

  # Variance decomposition
  var_total <- var(y)

  list(
    var_g = var(g * sqrt(h2)) / var_total,
    var_gxe = var(gxe) / var_total,
    var_e = var(e) / var_total,
    target_gxe = gxe_var
  )
}

# Test different GxE variances
gxe_levels <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25)
gxe_results <- data.frame()

for (gxe_var in gxe_levels) {
  for (rep in 1:20) {
    result <- simulate_gxe(G, qtn_idx, h2 = 0.3, gxe_var = gxe_var)
    gxe_results <- rbind(gxe_results, data.frame(
      target = gxe_var,
      var_g = result$var_g,
      var_gxe = result$var_gxe,
      var_e = result$var_e,
      rep = rep
    ))
  }
}

# Reshape for stacked bar
gxe_long <- data.frame()
for (gxe_var in gxe_levels) {
  subset_data <- gxe_results[gxe_results$target == gxe_var, ]
  gxe_long <- rbind(gxe_long, data.frame(
    target = gxe_var,
    component = "Genetic",
    variance = mean(subset_data$var_g)
  ))
  gxe_long <- rbind(gxe_long, data.frame(
    target = gxe_var,
    component = "GxE",
    variance = mean(subset_data$var_gxe)
  ))
  gxe_long <- rbind(gxe_long, data.frame(
    target = gxe_var,
    component = "Residual",
    variance = mean(subset_data$var_e)
  ))
}

gxe_long$component <- factor(gxe_long$component, levels = c("Residual", "GxE", "Genetic"))
gxe_long$target_label <- paste0(gxe_long$target * 100, "%")

panel_e <- ggplot(gxe_long, aes(x = factor(target), y = variance, fill = component)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Genetic" = "#E41A1C", "GxE" = "#FF7F00", "Residual" = "#CCCCCC"),
                    name = "Component") +
  scale_x_discrete(labels = paste0(gxe_levels * 100, "%")) +
  labs(
    title = "E. Gene-Environment Interaction",
    x = "Target GxE Variance",
    y = "Proportion of Variance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

# =============================================================================
# Panel F: Covariate Effects
# =============================================================================

cat("Panel F: Covariate effects...\n")

simulate_with_covariates <- function(G, qtn_idx, h2 = 0.4,
                                      age_effect = 0.1, sex_effect = 0.3,
                                      batch_var = 0.05, n_batches = 3) {
  n <- nrow(G)
  n_qtn <- length(qtn_idx)

  # Genetic component
  effects <- rnorm(n_qtn, 0, 0.5)
  G_qtn <- G[, qtn_idx]
  g <- as.numeric(G_qtn %*% effects)
  g <- scale(g)[, 1] * sqrt(h2)

  # Covariates
  age <- rnorm(n, 50, 15)
  age_scaled <- scale(age)[, 1]

  sex <- rbinom(n, 1, 0.5)

  batch <- sample(1:n_batches, n, replace = TRUE)
  batch_effects <- rnorm(n_batches, 0, sqrt(batch_var))

  # Covariate contributions
  cov_contrib <- age_scaled * age_effect + sex * sex_effect + batch_effects[batch]

  # Residual
  var_remaining <- 1 - h2 - age_effect^2 - sex_effect^2 * 0.25 - batch_var
  e <- rnorm(n, 0, sqrt(max(0.1, var_remaining)))

  # Phenotype
  y <- g + cov_contrib + e

  # Return data for plotting
  data.frame(
    y = y,
    g = g,
    age = age,
    sex = factor(sex, labels = c("Female", "Male")),
    batch = factor(batch)
  )
}

# Generate data with covariates
cov_data <- simulate_with_covariates(G, qtn_idx)

# Create subplot showing covariate effects
p_age <- ggplot(cov_data, aes(x = age, y = y)) +
  geom_point(alpha = 0.3, size = 1, color = "#377EB8") +
  geom_smooth(method = "lm", color = "#E41A1C", se = FALSE) +
  labs(x = "Age", y = "Phenotype") +
  theme_bw(base_size = 9) +
  theme(plot.margin = margin(2, 2, 2, 2))

p_sex <- ggplot(cov_data, aes(x = sex, y = y, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = c("#FF7F00", "#377EB8")) +
  labs(x = "Sex", y = "Phenotype") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))

p_batch <- ggplot(cov_data, aes(x = batch, y = y, fill = batch)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Batch", y = "Phenotype") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))

panel_f <- (p_age | p_sex | p_batch) +
  plot_annotation(title = "F. Covariate Effects",
                  theme = theme(plot.title = element_text(face = "bold", size = 11)))

# =============================================================================
# Combine Panels
# =============================================================================

cat("\nCombining panels...\n")

# Layout: 2 rows x 3 columns
fig_s2 <- (panel_a | panel_b | panel_c) / (panel_d | panel_e | panel_f) +
  plot_annotation(
    title = "Supplementary Figure S2: Advanced Phenotype Simulation Validation",
    subtitle = sprintf("Simulated population: %d samples, %d SNPs, %d QTNs", n_samples, n_snps, n_qtn),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "figS2_phenotype_simulation.pdf"), fig_s2,
       width = 14, height = 10, dpi = 300)
ggsave(file.path(out_dir, "figS2_phenotype_simulation.png"), fig_s2,
       width = 14, height = 10, dpi = 300)

cat("\nFigure S2 saved to:\n")
cat("  ", file.path(out_dir, "figS2_phenotype_simulation.pdf"), "\n")
cat("  ", file.path(out_dir, "figS2_phenotype_simulation.png"), "\n")

# =============================================================================
# Save Source Data
# =============================================================================

write.csv(h2_results, file.path(out_dir, "figS2_source_heritability.csv"), row.names = FALSE)
write.csv(var_data, file.path(out_dir, "figS2_source_architecture.csv"), row.names = FALSE)
write.csv(rg_results, file.path(out_dir, "figS2_source_correlation.csv"), row.names = FALSE)
write.csv(prev_results, file.path(out_dir, "figS2_source_binary.csv"), row.names = FALSE)
write.csv(gxe_results, file.path(out_dir, "figS2_source_gxe.csv"), row.names = FALSE)

cat("\nSource data saved.\n")
cat("\n=== Done ===\n")
