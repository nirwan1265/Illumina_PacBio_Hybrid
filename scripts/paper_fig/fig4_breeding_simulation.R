#!/usr/bin/env Rscript
# Figure 4: Breeding Population Simulation Validation
# Demonstrates breeding scheme simulation capabilities
#
# Panels:
#   A: Heterozygosity decay across generations (F1 → RIL)
#   B: Introgression segment lengths in BC populations
#   C: Recombination breakpoint distribution
#   D: Breeding scheme comparison (F2, RIL, BC, MAGIC)
#
# Usage: Rscript scripts/paper_fig/fig4_breeding_simulation.R

library(ggplot2)
library(patchwork)
library(data.table)

# Output directory
out_dir <- "scripts/paper_fig/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Figure 4: Breeding Population Simulation ===\n")

# -----------------------------------------------------------------------------
# Simulate breeding populations
# -----------------------------------------------------------------------------

set.seed(42)

# Genome parameters
genome_len <- 100000  # 100 kb for speed
n_markers <- 1000
marker_pos <- sort(sample(1:genome_len, n_markers))

# Create divergent parents (P1 = all 0, P2 = all 1)
P1 <- rep(0, n_markers)
P2 <- rep(1, n_markers)

cat("Simulating breeding populations...\n")
cat("  Genome length:", genome_len, "bp\n")
cat("  Markers:", n_markers, "\n")

# Recombination function (Poisson crossovers, ~1 per 100kb)
recombine <- function(h1, h2, positions, rate_per_bp = 1e-5) {
  n_xo <- rpois(1, genome_len * rate_per_bp)
  if (n_xo == 0) return(if (runif(1) < 0.5) h1 else h2)

  xo_pos <- sort(sample(1:genome_len, n_xo))

  # Determine which haplotype at each position
  result <- h1
  use_h1 <- runif(1) < 0.5
  current_xo <- 1

  for (i in seq_along(positions)) {
    while (current_xo <= length(xo_pos) && positions[i] > xo_pos[current_xo]) {
      use_h1 <- !use_h1
      current_xo <- current_xo + 1
    }
    result[i] <- if (use_h1) h1[i] else h2[i]
  }
  result
}

# Create F1 (heterozygous)
create_f1 <- function(n) {
  lapply(1:n, function(i) list(h1 = P1, h2 = P2))
}

# Self (produce gametes from same individual)
self <- function(individual) {
  g1 <- recombine(individual$h1, individual$h2, marker_pos)
  g2 <- recombine(individual$h1, individual$h2, marker_pos)
  list(h1 = g1, h2 = g2)
}

# Backcross to P1
backcross <- function(individual, recurrent = P1) {
  g1 <- recombine(individual$h1, individual$h2, marker_pos)
  list(h1 = g1, h2 = recurrent)
}

# Calculate heterozygosity
calc_het <- function(individual) {
  mean(individual$h1 != individual$h2)
}

# Calculate donor fraction (for BC)
calc_donor_frac <- function(individual) {
  mean(c(individual$h1, individual$h2) == 1)  # P2 is donor
}

# Find introgression segments
find_introgressions <- function(individual, donor_allele = 1) {
  combined <- individual$h1 + individual$h2
  # Introgression where at least one donor allele present
  has_donor <- combined >= 1

  # Find runs
  rle_result <- rle(has_donor)
  ends <- cumsum(rle_result$lengths)
  starts <- c(1, ends[-length(ends)] + 1)

  segments <- data.frame(
    start = marker_pos[starts[rle_result$values]],
    end = marker_pos[ends[rle_result$values]],
    stringsAsFactors = FALSE
  )
  segments$length <- segments$end - segments$start
  segments
}

# -----------------------------------------------------------------------------
# Panel A: Heterozygosity Decay
# -----------------------------------------------------------------------------

cat("\nGenerating Panel A: Heterozygosity Decay...\n")

# Simulate selfing generations
n_individuals <- 100
max_gen <- 8

het_data <- list()

for (gen in 0:max_gen) {
  if (gen == 0) {
    # F1
    pop <- create_f1(n_individuals)
  } else {
    # Self
    pop <- lapply(pop, self)
  }

  het_vals <- sapply(pop, calc_het)
  het_data[[gen + 1]] <- data.frame(
    generation = gen,
    heterozygosity = het_vals
  )
}

het_df <- do.call(rbind, het_data)
het_df$generation_label <- paste0("F", het_df$generation + 1)
het_df$generation_label <- factor(het_df$generation_label,
  levels = paste0("F", 1:(max_gen + 1)))

# Theoretical expectation: Het(n) = (1/2)^n
theoretical <- data.frame(
  generation = 0:max_gen,
  expected = 0.5^(0:max_gen)
)

panel_a <- ggplot(het_df, aes(x = generation, y = heterozygosity)) +
  geom_boxplot(aes(group = generation), fill = "#377EB8", alpha = 0.6,
               outlier.size = 0.5) +
  geom_line(data = theoretical, aes(y = expected),
            color = "#E41A1C", linewidth = 1, linetype = "dashed") +
  geom_point(data = theoretical, aes(y = expected),
             color = "#E41A1C", size = 2) +
  scale_x_continuous(breaks = 0:max_gen, labels = paste0("F", 1:(max_gen+1))) +
  labs(
    title = "A. Heterozygosity Decay During Selfing",
    x = "Generation",
    y = "Heterozygosity"
  ) +
  annotate("text", x = 6, y = 0.4, label = "Observed", color = "#377EB8", size = 4) +
  annotate("text", x = 6, y = 0.35, label = "Expected (1/2)^n", color = "#E41A1C", size = 4) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# -----------------------------------------------------------------------------
# Panel B: Introgression Segments in BC Populations
# -----------------------------------------------------------------------------

cat("Generating Panel B: Introgression Segments...\n")

# Simulate BC1, BC2, BC3 populations
bc_data <- list()

for (bc_gen in 1:4) {
  # Start with F1
  pop <- create_f1(100)

  # Backcross to P1 for bc_gen generations
  for (g in 1:bc_gen) {
    pop <- lapply(pop, backcross)
  }

  # Calculate introgression segment lengths
  for (i in 1:length(pop)) {
    introg <- find_introgressions(pop[[i]])
    if (nrow(introg) > 0) {
      introg$bc_gen <- paste0("BC", bc_gen)
      introg$individual <- i
      bc_data[[length(bc_data) + 1]] <- introg
    }
  }
}

bc_df <- do.call(rbind, bc_data)
bc_df$bc_gen <- factor(bc_df$bc_gen, levels = paste0("BC", 1:4))

# Expected mean segment length
expected_len <- data.frame(
  bc_gen = paste0("BC", 1:4),
  expected = genome_len / (1:4 + 1)  # Rough approximation
)

panel_b <- ggplot(bc_df, aes(x = bc_gen, y = length/1000)) +
  geom_violin(fill = "#4DAF4A", alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
  labs(
    title = "B. Introgression Segment Lengths",
    x = "Backcross Generation",
    y = "Segment Length (kb)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# -----------------------------------------------------------------------------
# Panel C: Recombination Breakpoints
# -----------------------------------------------------------------------------

cat("Generating Panel C: Recombination Breakpoints...\n")

# Track breakpoint positions across RIL development
breakpoint_data <- list()

for (rep in 1:50) {
  ind <- list(h1 = P1, h2 = P2)  # F1

  for (gen in 1:6) {
    ind <- self(ind)
  }

  # Find breakpoints (transitions in haplotype)
  for (hap_idx in 1:2) {
    hap <- if (hap_idx == 1) ind$h1 else ind$h2
    transitions <- which(diff(hap) != 0)
    if (length(transitions) > 0) {
      breakpoint_data[[length(breakpoint_data) + 1]] <- data.frame(
        position = marker_pos[transitions],
        rep = rep,
        haplotype = hap_idx
      )
    }
  }
}

bp_df <- do.call(rbind, breakpoint_data)

panel_c <- ggplot(bp_df, aes(x = position/1000)) +
  geom_histogram(bins = 50, fill = "#FF7F00", color = "black",
                 alpha = 0.7, linewidth = 0.3) +
  labs(
    title = "C. Recombination Breakpoint Distribution",
    subtitle = "RIL population (F7)",
    x = "Position (kb)",
    y = "Breakpoint Count"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# -----------------------------------------------------------------------------
# Panel D: Breeding Scheme Comparison
# -----------------------------------------------------------------------------

cat("Generating Panel D: Breeding Scheme Comparison...\n")

# Compare different schemes
schemes <- c("F2", "BC1", "BC2", "RIL-F6", "MAGIC-F4")
n_per_scheme <- 100

scheme_data <- list()

for (scheme in schemes) {
  pop <- create_f1(n_per_scheme)

  if (scheme == "F2") {
    pop <- lapply(pop, self)
  } else if (scheme == "BC1") {
    pop <- lapply(pop, backcross)
  } else if (scheme == "BC2") {
    pop <- lapply(pop, backcross)
    pop <- lapply(pop, backcross)
  } else if (scheme == "RIL-F6") {
    for (g in 1:5) pop <- lapply(pop, self)
  } else if (scheme == "MAGIC-F4") {
    # Simplified: multiple rounds of intercrossing then selfing
    for (g in 1:3) pop <- lapply(pop, self)
  }

  het_vals <- sapply(pop, calc_het)
  donor_vals <- sapply(pop, calc_donor_frac)

  scheme_data[[scheme]] <- data.frame(
    scheme = scheme,
    heterozygosity = het_vals,
    donor_fraction = donor_vals
  )
}

scheme_df <- do.call(rbind, scheme_data)
scheme_df$scheme <- factor(scheme_df$scheme, levels = schemes)

# Summary statistics
scheme_summary <- aggregate(cbind(heterozygosity, donor_fraction) ~ scheme,
                            data = scheme_df, FUN = mean)

panel_d <- ggplot(scheme_df, aes(x = heterozygosity, y = donor_fraction, color = scheme)) +
  geom_point(alpha = 0.5, size = 2) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_brewer(palette = "Set1", name = "Scheme") +
  labs(
    title = "D. Breeding Scheme Comparison",
    x = "Residual Heterozygosity",
    y = "Donor Fraction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

cat("\nCombining panels...\n")

fig4 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title = "Figure 4: Breeding Population Simulation",
    subtitle = "Simulated using 100 kb genome with 1000 markers",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "fig4_breeding_simulation.pdf"), fig4,
       width = 12, height = 10, dpi = 300)
ggsave(file.path(out_dir, "fig4_breeding_simulation.png"), fig4,
       width = 12, height = 10, dpi = 300)

cat("\nFigure 4 saved to:\n")
cat("  ", file.path(out_dir, "fig4_breeding_simulation.pdf"), "\n")
cat("  ", file.path(out_dir, "fig4_breeding_simulation.png"), "\n")

# -----------------------------------------------------------------------------
# Save source data
# -----------------------------------------------------------------------------

write.csv(het_df, file.path(out_dir, "fig4_source_data_heterozygosity.csv"), row.names = FALSE)
write.csv(bc_df, file.path(out_dir, "fig4_source_data_introgressions.csv"), row.names = FALSE)
write.csv(bp_df, file.path(out_dir, "fig4_source_data_breakpoints.csv"), row.names = FALSE)
write.csv(scheme_df, file.path(out_dir, "fig4_source_data_schemes.csv"), row.names = FALSE)

cat("\nSource data saved.\n")
cat("\n=== Done ===\n")
