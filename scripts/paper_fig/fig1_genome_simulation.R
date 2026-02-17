#!/usr/bin/env Rscript
# Figure 1: Genome Simulation Validation
# Demonstrates repeat-stress genome generation capabilities
#
# Panels:
#   A: Genome size comparison (original vs repeat-stress levels)
#   B: GC content distribution (sliding window)
#   C: Tandem repeat locations along chromosome
#   D: Repeat density heatmap
#
# Usage: Rscript scripts/paper_fig/fig1_genome_simulation.R

library(ggplot2)
library(patchwork)
library(data.table)

# Output directory
out_dir <- "scripts/paper_fig/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Figure 1: Genome Simulation ===\n")

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

gc_content <- function(seq) {
  seq <- toupper(seq)
  gc <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
  gc / nchar(seq)
}

gc_sliding_window <- function(seq, window = 1000, step = 500) {
  seq <- toupper(seq)
  len <- nchar(seq)
  positions <- seq(1, len - window + 1, by = step)
  gc_vals <- sapply(positions, function(i) {
    subseq <- substr(seq, i, i + window - 1)
    gc_content(subseq)
  })
  data.frame(position = positions + window/2, gc = gc_vals)
}

# Simple tandem repeat finder (look for exact repeats)
find_tandem_repeats <- function(seq, min_unit = 2, max_unit = 50, min_copies = 3) {
  seq <- toupper(seq)
  len <- nchar(seq)
  repeats <- list()

  # Sample positions to check (don't check every position for speed)
  check_positions <- seq(1, len - min_unit * min_copies, by = 100)

  for (pos in check_positions) {
    for (unit_len in min_unit:min(max_unit, (len - pos) / min_copies)) {
      unit <- substr(seq, pos, pos + unit_len - 1)
      # Count consecutive copies
      copies <- 1
      check_pos <- pos + unit_len
      while (check_pos + unit_len - 1 <= len) {
        if (substr(seq, check_pos, check_pos + unit_len - 1) == unit) {
          copies <- copies + 1
          check_pos <- check_pos + unit_len
        } else {
          break
        }
      }
      if (copies >= min_copies) {
        repeats[[length(repeats) + 1]] <- data.frame(
          start = pos,
          end = pos + unit_len * copies - 1,
          unit_len = unit_len,
          copies = copies,
          total_len = unit_len * copies
        )
      }
    }
  }

  if (length(repeats) == 0) {
    return(data.frame(start = integer(), end = integer(),
                      unit_len = integer(), copies = integer(),
                      total_len = integer()))
  }
  do.call(rbind, repeats)
}

# -----------------------------------------------------------------------------
# Generate test genomes
# -----------------------------------------------------------------------------

cat("Generating test genomes...\n")

# Use E. coli as base (or create random if not available)
ecoli_path <- "inst/extdata/ref_files/Ecoli_K12_MG1655.fa"

if (file.exists(ecoli_path)) {
  ref <- read_fasta(ecoli_path)
  # Use first 500kb for speed
  ref$seq <- substr(ref$seq, 1, min(500000, nchar(ref$seq)))
  ref$len <- nchar(ref$seq)
  cat("Using E. coli K12 (first 500kb)\n")
} else {
  # Generate random genome
  set.seed(42)
  ref <- list(
    name = "random_genome",
    seq = paste(sample(c("A", "C", "G", "T"), 500000, replace = TRUE,
                       prob = c(0.25, 0.25, 0.25, 0.25)), collapse = ""),
    len = 500000
  )
  cat("Generated random 500kb genome\n")
}

# Create repeat-stress genomes at different levels
create_repeat_stress <- function(seq, n_events, seg_len, copies) {
  result <- seq
  for (i in 1:n_events) {
    len <- nchar(result)
    if (len < seg_len * 2) next

    # Pick random position
    pos <- sample(1:(len - seg_len), 1)
    segment <- substr(result, pos, pos + seg_len - 1)

    # Create tandem array
    tandem <- paste(rep(segment, copies), collapse = "")

    # Insert
    result <- paste0(
      substr(result, 1, pos - 1),
      tandem,
      substr(result, pos + seg_len, nchar(result))
    )
  }
  result
}

set.seed(123)

# Low repeat stress
cat("Creating low repeat-stress genome...\n")
low_stress <- create_repeat_stress(ref$seq, n_events = 5, seg_len = 500, copies = 3)

# Medium repeat stress
cat("Creating medium repeat-stress genome...\n")
med_stress <- create_repeat_stress(ref$seq, n_events = 15, seg_len = 1000, copies = 5)

# High repeat stress
cat("Creating high repeat-stress genome...\n")
high_stress <- create_repeat_stress(ref$seq, n_events = 30, seg_len = 2000, copies = 8)

# Collect stats
genomes <- data.frame(
  type = c("Original", "Low Stress", "Medium Stress", "High Stress"),
  length = c(nchar(ref$seq), nchar(low_stress), nchar(med_stress), nchar(high_stress)),
  gc = c(gc_content(ref$seq), gc_content(low_stress),
         gc_content(med_stress), gc_content(high_stress))
)
genomes$type <- factor(genomes$type, levels = genomes$type)
genomes$increase_pct <- (genomes$length / genomes$length[1] - 1) * 100

cat("\nGenome statistics:\n")
print(genomes)

# -----------------------------------------------------------------------------
# Panel A: Genome Size Comparison
# -----------------------------------------------------------------------------

cat("\nGenerating Panel A: Genome Size Comparison...\n")

panel_a <- ggplot(genomes, aes(x = type, y = length / 1e6, fill = type)) +
geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("+%.1f%%", increase_pct)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C")) +
  labs(
    title = "A. Genome Size After Repeat Stress",
    x = "Repeat Stress Level",
    y = "Genome Size (Mb)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, max(genomes$length) / 1e6 * 1.15))

# -----------------------------------------------------------------------------
# Panel B: GC Content Distribution
# -----------------------------------------------------------------------------

cat("Generating Panel B: GC Content Distribution...\n")

# Calculate GC in sliding windows
gc_original <- gc_sliding_window(ref$seq, window = 5000, step = 2500)
gc_original$type <- "Original"

gc_med <- gc_sliding_window(med_stress, window = 5000, step = 2500)
gc_med$type <- "Medium Stress"

gc_data <- rbind(gc_original, gc_med)
gc_data$type <- factor(gc_data$type, levels = c("Original", "Medium Stress"))

panel_b <- ggplot(gc_data, aes(x = gc, fill = type)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("#4DAF4A", "#FF7F00")) +
  labs(
    title = "B. GC Content Distribution",
    x = "GC Content (5kb windows)",
    y = "Density",
    fill = "Genome"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.4, "cm")
  )

# -----------------------------------------------------------------------------
# Panel C: Repeat Locations Along Chromosome
# -----------------------------------------------------------------------------

cat("Generating Panel C: Repeat Locations...\n")

# Find repeats in medium stress genome
repeats_med <- find_tandem_repeats(med_stress, min_unit = 100, max_unit = 2000, min_copies = 3)

if (nrow(repeats_med) > 0) {
  repeats_med$mid <- (repeats_med$start + repeats_med$end) / 2

  panel_c <- ggplot() +
    # Chromosome backbone
    geom_rect(aes(xmin = 0, xmax = nchar(med_stress), ymin = -0.3, ymax = 0.3),
              fill = "gray80", color = "black", linewidth = 0.3) +
    # Repeat locations
    geom_segment(data = repeats_med,
                 aes(x = start, xend = end, y = 0, yend = 0,
                     color = log10(total_len), linewidth = copies),
                 lineend = "round") +
    scale_color_viridis_c(option = "plasma", name = "Log10(Length)") +
    scale_linewidth_continuous(range = c(1, 4), name = "Copies") +
    scale_x_continuous(labels = function(x) paste0(x/1e6, " Mb")) +
    labs(
      title = "C. Tandem Repeat Locations (Medium Stress)",
      x = "Genomic Position",
      y = ""
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    coord_cartesian(ylim = c(-1, 1))
} else {
  # Fallback if no repeats detected
  panel_c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Tandem repeats inserted\n(visualization pending)", size = 4) +
    labs(title = "C. Tandem Repeat Locations") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 12))
}

# -----------------------------------------------------------------------------
# Panel D: Repeat Parameters Effect
# -----------------------------------------------------------------------------

cat("Generating Panel D: Repeat Parameters Effect...\n")

# Simulate different parameter combinations
param_grid <- expand.grid(
  n_events = c(5, 10, 20),
  seg_len = c(500, 1000, 2000),
  copies = 4
)

set.seed(456)
param_results <- lapply(1:nrow(param_grid), function(i) {
  params <- param_grid[i, ]
  sim_seq <- create_repeat_stress(ref$seq, params$n_events, params$seg_len, params$copies)
  data.frame(
    n_events = params$n_events,
    seg_len = params$seg_len,
    copies = params$copies,
    final_length = nchar(sim_seq),
    increase_pct = (nchar(sim_seq) / nchar(ref$seq) - 1) * 100
  )
})
param_df <- do.call(rbind, param_results)
param_df$seg_len_label <- paste0(param_df$seg_len, " bp")

panel_d <- ggplot(param_df, aes(x = factor(n_events), y = increase_pct, fill = seg_len_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Blues", name = "Segment Length") +
  labs(
    title = "D. Effect of Repeat Parameters",
    x = "Number of Duplication Events",
    y = "Genome Size Increase (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = c(0.25, 0.80),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.4, "cm")
  )

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

cat("\nCombining panels...\n")

fig1 <- (panel_a | panel_b) / (panel_c) / (panel_d) +
  plot_annotation(
    title = "Figure 1: Genome Simulation and Repeat Stress Generation",
    subtitle = "Validation using E. coli K12 MG1655 (500 kb subset)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "fig1_genome_simulation.pdf"), fig1,
       width = 10, height = 12, dpi = 300)
ggsave(file.path(out_dir, "fig1_genome_simulation.png"), fig1,
       width = 10, height = 12, dpi = 300)

cat("\nFigure 1 saved to:\n")
cat("  ", file.path(out_dir, "fig1_genome_simulation.pdf"), "\n")
cat("  ", file.path(out_dir, "fig1_genome_simulation.png"), "\n")

# -----------------------------------------------------------------------------
# Save source data
# -----------------------------------------------------------------------------

write.csv(genomes, file.path(out_dir, "fig1_source_data_genomes.csv"), row.names = FALSE)
write.csv(gc_data, file.path(out_dir, "fig1_source_data_gc.csv"), row.names = FALSE)
write.csv(param_df, file.path(out_dir, "fig1_source_data_params.csv"), row.names = FALSE)

cat("\nSource data saved.\n")
cat("\n=== Done ===\n")
