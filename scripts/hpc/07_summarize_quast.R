#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: summarize_quast.R <quast_root_dir> <out_csv>")
}

root <- args[1]
out_csv <- args[2]

reports <- list.files(root, pattern = "report\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(reports) == 0) stop("No report.tsv found under: ", root)

read_quast <- function(path) {
  df <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  vals <- setNames(df[[2]], df[[1]])
  as.list(vals)
}

rows <- lapply(reports, function(p) {
  d <- read_quast(p)
  folder <- basename(dirname(p))
  ill <- sub("^ill([0-9]+)_pb([0-9]+)$", "\\1", folder)
  pb  <- sub("^ill([0-9]+)_pb([0-9]+)$", "\\2", folder)

  data.frame(
    ill_cov = as.integer(ill),
    pb_cov  = as.integer(pb),
    total_length = as.numeric(d[["Total length"]]),
    n_contigs = as.numeric(d[["# contigs"]]),
    n50 = as.numeric(d[["N50"]]),
    largest_contig = as.numeric(d[["Largest contig"]]),
    gc = as.numeric(d[["GC (%)"]]),
    genome_fraction = as.numeric(d[["Genome fraction (%)"]]),
    mismatches_per_100kbp = as.numeric(d[["# mismatches per 100 kbp"]]),
    indels_per_100kbp = as.numeric(d[["# indels per 100 kbp"]]),
    misassemblies = as.numeric(d[["# misassemblies"]]),
    stringsAsFactors = FALSE
  )
})

res <- do.call(rbind, rows)
res <- res[order(res$ill_cov, res$pb_cov), ]

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out_csv, row.names = FALSE)
