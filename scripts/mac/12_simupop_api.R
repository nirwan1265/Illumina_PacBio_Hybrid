#!/usr/bin/env Rscript

# R API wrapper around SimuPOP via reticulate

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("Usage: 12_simupop_api.R --config <path> --out_prefix <path>\n")
  cat("\nConfig JSON example:\n")
  cat("{\n")
  cat("  \"population\": {\n")
  cat("    \"size\": 100,\n")
  cat("    \"ploidy\": 2,\n")
  cat("    \"loci\": [100],\n")
  cat("    \"infoFields\": [\"ind_id\"]\n")
  cat("  },\n")
  cat("  \"mating\": {\n")
  cat("    \"scheme\": \"RandomMating\",\n")
  cat("    \"offspring\": 100\n")
  cat("  },\n")
  cat("  \"generations\": 5\n")
  cat("}\n")
  quit(status = 1)
}

get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- match(flag, args)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

cfg_path <- get_arg("--config")
out_prefix <- get_arg("--out_prefix")
if (is.null(cfg_path) || is.null(out_prefix)) usage()

if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)

suppressPackageStartupMessages({
  library(jsonlite)
  library(reticulate)
})

cfg <- fromJSON(cfg_path)

# Use simitall env if available
if ("simitall" %in% conda_list()$name) {
  use_condaenv("simitall", required = TRUE)
}

source_python("scripts/py/simupop_driver.py")

res <- run_simupop(cfg, out_prefix)
cat("Wrote:\n")
cat("  ", res$genotypes, "\n", sep = "")
cat("  ", res$meta, "\n", sep = "")
