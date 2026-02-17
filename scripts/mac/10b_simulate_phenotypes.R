#!/usr/bin/env Rscript
# Mac wrapper for advanced phenotype simulation
# Uses simplePHENOTYPES with additional custom features
#
# Usage: Rscript scripts/mac/10b_simulate_phenotypes.R --geno_file <vcf> --out_prefix <path> [options]

# Pass all arguments to the CLI script
cli_script <- file.path(
  dirname(dirname(dirname(sys.frame(1)$ofile))),
  "inst", "cli", "cli_10b_simulate_phenotypes.R"
)

# If running from different location, try relative path
if (!file.exists(cli_script)) {
  # Try from project root
  cli_script <- "inst/cli/cli_10b_simulate_phenotypes.R"
}

if (!file.exists(cli_script)) {
  stop("Cannot find CLI script. Run from project root or scripts/mac/ directory.")
}

source(cli_script)
