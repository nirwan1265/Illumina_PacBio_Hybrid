#!/usr/bin/env Rscript

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl('^--file=', cmd)]
  if (length(file_arg)) return(normalizePath(sub('^--file=', '', file_arg[1])))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile))
  stop('Unable to determine script path')
}

this_file <- get_script_path()
script_dir <- dirname(this_file)
source(normalizePath(file.path(script_dir, '..', '..', 'inst', 'cli', 'cli_10_simulate_gwas_cohort.R')))
main_10_simulate_gwas_cohort(commandArgs(trailingOnly = TRUE))
