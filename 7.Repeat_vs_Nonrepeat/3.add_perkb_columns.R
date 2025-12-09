#!/usr/bin/env Rscript

# add_perkb_columns.R
# Rscript add_perkb_columns.R t2t_repeats.bed t2t_nonrepeats.bed in_summary.tsv out_tsv [grch38_repeats.{bed}] [grch38_nonrepeats.{bed}]

library(data.table)

args <- commandArgs(trailingOnly=TRUE)

t2t_rep_bed <- args[1]
t2t_non_bed <- args[2]
in_tsv      <- args[3]
out_tsv     <- args[4]

have_grch38 <- length(args) >= 6
grch38_rep_arg <- if (have_grch38) args[5] else NA_character_
grch38_non_arg <- if (have_grch38) args[6] else NA_character_

# Function to get total bp from .bed path
bp_from_arg <- function(x) {
    # read BED and sum spans
    bed <- fread(x, col.names=c("chr","start","end"))
    sum(as.numeric(bed$end) - as.numeric(bed$start))
}

# T2T denominators (default for all rows)
t2t_rep_bp <- bp_from_arg(t2t_rep_bed)
t2t_non_bp <- bp_from_arg(t2t_non_bed)
t2t_rep_kb <- t2t_rep_bp / 1000
t2t_non_kb <- t2t_non_bp / 1000

# GRCh38 denominators (used only for Experimental GRO-cap)
grch38_rep_bp <- if (have_grch38) bp_from_arg(grch38_rep_arg) else NA_real_
grch38_non_bp <- if (have_grch38) bp_from_arg(grch38_non_arg) else NA_real_
grch38_rep_kb <- if (!is.na(grch38_rep_bp)) grch38_rep_bp / 1000 else NA_real_
grch38_non_kb <- if (!is.na(grch38_non_bp)) grch38_non_bp / 1000 else NA_real_

# Load summary table
dt <- fread(in_tsv, sep="\t")
# Set dataset as character
dt[, dataset := as.character(dataset)]

# Compute per case length, switching denominators to grch38 for Experimental GRO-cap
if (!is.na(grch38_rep_kb) && !is.na(grch38_non_kb)) {
  # Default T2T per-kb
  dt[, repeat_per_kb    := n_repeat    / t2t_rep_kb]
  dt[, nonrepeat_per_kb := n_nonrepeat / t2t_non_kb]
  # Override for Experimental GRO-cap
  dt[dataset == "Experimental_GROcap", `:=`(
    repeat_per_kb    = n_repeat    / grch38_rep_kb,
    nonrepeat_per_kb = n_nonrepeat / grch38_non_kb
  )]
  msg <- sprintf(
    paste0("Denominators used:\n",
           "  T2T:     repeat_kb=%f  nonrepeat_kb=%f\n",
           "  GRCh38*: repeat_kb=%f  nonrepeat_kb=%f  (Experimental GROcap only)\n"),
    t2t_rep_kb, t2t_non_kb, grch38_rep_kb, grch38_non_kb)
} else {
  # Use T2T for all rows if no grch38 input
  dt[, repeat_per_kb    := n_repeat    / t2t_rep_kb]
  dt[, nonrepeat_per_kb := n_nonrepeat / t2t_non_kb]
  msg <- sprintf(
    "Denominators used:\n  T2T: repeat_kb=%f  nonrepeat_kb=%f\n  (No GRCh38 denominators provided; used T2T for all rows)\n",
    t2t_rep_kb, t2t_non_kb)
}

# Write
fwrite(dt, out_tsv, sep="\t")
cat("Wrote:", out_tsv, "\n")
cat(msg)
