#!/usr/bin/env Rscript

# Load Libraries and Datasets

library(data.table)
library(dplyr)
library(IRanges)
library(ggplot2)
library(future.apply)
library(readr)
library(tidyr)
library(forcats)
library(scattermore)
library(reticulate)
library(scales)

np <- import("numpy")
hprt1 <- np$load("HPRT1.npy", allow_pickle = TRUE)[[1]]
bedgraph_path <- "enbw_plus_HPRT1.bedgraph"

# Convert loaded Puffin-D prediction data into an R data.frame
hprt1_df <- data.frame()

for (category in names(hprt1)) {
  arrays <- hprt1[[category]]$arrays
  filenames <- hprt1[[category]]$filenames
  
  for (i in seq_along(arrays)) {
    array <- arrays[[i]]
    filename <- filenames[i]
    
    array_2d <- array[1, , ]
    array_df <- as.data.frame(t(array_2d))
    colnames(array_df) <- c(
      "FANTOM_CAGE (forward strand)",
      "ENCODE_CAGE (forward strand)",
      "ENCODE_RAMPAGE (forward strand)",
      "GRO_CAP (forward strand)",
      "PRO_CAP (forward strand)",
      "PRO_CAP (reverse strand)",
      "GRO_CAP (reverse strand)",
      "ENCODE_RAMPAGE (reverse strand)",
      "ENCODE_CAGE (reverse strand)",
      "FANTOM_CAGE (reverse strand)"
    )
    
    array_df$Category <- category
    array_df$Filename <- filename
    hprt1_df <- bind_rows(hprt1_df, array_df)
  }
}

rm(array_2d, array_df, array, arrays)

# Add easy names and stats columns
hprt1_df <- hprt1_df %>%
  mutate(
    Easy_Name = case_when(
      Filename == "HPRT1.npy"                   ~ "HPRT1",
      Filename == "HPRT1R.npy"                  ~ "HPRT1R",
      Filename == "HPRT1RnoCpG_AllRightPad.npy" ~ "HPRT1RnoCpG",
      TRUE                                      ~ NA_character_
    ),
    Mean_10              = rowMeans(select(., ends_with(")"))),
    Mean_5_forwardstrand = rowMeans(select(., ends_with("forward strand)"))),
    Mean_5_reversestrand = rowMeans(select(., ends_with("reverse strand)"))),
    Max_10               = apply(select(., ends_with(")")), 1, max),
    Max_5_forwardstrand  = apply(select(., ends_with("forward strand)")), 1, max),
    Max_5_reversestrand  = apply(select(., ends_with("reverse strand)")), 1, max)
  )

# Filter check for HPRT1 variants and convert log10(s+1) to pseudo-reads
hprt1_df <- hprt1_df %>%
  filter(Category %in% c("HPRT1", "HPRT1R", "HPRT1RnoCpG")) %>%
  mutate(GRO_CAP_ReadCounts = 10^(`GRO_CAP (forward strand)`) - 1) %>%
  drop_na()

# Round read values <1 to 0 and round read counts to the nearest integer
hprt1_df$GRO_CAP_Forward_Strand_Binned <- ifelse(
  hprt1_df$GRO_CAP_ReadCounts < 1,
  0,
  round(hprt1_df$GRO_CAP_ReadCounts)
)

# Factor order: HPRT1, HPRT1R, HPRT1RnoCpG
hprt1_df$Easy_Name <- factor(
  hprt1_df$Easy_Name,
  levels = c("HPRT1", "HPRT1R", "HPRT1RnoCpG")
)

### Experimental GRO-cap ###

# HPRT1 1-based inclusive window
S1 <- 134429875L
E1 <- 134529874L
stopifnot((E1 - S1 + 1L) == 100000L)

# Load BEDGraph
grocap_exp_raw <- fread(
  bedgraph_path,
  col.names = c("Chrom", "Start0", "End1", "ReadCount")
)

# Convert to 1-based closed, clip to S1-E1 window, and drop empty rows
grocap_exp_raw <- grocap_exp_raw %>%
  mutate(
    start1 = Start0 + 1L,
    end1   = End1,
    clip_start = pmax.int(start1, S1),
    clip_end   = pmin.int(end1,   E1),
    clip_width = pmax.int(0L, clip_end - clip_start + 1L)
  ) %>%
  filter(clip_width > 0L)

# Map clipped coordinates to the 1-100000 window
rng_win <- IRanges(
  start = grocap_exp_raw$clip_start - S1 + 1L,
  end   = grocap_exp_raw$clip_end   - S1 + 1L
)

# Force length to 100000 so uncovered positions are zeros
cov_rle <- coverage(rng_win, weight = grocap_exp_raw$ReadCount, width = 100000L)
grocap_vec <- as.integer(cov_rle)  # length == 100000

stopifnot(length(grocap_vec) == 100000L)

grocap_exp <- data.frame(
  Easy_Name = "HPRT1 Experimental",
  Category  = "HPRT1 Experimental",
  GRO_CAP_Forward_Strand_Binned = grocap_vec,
  stringsAsFactors = FALSE
)

# Combine with Puffin-D
puffin_plot <- hprt1_df %>%
  filter(Easy_Name %in% c("HPRT1", "HPRT1R", "HPRT1RnoCpG")) %>%
  select(Easy_Name, Category, GRO_CAP_Forward_Strand_Binned)

combined_df <- bind_rows(
  grocap_exp %>% select(Easy_Name, Category, GRO_CAP_Forward_Strand_Binned),
  puffin_plot
)

combined_df$Easy_Name <- factor(
  combined_df$Easy_Name,
  levels = c("HPRT1 Experimental", "HPRT1", "HPRT1R", "HPRT1RnoCpG")
)

# Plot raw read counts (forward strand)
prediction_method_plotname <- "ENCODE GRO_CAP"
prediction_method_outname  <- "ENCODE_GROcap"

ggplot(combined_df,
             aes(x = Easy_Name,
                 y = GRO_CAP_Forward_Strand_Binned,
                 colour = Category)) +
  geom_jitter(width = 0.4, height = 0, size = 0.4) +
  scale_y_continuous(breaks = pretty_breaks(n = 8)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none") +
  labs(x = NULL, y = "Read Count")

ggsave(
  filename = "Puffin_HPRT1_vs_Experimental_GROcap_Forward_BINNED.png",
  path = "./Puffin_Summary_Plots/post_review_HPRT1only",
  device = "png",
  dpi = 1200
)

# Get counts of no. of reads >0
combined_df %>%
  filter(Easy_Name %in% c("HPRT1 Experimental", "HPRT1")) %>%
  group_by(Easy_Name) %>%
  summarise(Count_gt1 = sum(GRO_CAP_Forward_Strand_Binned > 0),
            Total_Positions = n(),
            Fraction_gt1 = Count_gt1 / Total_Positions)

# View frequency of counts
table(combined_df$GRO_CAP_Forward_Strand_Binned[combined_df$Easy_Name == "HPRT1 Experimental"])
table(combined_df$GRO_CAP_Forward_Strand_Binned[combined_df$Easy_Name == "HPRT1"])


