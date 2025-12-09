#!/usr/bin/env Rscript

### Setup ###
library(data.table)
library(dplyr)
library(IRanges)
library(ggplot2)
library(readr)
library(viridis)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

data_name <- "Experimental_GROcap"

# Set threshold and max gap
threshold <- args
max_gap <- 25

# Create output directory
output_directory <- paste0("./TSS_analysis_outputs/", data_name, "/", data_name, "_thresh_", threshold, "_maxgap_", max_gap, "_output/")
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# Create plot directory
plot_directory <- paste0(output_directory, "new_plots/")
dir.create(plot_directory, showWarnings = FALSE, recursive = TRUE)

message("Loading GROcap bedgraph data")

bedgraph_df <- read_tsv("enbw_plus_perbp.bedgraph", col_names=c("Chromosome","Start","Stop","Transcription"))

message("GROcap data loaded")

### Function to process bedgraph data ###

process_bedgraph <- function(bedgraph_df) {
  # Ensure the dataframe has a Chromosome column
  if (!"Chromosome" %in% colnames(bedgraph_df)) {
    stop("The dataframe must contain a 'Chromosome' column")
  }

  # Define human chromosomes
  standard_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

  # Filter to keep only these chromosomes
  bedgraph_df <- bedgraph_df %>%
    filter(Chromosome %in% standard_chromosomes)

  # Extract numeric part of Chromosome name
  bedgraph_df$Chromosome_Num <- ifelse(
    bedgraph_df$Chromosome == "chrX", 23,
    ifelse(
      bedgraph_df$Chromosome == "chrY", 24,
      as.numeric(gsub("chr", "", bedgraph_df$Chromosome))
    )
  )

  # Change Chromosome to a factor and order by numeric column
  bedgraph_df$Chromosome <- factor(
    bedgraph_df$Chromosome,
    levels = unique(bedgraph_df$Chromosome[order(bedgraph_df$Chromosome_Num)])
  )
  
  # Order by Chromosome, then Start
  bedgraph_df <- bedgraph_df %>%
    arrange(Chromosome, Start)

  # Add log10(s+1) transformation for thresholding later
  if ("Transcription" %in% colnames(bedgraph_df)) {
    bedgraph_df <- bedgraph_df %>%
      mutate(`Log10(s+1)` = log10(Transcription + 1))
  }
  
  return(bedgraph_df)
}

bedgraph_processed <- process_bedgraph(bedgraph_df)

# Remove temp column (chromosome_num)
bedgraph_processed <- bedgraph_processed %>%
  select(Chromosome, Start, Stop, Transcription, `Log10(s+1)`)

# Initialise variables
chromosomes <- levels(bedgraph_processed$Chromosome)

all_tss_magnitudes <- numeric(0)
all_tss_distances <- numeric(0)
all_cluster_widths <- numeric(0)
all_cluster_signif_widths <- numeric(0)
all_cluster_distances <- numeric(0)

total_tss_significant <- 0
total_clusters <- 0
total_cluster_width <- 0
total_signif_width <- 0
total_distance_between_clusters <- 0
num_cluster_distances <- 0


# Main loop per chromosome
for (chrom in chromosomes) {
  message(paste("Processing Chromosome", chrom))

  chrom_data <- bedgraph_processed %>% 
    filter(Chromosome == chrom)

  # Filter positions below the signal threshold
  significant_data <- chrom_data %>% 
    filter(`Log10(s+1)` >= threshold)

  # Positions of significant TSSs
  signif_positions <- significant_data$Start

  # Collect TSS magnitudes & update total significant TSS count
  tss_magnitudes <- significant_data$`Log10(s+1)`
  total_tss_significant <- total_tss_significant + length(signif_positions)

  # Distances between significant TSSs
  if (length(signif_positions) > 1) {
    tss_positions_sorted <- sort(signif_positions)
    tss_dists <- diff(tss_positions_sorted)
    all_tss_distances <- c(all_tss_distances, tss_dists)
  }

  # Append TSS magnitudes
  all_tss_magnitudes <- c(all_tss_magnitudes, tss_magnitudes)

  # Clustering
  if (length(signif_positions) > 0) {
    signif_ranges <- IRanges(start = signif_positions, width = 1)
    clusters <- reduce(signif_ranges, min.gapwidth = max_gap + 1)

    # Cluster widths
    chrom_cluster_widths <- width(clusters)
    all_cluster_widths <- c(all_cluster_widths, chrom_cluster_widths)

    # Distances between clusters
    if (length(clusters) > 1) {
      cluster_starts <- start(clusters)
      cluster_ends <- end(clusters)
      chrom_cluster_distances <- cluster_starts[-1] - cluster_ends[-length(cluster_ends)]
      all_cluster_distances <- c(all_cluster_distances, chrom_cluster_distances)
    }

    # Calculate significant widths within clusters
    overlaps <- findOverlaps(clusters, signif_ranges)
    cluster_indices <- as.data.table(overlaps)
    cluster_indices[, width := 1]
    chrom_cluster_signif_widths <- cluster_indices[, .(signif_width = sum(width)), by = queryHits]$signif_width

    all_cluster_signif_widths <- c(all_cluster_signif_widths, chrom_cluster_signif_widths)

    # Update genome-wide stats
    total_clusters <- total_clusters + length(clusters)
    total_cluster_width <- total_cluster_width + sum(chrom_cluster_widths, na.rm = TRUE)
    total_signif_width <- total_signif_width + sum(chrom_cluster_signif_widths, na.rm = TRUE)
    total_distance_between_clusters <- total_distance_between_clusters + sum(chrom_cluster_distances, na.rm = TRUE)
    num_cluster_distances <- num_cluster_distances + length(chrom_cluster_distances)
  }
}



### Genome-wide statistics ###

average_cluster_width <- total_cluster_width / total_clusters
average_signif_width <- total_signif_width / total_clusters
average_distance_between_clusters <- total_distance_between_clusters / num_cluster_distances
average_tss_distance <- mean(all_tss_distances, na.rm = TRUE)

message(paste("Total significant TSSs:", total_tss_significant))
message(paste("Total clusters:", total_clusters))
message(paste("Average cluster width (span of clusters):", round(average_cluster_width, 2), "bp"))
message(paste("Average significant cluster width (number of significant TSSs per cluster):", round(average_signif_width, 2)))
message(paste("Average TSS distance:", round(average_tss_distance, 2), "bp"))
message(paste("Average distance between clusters:", round(average_distance_between_clusters, 2), "bp"))

genome_stats <- data.frame(
  total_significant_TSSs = total_tss_significant,
  total_clusters = total_clusters,
  average_cluster_width = average_cluster_width,
  average_significant_cluster_width = average_signif_width,
  average_tss_distance = average_tss_distance,
  average_distance_between_clusters = average_distance_between_clusters
)

write_tsv(genome_stats, file = paste0(output_directory, data_name, "_genome_wide_statistics.tsv"))

count_above_1e6 <- sum(all_cluster_distances > 1e6)
message(paste("Count of cluster distances above 1x10^6:", count_above_1e6))

sum_cluster_distances <- sum(all_cluster_distances)
message(paste("Sum of all TSS Cluster distances:", sum_cluster_distances))

min_cluster_distances <- min(all_cluster_distances)
message(paste("Min of all TSS Cluster distances:", min_cluster_distances))

max_cluster_distances <- max(all_cluster_distances)
message(paste("Max of all TSS Cluster distances:", max_cluster_distances))


### Save aggregated data ###
message("Writing aggregated data")

tss_magnitudes_df <- data.frame(magnitude = all_tss_magnitudes)
write_tsv(tss_magnitudes_df, file = paste0(output_directory, data_name, "_all_tss_magnitudes.tsv"))

tss_distances_df <- data.frame(distance = all_tss_distances)
write_tsv(tss_distances_df, file = paste0(output_directory, data_name, "_all_tss_distances.tsv"))

cluster_widths_df <- data.frame(width = all_cluster_widths)
write_tsv(cluster_widths_df, file = paste0(output_directory, data_name, "_all_cluster_widths.tsv"))

cluster_signif_widths_df <- data.frame(signif_width = all_cluster_signif_widths)
write_tsv(cluster_signif_widths_df, file = paste0(output_directory, data_name, "_all_cluster_signif_widths.tsv"))

cluster_distances_df <- data.frame(distance = all_cluster_distances)
write_tsv(cluster_distances_df, file = paste0(output_directory, data_name, "_all_cluster_distances.tsv"))

message("Finished writing aggregated data")

### Make Plots ###

# Distances Between TSSs
plot1 <- ggplot(data.frame(distance = all_tss_distances), aes(x = distance)) +
  geom_histogram(binwidth = 0.05, fill = palette[3], color = "black") +
  scale_x_log10(labels = comma) +
  labs(
    x = "Distance to Next TSS (bp, log scale)",
    y = "Frequency",
    caption = "Distribution of Distances Between TSSs."
  ) +
  theme(
    plot.caption = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  filename = paste0("Puffin_", data_name, "_TSSdistances_logged.png"),
  plot = plot1,
  path = plot_directory,
  width = 6, height = 4, units = "in", dpi = 300,
  device = "png"
)

message("TSS distances plot saved")

# Distances Between Clusters
plot2 <- ggplot(data.frame(distance = all_cluster_distances), aes(x = distance)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
  scale_x_log10(
    labels = comma,
    breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000),
    limits = c(25, 40005000)
  ) +
  labs(
    x = "Distance to Next Cluster (bp, log scale)",
    y = "Frequency"
  ) +
  theme(
    axis.title = element_text(face = "bold")
  )

 ggsave(
   filename = paste0("Puffin_", data_name, "_ClusterDistances_logged.pdf"),
   plot = plot2,
   path = plot_directory,
   device = "pdf"
 )

message("Cluster distances plot saved")
message("Done")
