#!/usr/bin/env Rscript

### T2T Analysis Script ###
# Rscript T2T_Analysis.R data_name bed_file [threshold=0.1] [max_gap=25]
# E.g.: Rscript T2T_Analysis.R T2T_Forward './puffin_outputs/T2T Forward/coordinates.tsv' 0.05 10

### Setup ###

library(data.table)
library(dplyr)
library(IRanges)
library(ggplot2)
library(future.apply)
library(readr)
library(viridis)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

data_name <- args[1]
bed_file <- args[2]

# Defaults
threshold <- 0.1
max_gap <- 25

# Set parallel processing workers (edit this based on CPUs available)
workerplan <- 23

if (length(args) >= 3) {
  threshold <- as.numeric(args[3])
}
if (length(args) >= 4) {
  max_gap <- as.numeric(args[4])
}

message(paste0("Data Name = ", data_name))
message(paste0("Bed File = ", bed_file))
message(paste0("Threshold = ", threshold))
message(paste0("Max Gap = ", max_gap))

# Create Output Directories
output_directory <- paste0("./TSS_analysis_outputs/",data_name, "_thresh_", threshold, "_maxgap_", max_gap, "_output/")
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

plotting_data_directory <- paste0(output_directory, "plotting_data/")
dir.create(plotting_data_directory, showWarnings = FALSE, recursive = TRUE)

plot_directory <- paste0(output_directory, "new_plots/")
dir.create(plot_directory, showWarnings = FALSE, recursive = TRUE)


# Load T2T Bed Data
message("Loading T2T_bed data...")
T2T_bed <- fread(bed_file)

# Convert columns to numeric
T2T_bed[, Start := as.numeric(Start)]
T2T_bed[, Stop  := as.numeric(Stop)]

# Sort
setorder(T2T_bed, Chromosome, Start)

T2T_bed[, Start := as.numeric(Start)]
T2T_bed[, Stop := as.numeric(Stop)]

# Sort T2T_bed by 'Chromosome' and 'Start'
setorder(T2T_bed, Chromosome, Start)
message("T2T_bed data loaded and sorted.")

head(T2T_bed)

# Chromosome List
chromosomes <- unique(T2T_bed$Chromosome)
# For testing purposes, limit to one chromosome
#chromosomes <- chromosomes[20:21]

### Function to process a single chromosome ###

process_chromosome <- function(chrom) {
  message(paste0("Processing Chromosome: ", chrom))
  
  chrom_filenames <- T2T_bed[Chromosome == chrom, Filename]
  
  # Initialise variables for this chromosome
  total_tss_significant_chrom <- 0
  total_clusters_chrom        <- 0
  total_cluster_width_chrom   <- 0
  total_signif_width_chrom    <- 0
  total_distance_between_clusters_chrom <- 0
  num_cluster_distances_chrom <- 0

  tss_distances       <- numeric(0)
  tss_magnitudes      <- numeric(0)
  all_signif_positions <- numeric(0)
  cluster_widths      <- numeric(0)
  cluster_signif_widths <- numeric(0)
  cluster_distances   <- numeric(0)
  
  # Identify chunk files
  chunk_files <- list.files(
    path = paste0(data_name, "_chunks"),
    pattern = paste0(data_name, "_chunk_\\d+\\.rds"),
    full.names = TRUE
  )
  
  # Map each chunk file to its minimal start position for 'chrom'
  chunk_positions <- rbindlist(lapply(chunk_files, function(cf) {
    chunk_data <- readRDS(cf)
    filenames_in_chunk <- unique(chunk_data$Filename)
    
    # Merge with T2T_bed to get Start positions
    merged_data <- merge(
      data.table(Filename = filenames_in_chunk),
      T2T_bed[, .(Filename, Start, Chromosome)],
      by = "Filename",
      all.x = TRUE
    )
    # Determine if the chunk contains data for the chromosome
    if (!any(merged_data$Chromosome == chrom)) {
      return(NULL) # Skip this chunk
    }
    # Get the minimum Start position for the chromosome in this chunk
    min_start <- min(merged_data[Chromosome == chrom, Start], na.rm = TRUE)
    data.table(chunk_file = cf, min_start = min_start)
  }))
  
  # Remove NULL entries
  chunk_positions <- chunk_positions[!is.na(chunk_file)]
  
  # Sort chunks by min_start
  setorder(chunk_positions, min_start)
  
  # Get sorted chunk files
  chrom_chunk_files <- chunk_positions$chunk_file
  
  message(paste0("Chromosome ", chrom, ": chunk file count = ", length(chrom_chunk_files)))

  # Process chunks in genomic order
  for (chunk_file in chrom_chunk_files) {
    message(paste0("Processing chunk file: ", chunk_file))

    chunk_data <- readRDS(chunk_file)     # Read chunk file
    chunk_data[, index := seq_len(.N), by = Filename]     # Add index to chunk_data
    
    # Merge chunk_data with T2T_bed on Filename
    merged_data <- merge(chunk_data, T2T_bed, by = "Filename", all.x = TRUE)
    merged_data[, position := Start + index - 1] # Assign positions within each chunk
    
    # Filter for the chromosome
    chrom_data <- merged_data[Chromosome == chrom]
    message(paste0("Number of rows in chrom_data for chromosome ", chrom, ": ", nrow(chrom_data)))

    if (nrow(chrom_data) == 0) next  # Skip empty chromosomes
    
    # Identify significant positions
    chrom_data[, is_significant := `GRO_CAP (forward strand)` >= threshold]
    tss_significant <- chrom_data[is_significant == TRUE]
    
    # Update total significant TSS count
    total_tss_significant_chrom <- total_tss_significant_chrom + nrow(tss_significant)
    
    # Collect TSS magnitudes
    tss_magnitudes <- c(tss_magnitudes, tss_significant$`GRO_CAP (forward strand)`)
    
    # Collect significant positions
    signif_positions <- tss_significant$position
    
    # Append to the vector storing all significant positions for this chromosome
    all_signif_positions <- c(all_signif_positions, signif_positions)
    
    # Remove variables
    rm(chunk_data, merged_data, chrom_data, tss_significant)
  }
  
  # Calculate distances
  all_signif_positions <- unique(all_signif_positions)
  all_signif_positions <- sort(all_signif_positions)
  
  # Calculate distances between all significant positions
  if (length(all_signif_positions) > 1) {
    tss_distances <- diff(all_signif_positions)
  } else {
    tss_distances <- numeric(0)
  }

  # Clustering
  if (length(all_signif_positions) > 0) {
    signif_ranges <- IRanges(start = all_signif_positions, width = 1)
    clusters <- reduce(signif_ranges, min.gapwidth = max_gap + 1)
    
    # Collect cluster widths
    chrom_cluster_widths <- width(clusters)
    cluster_widths <- c(cluster_widths, chrom_cluster_widths)
    
    # Calculate distances between clusters
    if (length(clusters) > 1) {
      cluster_starts <- start(clusters)
      cluster_ends   <- end(clusters)
      chrom_cluster_distances <- cluster_starts[-1] - cluster_ends[-length(cluster_ends)]
      cluster_distances <- c(cluster_distances, chrom_cluster_distances)
    } else {
      chrom_cluster_distances <- numeric(0)
    }
    
    # Calculate significant widths within clusters
    overlaps <- findOverlaps(clusters, signif_ranges)
    cluster_indices <- as.data.table(overlaps)
    cluster_indices[, width := 1]  # Each significant position has a width of 1
    
    # Sum widths per cluster
    chrom_cluster_signif_widths <- cluster_indices[, .(signif_width = sum(width)), by = queryHits]$signif_width
    cluster_signif_widths <- c(cluster_signif_widths, chrom_cluster_signif_widths)
    
    # Update cluster metrics for genome-wide statistics
    total_clusters_chrom <- total_clusters_chrom + length(clusters)
    total_cluster_width_chrom <- total_cluster_width_chrom + sum(chrom_cluster_widths, na.rm = TRUE)
    total_signif_width_chrom <- total_signif_width_chrom + sum(chrom_cluster_signif_widths, na.rm = TRUE)
    total_distance_between_clusters_chrom <- total_distance_between_clusters_chrom + sum(chrom_cluster_distances, na.rm = TRUE)
    num_cluster_distances_chrom <- num_cluster_distances_chrom + length(chrom_cluster_distances)
  }
  
  # Save plotting data for this chromosome
  plotting_data <- list(
    tss_magnitudes      = tss_magnitudes,
    tss_distances       = tss_distances,
    cluster_widths      = cluster_widths,
    cluster_signif_widths = cluster_signif_widths,
    cluster_distances   = cluster_distances
  )
  saveRDS(plotting_data, file = paste0(plotting_data_directory, "plotting_data_chromosome_", chrom, ".rds"))
  
  message(paste0("Chromosome ", chrom, " processed and plotting data saved."))

  # Return chromosome-level statistics
  list(
    total_tss_significant = total_tss_significant_chrom,
    total_clusters        = total_clusters_chrom,
    total_cluster_width   = total_cluster_width_chrom,
    total_signif_width    = total_signif_width_chrom,
    total_distance_between_clusters = total_distance_between_clusters_chrom,
    num_cluster_distances = num_cluster_distances_chrom,
    positions             = all_signif_positions
  )
  
}

# Set up parallel processing plan
message("Setting up parallel processing plan")
plan(multisession, workers = workerplan)

# Process All Chromosomes
message("Starting chromosome processing")
chromosome_results <- future_lapply(chromosomes, process_chromosome)

### Aggregate Genome-Wide Statsistics ###

message("Aggregating genome-wide statistics")

# Aggregate position lists
all_positions <- unlist(lapply(chromosome_results, `[[`, "positions"))
write.table(
  all_positions,
  file = paste0(output_directory, data_name, "_positions.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
message("Wrote genome-wide TSS positions to:",
        paste0(output_directory, data_name, "_positions.txt"))

total_tss_significant <- 0
total_clusters        <- 0
total_cluster_width   <- 0
total_signif_width    <- 0
total_distance_between_clusters <- 0
num_cluster_distances <- 0

all_tss_magnitudes       <- numeric(0)
all_tss_distances        <- numeric(0)
all_cluster_widths       <- numeric(0)
all_cluster_signif_widths<- numeric(0)
all_cluster_distances    <- numeric(0)

# Aggregate results from all chromosomes (from 'chromosome_results')
for (result in chromosome_results) {
  total_tss_significant <- total_tss_significant + result$total_tss_significant
  total_clusters        <- total_clusters + result$total_clusters
  total_cluster_width   <- total_cluster_width + result$total_cluster_width
  total_signif_width    <- total_signif_width + result$total_signif_width
  total_distance_between_clusters <- total_distance_between_clusters + result$total_distance_between_clusters
  num_cluster_distances <- num_cluster_distances + result$num_cluster_distances
}

# Also aggregate the per-chromosome plotting data from .rds files
message("Aggregating plotting data from all chromosomes")

for (chrom in chromosomes) {
  plotting_data_file <- paste0(plotting_data_directory, "plotting_data_chromosome_", chrom, ".rds")
  if (file.exists(plotting_data_file)) {
    plotting_data <- readRDS(plotting_data_file)
    all_tss_magnitudes    <- c(all_tss_magnitudes,    plotting_data$tss_magnitudes)
    all_tss_distances     <- c(all_tss_distances,     plotting_data$tss_distances)
    all_cluster_widths    <- c(all_cluster_widths,    plotting_data$cluster_widths)
    all_cluster_signif_widths <- c(all_cluster_signif_widths, plotting_data$cluster_signif_widths)
    all_cluster_distances <- c(all_cluster_distances, plotting_data$cluster_distances)
    
    rm(plotting_data)
    gc()
    message(paste("Aggregated data from chromosome", chrom))
  } else {
    message(paste("Plotting data file not found for chromosome", chrom))
  }
}

message("Aggregation of plotting data completed.")

# Calculate genome-wide statistics
average_cluster_width              <- total_cluster_width / total_clusters
average_signif_width               <- total_signif_width / total_clusters
average_tss_distance               <- mean(all_tss_distances, na.rm = TRUE)
average_distance_between_clusters  <- total_distance_between_clusters / num_cluster_distances


# Print the genome-wide statistics
message(paste("Total significant TSSs:", total_tss_significant))
message(paste("Total clusters:", total_clusters))
message(paste("Average cluster width (span of clusters, bp):", round(average_cluster_width, 2), "bp"))
message(paste("Average significant cluster width (number of significant TSSs per cluster):", round(average_signif_width, 2)))
message(paste("Average TSS distance (bp):", round(average_tss_distance, 2), "bp"))
message(paste("Average distance between clusters (bp):", round(average_distance_between_clusters, 2), "bp"))

# Save the genome-wide statistics to a file
genome_stats <- data.frame(
  total_significant_TSSs             = total_tss_significant,
  total_clusters                     = total_clusters,
  average_cluster_width              = average_cluster_width,
  average_significant_cluster_width  = average_signif_width,
  average_tss_distance               = average_tss_distance,
  average_distance_between_clusters  = average_distance_between_clusters
)

write_tsv(genome_stats, file = paste0(output_directory, data_name, "_genome_wide_statistics.tsv"))
message("Genome-wide statistics saved")

### Save aggregated data ###
message("Writing aggregated data")

# Save all_tss_magnitudes
tss_magnitudes_df <- data.frame(magnitude = all_tss_magnitudes)
write_tsv(tss_magnitudes_df, file = paste0(output_directory, data_name, "_all_tss_magnitudes.tsv"))

# Save all_tss_distances
tss_distances_df <- data.frame(distance = all_tss_distances)
write_tsv(tss_distances_df, file = paste0(output_directory, data_name, "_all_tss_distances.tsv"))

# Save all_cluster_widths
cluster_widths_df <- data.frame(width = all_cluster_widths)
write_tsv(cluster_widths_df, file = paste0(output_directory, data_name, "_all_cluster_widths.tsv"))

# Save all_cluster_signif_widths
cluster_signif_widths_df <- data.frame(signif_width = all_cluster_signif_widths)
write_tsv(cluster_signif_widths_df, file = paste0(output_directory, data_name, "_all_cluster_signif_widths.tsv"))

# Save all_cluster_distances
cluster_distances_df <- data.frame(distance = all_cluster_distances)
write_tsv(cluster_distances_df, file = paste0(output_directory, data_name, "_all_cluster_distances.tsv"))

message("Finished writing aggregated plotting data")

# Print total values
total_tss_distances <- length(all_tss_distances)
message(paste("Total number of TSS distances:", total_tss_distances))
message(paste("Total significant TSSs:", total_tss_significant))

expected_total_tss_distances <- total_tss_significant - length(chromosomes)
message(paste("Expected total number of TSS distances:", expected_total_tss_distances))

sum_tss_distances <- sum(all_tss_distances)
message(paste("Sum of all TSS distances:", sum_tss_distances))

genome_size <- sum(T2T_bed$Stop - T2T_bed$Start)
message(paste("Total genome size:", genome_size))
message(paste("Ratio of sum of TSS distances to genome size:", sum_tss_distances / genome_size))

### Make Plots ###

message("Generating plots")

# Set theme and palette
theme_set(theme_classic(base_size = 14))
palette <- viridis_pal(option = "D")(5)

# Distances Between TSSs
ggplot(data.frame(distance = all_tss_distances), aes(x = distance)) +
  geom_histogram(binwidth = 0.05, fill = palette[3], color = "black") +
  scale_x_log10(labels = comma) +
  labs(
    x = "Distance to Next TSS (bp, log scale)",
    y = "Frequency",
    caption = "Figure 3: Distribution of Distances Between TSSs."
  ) +
  theme(
    plot.caption = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  filename = paste0("Puffin_", data_name, "_TSSdistances_logged.png"),
  path = plot_directory,
  width = 6, height = 4, units = "in", dpi = 300
)

message("TSS distances plot saved")

# 5) Distances Between Clusters
ggplot(data.frame(distance = all_cluster_distances), aes(x = distance)) +
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
  path = plot_directory,
  device = "pdf"
)

message("Cluster distances plot saved")
message("Done")

