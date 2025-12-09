#!/bin/bash -l
#SBATCH --output=/get_beds%j.out --cpus-per-task=8 --mem=64GB --time=0-4:00:00

# Get bed files for each Puffin-D predicted dataset from R chunks created with rds_chunker.rscript.
# Experimental GROcap bed file is created separately within lift_grocap_to_t2t_positions.sh.
# These will populate the datasets_positions.tsv file used as input for summarise_repeat_nonrepeat.rscript to create the repeat/non-repeat summary table.
# Requires dataset-specific coordinates.tsv files from make_output_bed.sh.

Rscript chunks_to_positions_and_clusters.R   "T2T_Forward" "/puffin_outputs/T2T Forward/coordinates.tsv"  "/T2T_Forward_chunks" 0.10 25 forward 8

Rscript chunks_to_positions_and_clusters.R   "Reversed_Repeats"   "/puffin_outputs/T2T Reversed Repeats/coordinates.tsv"   "/Reversed_Repeats_chunks"   0.10 25 rev_repeats 8
Rscript chunks_to_positions_and_clusters.R   "Reversed_Nonrepeats"   "/puffin_outputs/T2T Reversed Non-repeats/coordinates.tsv"   "/Reversed_Nonrepeats_chunks"   0.10 25 rev_nonrepeats 8

Rscript chunks_to_positions_and_clusters.R   "Mononuc_Glob_Repeats"   "/puffin_outputs/Mononucleotide Global Shuffled Repeats/coordinates.tsv"   "/Mononuc_Glob_Repeats_chunks"   0.10 25 mononuc_glob_repeats 8
Rscript chunks_to_positions_and_clusters.R   "Mononuc_Glob_Nonrepeats"   "/puffin_outputs/Mononucleotide Global Shuffled Non-repeats/coordinates.tsv"   "/Mononuc_Glob_Nonrepeats_chunks"   0.10 25 mononuc_glob_nonrepeats 8

Rscript chunks_to_positions_and_clusters.R   "Global_Mononucleotide_Shuffled" "/puffin_outputs/Global Mononucleotide Shuffled/coordinates.tsv" "/Global_Mononucleotide_Shuffled_chunks"   0.10 25 glob_mononuc_shuffled 8


