#!/usr/bin/env bash
set -euo pipefail

# Create 'coordinates.tsv' files for each Puffin-D output

#E.g:
#bash make_output_bed.sh "Local Mononucleotide Shuffled"
#bash make_output_bed.sh "Local Dinucleotide Shuffled"
#bash make_output_bed.sh "T2T Reversed Repeats"
#bash make_output_bed.sh "Mononucleotide Global Shuffled Repeats"
#bash make_output_bed.sh "Mononucleotide Global Shuffled Non-repeats"
#bash make_output_bed.sh "T2T Holdout Forward"
#bash make_output_bed.sh "T2T Holdout Reversed"


run="$1"

cd "./puffin_outputs/${run}/"

out="coordinates.tsv"
printf 'Filename\tChromosome\tStart\tStop\n' > "$out"

for filename in *.npy; do
  base="${filename%.npy}"        # Strip extension
  chrom="${base%%_sliding_*}"    # Prefix before “_sliding_”
  coords="${base#*_sliding_}"    # Suffix after “_sliding_”
  coords="${coords%%_*}"         # Drop any “_…_shuffled” suffix
  start="${coords%-*}"           # Before the hyphen
  stop="${coords#*-}"            # After the hyphen

  printf '%s\t%s\t%s\t%s\n' \
    "$filename" "$chrom" "$start" "$stop" \
    >> "$out"
done


