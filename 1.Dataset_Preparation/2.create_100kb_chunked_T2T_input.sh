#!/usr/bin/env bash
#SBATCH --job-name=prep_t2t_split
#SBATCH --output=logs/prep_%j.out
#SBATCH --error=logs/prep_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail
seqkit="$HOME/bin/seqkit"

# Create 100kb/single line fasta
"$seqkit" sliding -s 100000 -W 100000 GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz | "$seqkit" seq -w 0 > T2T_100kbchunks.fasta
  
# Split by chr
mkdir -p T2T_fastas_by_chr
awk '/^>/ {split($1,a,"_sliding"); OUT="T2T_fastas_by_chr/" substr(a[1],2) ".fasta"} {print >> OUT}' T2T_100kbchunks.fasta

# Manifest: absolute paths, one per line
ls -1 /NC_*.fasta | sort > T2T_fastas_by_chr_manifest.txt


