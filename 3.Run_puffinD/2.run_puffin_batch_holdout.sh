#!/bin/bash -l
#SBATCH --output=/puffin/logs/run_puffin%j.out --ntasks=1 --cpus-per-task=8 --mem=64GB --time=0-5:00:00 --mail-type=END,FAIL

# Run Puffin on a single 100kb region in a Puffin-D holdout chromosome
python puffin_D.py sequence chr8:89,892,879-89,992,878.fasta
mv chr8_89,892,879-89,992,878.npy "./puffin_outputs/Chr8 Holdout"
