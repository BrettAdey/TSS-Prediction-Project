#!/usr/bin/env bash
set -euo pipefail

# repeats.bed and nonrepeats.bed generated using make_repeat_masks.py

bedtools sort -i repeats.bed | bedtools merge -i - > t2t_repeats.merged.bed
bedtools sort -i nonrepeats.bed | bedtools merge -i - > t2t_nonrepeats.merged.bed

bgzip -cd GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz > GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
samtools faidx GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
cut -f1,2 GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.fai > grch38.chrom.sizes

# Change from NC to chr format for liftover
zcat GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz \
| awk '
  BEGIN{OFS="\t"}
  /^>/{
    acc=$1; sub(/^>/,"",acc)
    chr=""
    if ($0 ~ /chromosome [0-9XY]+/) {
      match($0,/chromosome ([0-9XY]+)/,m); chr="chr" m[1]
    }
    if (chr!="") print acc, chr
  }' \
  > t2t_nc_to_chr.tsv

# Check conversion TSV file
head t2t_nc_to_chr.tsv

# Create ucsc style (chrx) beds
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {map[$1]=$2; next}
     ($1 in map) {print map[$1], $2, $3}' t2t_nc_to_chr.tsv t2t_repeats.merged.bed > t2t_repeats.ucsc.bed

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {map[$1]=$2; next}
     ($1 in map) {print map[$1], $2, $3}'     t2t_nc_to_chr.tsv     t2t_nonrepeats.merged.bed > t2t_nonrepeats.ucsc.bed

# Sort
bedtools sort -i t2t_nonrepeats.ucsc.bed > t2t_nonrepeats.ucsc.sorted.bed
bedtools sort -i t2t_repeats.ucsc.bed    > t2t_repeats.ucsc.sorted.bed

# Liftover repeats & non-repeats
./liftOver ../reference/T2T/t2t_repeats.ucsc.sorted.merged.bed chm13v2-hg38.over.chain.gz repeats.hg38.raw.bed repeats.hg38.unmapped.bed
./liftOver ../reference/T2T/t2t_nonrepeats.ucsc.sorted.merged.bed chm13v2-hg38.over.chain.gz nonrepeats.hg38.raw.bed nonrepeats.hg38.unmapped.bed

# Sort and merge hg38 repeats & non-repeats
bedtools sort -i repeats.hg38.raw.bed | bedtools merge -i - > repeats.hg38.raw.merged.bed
bedtools sort -i nonrepeats.hg38.raw.bed | bedtools merge -i - > nonrepeats.hg38.raw.merged.bed

