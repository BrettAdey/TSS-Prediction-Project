#!/usr/bin/env bash
set -euo pipefail

# inputs
BG=./enbw_plus.bedgraph
CHAIN=./hg38-chm13v2.over.chain.gz
CHR2NC=./chr_to_nc.tsv

# outputs
BG_PERBP=./enbw_plus_perbp.bedgraph
HG38_BED=GROcap_hg38_pos.bed
LIFTED_CHM13=GROcap_chm13_chr.bed
UNMAPPED=unmapped.txt
LIFTED_NC=GROcap_chm13_nc.bed
LIFTED_NC_SORT=GROcap_chm13_nc_sorted.bed
LIFTED_NC_PERBP=GROcap_chm13_nc_perbp.bed

# Expand to per-bp positions (Note: this is also used for Experimental GROcap TSS distance analysis)
awk '{for(i=$2;i<$3;i++) print $1"\t"i"\t"i+1"\t"$4}' "$BG" > "$BG_PERBP"

# Convert to BED
awk '{print $1"\t"$2"\t"$3}' "$BG_PERBP" > "$HG38_BED"

# liftOver to T2T
./liftOver "$HG38_BED" "$CHAIN" "$LIFTED_CHM13" "$UNMAPPED"

# Convert chr -> NC_ ids
awk 'NR==FNR{map[$1]=$2; next} {if($1 in map){$1=map[$1]} print}' OFS='\t' \
    "$CHR2NC" "$LIFTED_CHM13" > "$LIFTED_NC"

# Sort
sort -k1,1 -k2,2n -k3,3n "$LIFTED_NC" > "$LIFTED_NC_SORT"
