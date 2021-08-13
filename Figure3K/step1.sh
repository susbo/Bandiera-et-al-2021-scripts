#!/bin/bash

mkdir -p tmp

# Pool tpm for the three replicates
cat data/hg19.cage_peak_phase1and2combined_tpm_ann_decoded.osc.txt.gz.extract.tsv | sed 1d | cut -f1,8-10 | sed -e 's/:/\t/g' | sed -e 's/\.\./\t/g' | sed -e 's/,/\t/g' | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"1",$5+$6+$7,$4}' > tmp/hg19.fantom5.bed

# Intersect to TSSs
bedtools window -a data/hg19.tss.bed -b tmp/hg19.fantom5.bed -Sm -w 100 > tmp/hg19.tss+fantom.bed
bedtools window -a data/hg19.tss.bed -b tmp/hg19.fantom5.bed -Sm -w 100 -v > tmp/hg19.tss-fantom.bed

# Combine tpm for separate peaks near same promoter
awk '{ c[$4]+=$11; } END {for (cc in c) print cc,c[cc]}' tmp/hg19.tss+fantom.bed > tmp/hg19.tss+fantom.cnt

# Add genes with no CAGE tags
awk '{ print $4,"0" }' tmp/hg19.tss-fantom.bed >> tmp/hg19.tss+fantom.cnt

