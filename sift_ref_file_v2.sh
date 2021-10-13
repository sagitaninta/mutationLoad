#!/bin/bash

## Change folder to run script in:
CHR=chr${1:48:-8} # Extract chromosome name

mkdir $CHR
mv $1 $CHR
cd $CHR

## Filter out cds and intersect with ancestral alleles:
gunzip -c $1 | awk '$3~"CDS"' | grep -v ncRNA_gene | grep -v pseudo | \
sort -k1,1 -k4,4n | bedtools merge -i - | bedops --chop - | sed 's/^/chr/' > ${CHR}_cds_positions.bed

awk '{print $1 "\t" $2 "\t" $3 "\t"}' ${CHR}_cds_positions.bed | \
bedtools getfasta -fi ../../anc_sequences/AndeanFox_canFam3.1.fa -bed - -bedOut | \
bedtools getfasta -fi ../../anc_sequences/BlackBackJackal_canFam3.1.fa -bed - -bedOut | \
bedtools getfasta -fi ../../anc_sequences/Cat_canFam3.1.fa -bed - -bedOut | \
bedtools getfasta -fi ../../anc_sequences/Dhole_canFam3.1.fa -bed - -bedOut | \
bedtools getfasta -fi ../../anc_sequences/EthiopianWolf_canFam3.1.fa -bed - -bedOut | \
bedtools getfasta -fi ../../anc_sequences/GrayFox_canFam3.1.fa -bed - -bedOut | \
bedtools getfasta -fi ../../anc_sequences/PolarBear_canFam3.1.fa -bed - -bedOut | \
bedtools getfasta -fi ../../anc_sequences/RedFox_canFam3.1.fa -bed - -bedOut > ${CHR}_all_outgroup_cds.bed

python ../../anc_sequences/anc_seq_v1.py ${CHR}_all_outgroup_cds.bed ${CHR}_anc_cds.bed

## Get sift scores for each alternative allele:
for i in A C G T;
do \
# Convert cds positions to vep format file:
sed 's/chr//' ${CHR}_anc_cds.bed | awk '{print $1 "\t" $2+1 "\t" $3 "\t" $4"/"}' | sed "s/$/$i/"  > ${CHR}_${i}_cds.bed

# Run vep:
time vep -i ${CHR}_${i}_cds.bed --offline --cache --dir ../ \
--species "canis_lupus_familiaris" --force_overwrite --sift b --tab -o ${CHR}_${i}_vep.txt \
--fields "Location,Allele,Consequence,SIFT"

# Extract high confidence sift scores from vep output file:
grep -v "#" ${CHR}_${i}_vep.txt | grep -v "-" | sed 's/:/\t/' | sed 's/(/\t/' | sed 's/)//' | \
grep -v "low_confidence" | grep -v "tolerated" | \
awk '{print $1 "\t" $2-1 "\t" $2 "\t" $3  "\t" $6}' > ${CHR}_${i}_sift_high.bed

# Extract lowest sift score (most deleterious) for each position:
awk '{print $1 "\t" $2 "\t" $3}' ${CHR}_${i}_sift_high.bed | uniq | \
bedtools map -a - -b ${CHR}_${i}_sift_high.bed -o absmin > ${CHR}_${i}_sift_min.bed;
done

## Combine sift scores into one file:

# Create bed file with all sift positions:
bedops -u *min.bed | awk '{print $1 "\t" $2 "\t" $3}' | uniq > ${CHR}_all_sift_positions.bed

# Overlap sift scores for each allele with all positions (replace columns with no score with 1
# - will be 0 in calculation:

sed 's/^chr//' ${CHR}_anc_cds.bed | bedtools intersect -a ${CHR}_all_sift_positions.bed -b - -sorted -wb | cut -f 1-3,7 | \
bedtools intersect -a - -b ${CHR}_A_sift_min.bed -loj -sorted | cut -f 1-4,8 | \
bedtools intersect -a - -b ${CHR}_C_sift_min.bed -loj -sorted | cut -f 1-5,9 | \
bedtools intersect -a - -b ${CHR}_G_sift_min.bed -loj -sorted | cut -f 1-6,10 | \
bedtools intersect -a - -b ${CHR}_T_sift_min.bed -loj -sorted | cut -f 1-7,11 | \
sed 's/\t\./\t1/g' > ${CHR}_sift_scores.bed

# Remove tmp files:
rm cds.bed
rm A_*
rm C_*
rm G_*
rm T_*
