import time
import sys
import csv
import numpy as np
#import regex as re

#### Code ####

## Specify input file and extract sample name:
load_file = sys.argv[1] # bed file containing position in chromosome, base quality, consensus allele, ancestral allele, and phylop and phastcons scores
sift_file = sys.argv[2]
#sample_input = re.split("\.", bed_file) # Split on dots in name
#sample_name = re.split("\_", sample_input[0])
#sample_name = sample_name[0]
out_file = sys.argv[3] # output file prefix

## -----------------------------------------------------------
## Calculation with the file with phylop and phastcons scores
## -----------------------------------------------------------

## Set total positions to 0:
pr_anc=0
pr_der=0
total_anc=0
total_der=0
total_phylop_pos_anc=0
total_phylop_pos_der=0
total_phastcons_pos_anc=0
total_phastcons_pos_der=0
total_positions=0
## Set load score to 0:
total_phastcons_anc=0
total_phylop_anc=0
total_phastcons_der=0
total_phylop_der=0

## Calculate total positions of homozygous and heterozygous LoF mutations based on genotype probability
with open(load_file, "r") as file:
    bed = csv.reader(file, delimiter = "\t")
    for line in bed:
        hap=line[3] # haploid allele
        anc=line[4] # ancestral allele as inferred from prev pipeline
        phylop=float(line[5]) 
        phastcons=float(line[6])
        # Calculating the probability that a position is ancestral or derived
        # When the random haploid allele is ancestral
        if hap==anc:
            pr_anc = 1
            pr_der = 0
            total_anc += 1
        # when the random haploid allele is not ancestral
        else:
            pr_anc = 0
            pr_der = 1
            total_der += 1
        
        if phylop != 0:
            total_phylop_anc += phylop * pr_anc
            total_phylop_der += phylop * pr_der
            total_phylop_pos_der += pr_der
            total_phylop_pos_anc += pr_anc

        if phastcons != 0:
            total_phastcons_anc += phastcons * pr_anc 
            total_phastcons_der += phastcons * pr_der
            total_phastcons_pos_der += pr_der
            total_phastcons_pos_anc += pr_anc

        total_positions += 1

# Write output file names
load_out_file = out_file + "_phy_pha_angsdHapload.txt"

## Writing output for phastcons, phylop, and other single-score based load
with open(load_out_file, "w") as file:
        (file.write(load_file + "\t" + str(total_positions) + "\t"  + str(total_anc) + "\t" + str(total_der) + "\t" +
            str(total_phylop_anc) + "\t" + str(total_phylop_der) + "\t" + str(total_phastcons_anc) + "\t" + str(total_phastcons_der) + "\t" +
            str(total_phylop_der/(total_phylop_pos_der+total_phylop_pos_anc)) + "\t" + str(total_phastcons_der/(total_phastcons_pos_der+total_phastcons_pos_anc)) + "\n"))

        
## -----------------------------------------------------
## Calculation with the file with SIFT scores
## -----------------------------------------------------

## Set total load scores and positions to 0:
total_pr_anc = 0 # All ancestral allele
total_pr_der = 0 # All variant allele 
total_sifts_anc = 0 # total SIFT score for ancestral allele
total_sifts_der = 0 # total SIFT score for variant allele
total_positions = 0 # All positions

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:

with open(sift_file, "r") as file:
        bed = csv.reader(file, delimiter = "\t")
        for line in bed:
                # Split line into scores, ancestral allele, and haploid allele
                sift = [float(i) for i in line[5:]]
                hap = line[3]
                anc = line[4]
                # Calculate load score for ancestral allele:
                # when the consensus is ancestral allele ...
                if hap == anc:
                    pr_der = 0
                    pr_anc = 1
                    total_pr_anc += 1
                    total_pr_der += 0
                # when the consensus is not an ancestral allele ...
                else:
                    pr_der = 1
                    pr_anc = 0
                    total_pr_anc += 0
                    total_pr_der += 1

                if hap == "A":
                    total_sifts_der += pr_der * (1-sift[0])
                    total_sifts_anc += pr_anc * (1-sift[0])
                if hap == "C":
                    total_sifts_der += pr_der * (1-sift[1])
                    total_sifts_anc += pr_anc * (1-sift[1])
                if hap == "G":
                    total_sifts_der += pr_der * (1-sift[2])
                    total_sifts_anc += pr_anc * (1-sift[2])
                if hap == "T":
                    total_sifts_der += pr_der * (1-sift[3])
                    total_sifts_anc += pr_anc * (1-sift[3])
                
                total_positions += 1

# Write output filename
sift_out_file = out_file + "_sift_angsdHapload.txt"

## Writing output for SIFT
with open(sift_out_file, "w") as file:
        (file.write(sift_file + "\t" + str(total_positions) + "\t" + str(total_pr_anc) + "\t" + str(total_pr_der) + "\t" +
        str(total_sifts_der) + "\t" + str(total_sifts_anc) + "\t" + str(total_sifts_der/(total_pr_anc + total_pr_der)) + "\n"))
