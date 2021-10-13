import time
import sys
import csv
import numpy as np
#import regex as re

#### Code ####

## Specify input file and extract sample name:
bed_file = sys.argv[1]
#sample_name = re.split("\.", bed_file) # Split on dots in name
sample_name = bed_file #sample_name[0]
out_file = sys.argv[2]

## Set total load scores and probability of homozygous positions to 0:
# Homozygous positions:
total_hom_pr = 0.0 # All homozygous positions
total_anc = 0.0 # Anc genotype
total_tv = 0.0 # All hom transversions
total_positions = 0.0 # All positions

# Load scores:
total_sift = 0.0

## Dictionary of genotypes to get anc allele probability:
genotypes = {"A": 0, "C": 4, "G": 7, "T":9}

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:

with open(bed_file, "r") as file:
        bed = csv.reader(file, delimiter = "\t")
        for line in bed:
                # Split line into scores, ancestral allele and posterior probabilities:
                post_prs = [float(i) for i in line[8:18]]
                sift = [float(i) for i in line[4:8]]
                anc = line[3]
                sifts = 0.0
                # Calculate load score for homozygous transversions:
                if anc =="A" or anc == "G":
                        sifts += (1 - sift[1]) * post_prs[4]
                        sifts += (1 - sift[3]) * post_prs[9]
                        pr_tv = post_prs[4] + post_prs[9]
                elif anc == "C" or anc == "T":
                        sifts += (1 - sift[0]) * post_prs[0]
                        sifts += (1 - sift[2]) * post_prs[7]
                        pr_tv = post_prs[0] + post_prs[7]
                total_sift += sifts
                pr_anc =  post_prs[genotypes[anc]]
                total_hom_pr += pr_anc + pr_tv
                total_anc += pr_anc
                total_tv += pr_tv
                total_positions += 1.0
                #if sifts != 0:
                        #print(line + [str(sifts)])
                #print(total_sift)
                #print(total_hom_pr)

#print(total_sift)
#print(total_hom_pr)
#print(total_sift/total_hom_pr)

with open(out_file, "w") as file:
        (file.write(sample_name + "\t" + str(total_positions) + "\t"  + str(total_hom_pr) + "\t" +
        str(total_anc) + "\t" + str(total_tv) + "\t" + str(total_sift) + "\t" +
        str(total_sift/total_hom_pr) + "\n"))
