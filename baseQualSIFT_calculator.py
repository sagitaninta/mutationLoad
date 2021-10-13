import time
import sys
import csv
import numpy as np
#import regex as re

#### Code ####

## Specify input file and extract sample name:
bed_file = sys.argv[1]
#sample_name = re.split("\.", bed_file) # Split on dots in name
#sample_name = sample_name[0]
out_file = sys.argv[2]

## Set total load scores and positions to 0:
total_pr_anc = 0 # All ancestral allele quality prob
total_pr_der = 0 # All variant allele quality prob
total_sifts_anc = 0 # total SIFT score for ancestral allele
total_sifts_der = 0 # total SIFT score for variant allele
total_positions = 0 # All positions
total_anc = 0 # ancestral position
total_der = 0 # variant position

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:

with open(bed_file, "r") as file:
        bed = csv.reader(file, delimiter = "\t")
        for line in bed:
                # Split line into scores, ancestral allele and posterior probabilities:
                sift = [float(i) for i in line[10:]]
                ref = line[3]
                cns = line[4]
                ref_read = float(line[5])
                var_read = float(line[6])
                ref_qual = float(line[7])
                var_qual = float(line[8])
                anc = line[9]
                # Calculate load score for ancestral allele:
                # Dropping a position if consensus is N
                if cns == "N":
                    pr_der = 0
                    pr_anc = 0
                # when the consensus is ancestral allele ...
                elif cns == anc:
                    # ... look first if we can use reads supporting ref to add to the pr_anc
                    if cns == ref:
                        pr_anc = ( 1 - np.power(10,-np.divide(var_qual,10)) * np.divide(var_read,(ref_read+var_read)) ) + (1 - np.power(10,-np.divide(ref_qual,10)) * np.divide(ref_read,(ref_read+var_read)) ) 
                    # if the ref cannot be considered, then just consider the reads supporting the variant
                    else:
                        pr_anc = 1 - np.power(10,-np.divide(var_qual,10))
                    pr_der = 0
                    total_anc += 1
                # when the consensus is not an ancestral allele ...
                else:
                    # ... still check if the reads from the refs can be used
                    if cns == ref:
                        pr_der = ( 1 - np.power(10,-np.divide(var_qual,10)) * np.divide(var_read,(ref_read+var_read)) ) + (1 - np.power(10,-np.divide(ref_qual,10)) * np.divide(ref_read,(ref_read+var_read)) )
                    # if the ref cannot be considered, then just consider the reads supporting the variant
                    else:
                        pr_der = 1 - np.power(10,-np.divide(var_qual,10))
                    pr_anc = 0
                    total_der += 1

                if cns == "A":
                    total_sifts_der += pr_der * (1-sift[0])
                    total_sifts_anc += pr_anc * (1-sift[0])
                if cns == "C":
                    total_sifts_der += pr_der * (1-sift[1])
                    total_sifts_anc += pr_anc * (1-sift[1])
                if cns == "G":
                    total_sifts_der += pr_der * (1-sift[2])
                    total_sifts_anc += pr_anc * (1-sift[2])
                if cns == "T":
                    total_sifts_der += pr_der * (1-sift[3])
                    total_sifts_anc += pr_anc * (1-sift[3])
                total_pr_anc += pr_anc
                total_pr_der += pr_der
                total_positions += 1

with open(out_file, "w") as file:
        (file.write(bed_file + "\t" + str(total_positions) + "\t"  + str(total_anc) + "\t" + str(total_der) + "\t" + str(total_pr_anc) + "\t" + str(total_pr_der) + "\t" +
        str(total_sifts_der) + "\t" + str(total_sifts_anc) + "\t" + str(total_sifts_der/(total_pr_anc + total_pr_der)) + "\n"))
