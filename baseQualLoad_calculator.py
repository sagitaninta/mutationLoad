import time
import sys
import csv
import numpy as np
#import regex as re

#### Code ####

## Specify input file and extract sample name:
bed_file = sys.argv[1] # bed file containing position in chromosome, base quality, consensus allele, ancestral allele, and phylop and phastcons scores
#sample_input = re.split("\.", bed_file) # Split on dots in name
#sample_name = re.split("\_", sample_input[0])
#sample_name = sample_name[0]
out_file = sys.argv[2]

## Set total positions to 0:
total_anc=0
total_der=0
total_phylop=0
total_phastcons_anc=0
total_phylop_anc=0
total_phastcons_der=0
total_phylop_der=0
total_phylop_pos_anc=0
total_phylop_pos_der=0
total_phastcons_pos_anc=0
total_phastcons_pos_der=0
total_pr_anc=0
total_pr_der=0
total_positions=0

## Calculate total positions of homozygous and heterozygous LoF mutations based on genotype probability
with open(bed_file, "r") as file:
    bed = csv.reader(file, delimiter = "\t")
    for line in bed:
        ref=line[3] # ref allele
        cns=line[4] # consensus allele
        ref_read=float(line[5]) # count of reads supporting ref allele in that position
        var_read=float(line[6]) # count of reads supporting variant allele in that position
        ref_qual=float(line[7]) # average base quality from reads supporting ref allele
        var_qual=float(line[8]) # average base quality from reads supporting variant allele
        anc=line[9] # ancestral allele as inferred from prev pipeline
        phylop=float(line[10]) 
        phastcons=float(line[11])
        # Calculating the probability that a position is ancestral or derived
        # Dropping a position if consensus is N
        if cns=="N":
            pr_anc=0
            pr_der=0
        # when the consensus is ancestral allele
        elif cns==anc:
            # look first if we can use reads supporting ref to add to the pr_anc
            if ref==cns:               
                pr_anc = ( 1 - np.power(10,-np.divide(var_qual,10)) * np.divide(var_read,(ref_read+var_read)) ) + (1 - np.power(10,-np.divide(ref_qual,10))  * np.divide(ref_read,(ref_read+var_read)) )
            # if there are no reads supporting ref then just consider the reads supporting the variant
            else: 
                pr_anc = 1 - np.power(10,-np.divide(var_qual,10)) 
            total_anc += 1
            total_pr_anc += pr_anc
            pr_der = 0
        # when the consensus is not ancestral allele
        else:
            # look first if reads from ref allele can support qual score (ref allele that is non-ancestral is unlikely, but just in case)
            if ref==cns:
                pr_der = ( 1 - np.power(10,-np.divide(var_qual,10)) * np.divide(var_read,(ref_read+var_read)) ) + ( 1 - np.power(10,-np.divide(ref_qual,10)) * np.divide(ref_read,(ref_read+var_read)) )
            # if ref allele is not the same, any reads supporting ref will just reduce the pr_anc
            else: 
                pr_der = ( 1 - np.power(10,-np.divide(var_qual,10)) ) 
            pr_anc = 0 
            total_der += 1
            total_pr_der += pr_der
        
        if phylop != 0:
            total_phylop_anc += phylop * pr_anc
            total_phylop_der += phylop * pr_der
            total_phylop_pos_anc += 1 * pr_anc
            total_phylop_pos_der += 1 * pr_der

        if phastcons != 0:
            total_phastcons_anc += phastcons * pr_anc 
            total_phastcons_der += phastcons * pr_der
            total_phastcons_pos_anc += 1 * pr_anc
            total_phastcons_pos_der += 1 * pr_der

        total_positions += 1


with open(out_file, "w") as file:
        (file.write(bed_file + "\t" + str(total_positions) + "\t"  + str(total_anc) + "\t" + str(total_der) + "\t" + str(total_pr_anc) + "\t" + str(total_pr_der) + "\t" +
            str(total_phylop_anc) + "\t" + str(total_phylop_der) + "\t" + str(total_phastcons_anc) + "\t" + str(total_phastcons_der) + "\t" +
            str(total_phylop_der/(total_phylop_pos_der+total_phylop_pos_anc)) + "\t" + str(total_phastcons_der/(total_phastcons_pos_der+total_phastcons_pos_anc)) + "\n"))
