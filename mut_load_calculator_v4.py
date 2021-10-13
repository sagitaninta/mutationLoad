import time
import sys
import csv
import numpy as np
import regex as re

#### Code ####

## Specify input file and extract sample name:
bed_file = sys.argv[1]
sift_bed = sys.argv [2]
sample_name = re.split("\.", bed_file) # Split on dots in name
sample_name = sample_name[0] # Take name before first dot to get sample name
out_file = sys.argv[3]
h = float(sys.argv[4])

## Create dictionary of positions of homozygous transversions:
## these are positions in gpf file where genotype is homozygous and transversed from ancestral
ATV = (4, 9) # CC and TT
CTV = (0, 7) # AA and GG
GTV = (4, 9) # CC and TT
TTV = (0, 7) # AA and GG

transversions = {"A":ATV, "C":CTV, "G":GTV, "T":TTV}

# hom derived
AhomDer = (4, 7, 9) # CC, GG, TT
ChomDer = (0, 7, 9) # AA, GG, TT
GhomDer = (0, 4, 9) # AA, CC, TT
ThomDer = (0, 4, 7) # AA, CC, GG

homDer = {"A":AhomDer, "C":ChomDer, "G":GhomDer, "T":ThomDer}
        
# het derived
AhetDer = (5, 6, 8) # CG, CT, GT
ChetDer = (2, 3, 8) # AG, AT, GT
GhetDer = (1, 3, 6) # AC, AT, CT
ThetDer = (1, 2, 5) # AC, AG, CG

hetDer = {"A":AhetDer, "C":ChetDer, "G":GhetDer, "T":ThetDer}

# het anc
AhetAnc = (1, 2, 3) # AC, AG, AT
ChetAnc = (1, 5, 6) # AC, CG, CT
GhetAnc = (2, 5, 8) # AG, CG, GT
ThetAnc = (3, 6, 8) # AT, CT, GT

hetAnc = {"A":AhetAnc, "C":ChetAnc, "G":GhetAnc, "T":ThetAnc}

## Create dictionary of position of each homozygous genotype:
hom_genotypes = {"A": 0, "C": 4, "G": 7, "T":9}
het_genotypes = {"A": [1, 2, 3], "C": [5, 6, 7], "G": [2, 5, 8], "T": [3, 6, 8]}

def geno_prob(geno,state,anc):
    pr=0
    for geno in state[anc]:
        pr += post_prs[geno]
    return(pr)

## Set total load scores and probability of homozygous positions to 0:
# Homozygous positions:
total_hom_pr = 0.0 # All homozygous positions
total_hom_pr_2 = 0
total_het_pr = 0.0
total_hom_der = 0.0
phylop_hom_pr = 0.0 # Homozygous positions with phylop score
phastcons_hom_pr = 0.0 # Homozygous positions with phastcons score

phylop_anc = 0.0 # Homozygous anc
phastcons_anc = 0.0

phylop_tv = 0.0
phastcons_tv = 0.0 # Homozygous transversion

# Heterozygous positions and total number of positions:
total_pos = 0.0
total_het = 0.0
total_anc = 0.0
total_tv = 0.0

total_phylop_pos = 0.0
total_phastcons_pos = 0.0

total_phylop_het = 0.0
total_phastcons_het = 0.0

# Load scores:
total_phylop = 0.0
total_phastcons = 0.0
total_phylop_hom = 0.0
total_phastcons_hom = 0.0
total_phylop_het = 0.0
total_phastcons_het = 0.0

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:

with open(bed_file, "r") as file:
        bed = csv.reader(file, delimiter = "\t")
        for line in bed:
                if line[3] in ["A", "C", "G", "T"]:
                # Split line into scores, ancestral allele and posterior probabilities:
                        post_prs = [float(i) for i in line[6:16]] # Split out posterior probabilities
                        phylop = float(line[4]) # Extract phylop score
                        phastcons = float(line[5]) # Extract phastcons score
                        anc = line[3] # Extract ancestral genotype

                # Calculate load assuming ancestral genotypes are homozygous:
                        pr_transv = geno_prob(hom_genotypes, transversions, anc) # Add probability of homozygous transversions
                        pr_hom_der = geno_prob(hom_genotypes, homDer, anc) # probability of homozygous derived (should be more than hom transversions)
                        pr_het_der = geno_prob(het_genotypes, hetDer, anc)
                        pr_het_anc = geno_prob(het_genotypes, hetAnc, anc)

                        phylop_load = phylop * pr_transv
                        phastcons_load = phastcons * pr_transv

                        phylop_homLoad = phylop * pr_hom_der
                        phastcons_homLoad = phastcons * pr_hom_der

                        phylop_hetLoad = h * ((phylop * pr_het_der) + (phylop * pr_het_anc * 0.5))
                        phastcons_hetLoad = h * ((phastcons * pr_het_der) + (phastcons * pr_het_anc * 0.5))

                # Calculate probability of position being ancestral homozygous or tranversion homozygous:
                        hom_anc = hom_genotypes[anc] # Co-ordinate of ancestral genotype in post_pr

                        pr_anc = post_prs[hom_anc] # Probability of ancestral homozygous genotype
                        pr_hom = pr_transv + pr_anc # Probability transv or anc homozygous
                        pr_hom_2 = pr_hom_der + pr_anc 
                        pr_het = pr_het_anc + pr_het_der

                        # Calculate probability of position

                # Add load scores and probability of being homozygous to total:
                        total_phylop += phylop_load
                        total_phastcons += phastcons_load
                        total_phylop_het += phylop_hetLoad
                        total_phastcons_het += phastcons_hetLoad
                        total_phylop_hom += phylop_homLoad
                        total_phastcons_hom += phastcons_homLoad
                        if phylop != 0:
                                phylop_hom_pr += pr_hom  # Add pr if has a phylop score for position
                                total_phylop_pos += 1
                                phylop_anc += pr_anc
                                phylop_tv += pr_transv
                        if phastcons != 0:
                                phastcons_hom_pr += pr_hom # Add pr if has a phastcons score for position
                                total_phastcons_pos += 1
                                phastcons_anc += pr_anc
                                phastcons_tv += pr_transv
                        total_hom_pr += pr_hom # Add pr hom to hom total
                        total_hom_pr_2 += pr_hom_2
                        total_het_pr += pr_het
                        total_hom_der += pr_hom_der

                # Add 1 to total positions covered count and add heterozygous and transitions count:
                        total_pos += 1.0
                        total_het += np.sum(post_prs[1:4] + post_prs[5:7] + post_prs[8:9])
                        total_anc += pr_anc
                        total_tv += pr_transv
## Write to file:

phylop_out_file = out_file + "_phylop_scores.txt"
phastcons_out_file = out_file + "_phastcons_scores.txt"
summary_out_file = out_file + "_summary_cons_scores.txt"

with open(phylop_out_file, "w") as file:
    (file.write(sample_name + "\t" + str(total_phylop_pos) + "\t" +  str(phylop_hom_pr) + "\t" + str(phylop_anc) + "\t" + str(phylop_tv) + "\t" + str(total_phylop_het) + "\t" + str(total_phylop) + "\t" + str(total_phylop/phylop_hom_pr) + 
        "\t" + str(total_phylop_het/total_het_pr) + "\t" + str(total_phylop_hom/total_hom_pr_2) + "\t" + str((total_phylop_het+total_phylop_hom)/(total_het_pr+total_hom_pr_2)) + "\t" + str(h) + "\n"))

with open(phastcons_out_file, "w") as file:
    (file.write(sample_name + "\t" + str(total_phastcons_pos)  + "\t" + str(phastcons_hom_pr) + "\t" + str(phastcons_anc) + "\t" + str(phastcons_tv) + "\t" + str(total_phastcons_het) + "\t" + str(total_phastcons) + "\t" + str(total_phastcons/phastcons_hom_pr) + 
        "\t" + str(total_phastcons_het/total_het_pr) + "\t" + str(total_phastcons_hom/total_hom_pr_2) + "\t" + str((total_phastcons_het+total_phastcons_hom)/(total_het_pr+total_hom_pr_2)) + "\t" + str(h) + "\n"))

with open(summary_out_file, "w") as file:
    (file.write(sample_name + "\t" + str(total_pos) + "\t" + str(total_hom_pr) + "\t" + str(total_hom_pr_2) + "\t" + str(total_hom_der) + "\t" +
        str(total_anc) + "\t" + str(total_tv) + "\t" + str(total_het) + "\n"))

## ------------------------------------------------------------------------------------- ##
## CALCULATING SIFT SCORES
## ------------------------------------------------------------------------------------- ##

# Set total load scores and probability of homozygous positions to 0:
# Homozygous positions:
total_hom_pr = 0.0 # All homozygous positions
total_anc = 0.0 # Anc genotype
total_tv = 0.0 # All hom transversions
total_positions = 0.0 # All positions
total_hom_pr_2 = 0.0
total_het_pr = 0.0
total_sift_het = 0.0
total_sift_hom = 0.0

# Load scores:
total_sift = 0.0

sift_der = {"A":(1,2,3), "C":(0,2,3), "G":(0,1,3), "T":(0,1,2)}

def sift_load(allele,state,anc):
    s=0
    for allele in state[anc]:
        s += sift[allele]
    return(s)

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:

with open(sift_bed, "r") as file:
        bed = csv.reader(file, delimiter = "\t")
        for line in bed:
                # Split line into scores, ancestral allele and posterior probabilities:
                post_prs = [float(i) for i in line[8:18]]
                sift = [float(i) for i in line[4:8]]
                anc = line[3]
                sifts = 0.0
                
                # Calculate load score for homozygous transversions:
                if anc =="A" or anc == "G":
                        sifts += (1 - sift[1]) * post_prs[4] # C's sift scores times CC geno prob
                        sifts += (1 - sift[3]) * post_prs[9] # T's sift scores times TT geno prob
                        pr_tv = post_prs[4] + post_prs[9]
                elif anc == "C" or anc == "T":
                        sifts += (1 - sift[0]) * post_prs[0] # A's sift scores times AA geno prob
                        sifts += (1 - sift[2]) * post_prs[7] # G's sift scores times GG geno prob
                        pr_tv = post_prs[0] + post_prs[7]
                total_sift += sifts
                pr_anc =  post_prs[hom_genotypes[anc]] # probability of homozygous ancestral
                pr_transv = geno_prob(hom_genotypes, transversions, anc) # Add probability of homozygous transversions
                pr_hom_der = geno_prob(hom_genotypes, homDer, anc) # probability of homozygous derived (should be more than hom transversions)
                pr_het_der = geno_prob(het_genotypes, hetDer, anc)
                pr_het_anc = geno_prob(het_genotypes, hetAnc, anc)
                
                if anc == "A":
                    sift_score = sift_load(sift, sift_der, "A")
                    sift_homLoad = pr_hom_der * (1-sift_score)
                    sift_hetLoad = h * (((1-sift_score) * pr_het_der) + ((1-sift_score) * pr_het_anc * 0.5))
                elif anc == "C":
                    sift_score = sift_load(sift, sift_der, "C")
                    sift_homLoad = pr_hom_der * (1-sift_score)
                    sift_hetLoad = h * (((1-sift_score) * pr_het_der) + ((1-sift_score) * pr_het_anc * 0.5))
                elif anc == "G":
                    sift_score = sift_load(sift, sift_der, "G")
                    sift_homLoad = pr_hom_der * (1-sift_score)
                    sift_hetLoad = h * (((1-sift_score) * pr_het_der) + ((1-sift_score) * pr_het_anc * 0.5))
                elif anc == "T":
                    sift_score = sift_load(sift, sift_der, "T")
                    sift_homLoad = pr_hom_der * (1-sift_score)
                    sift_hetLoad = h * (((1-sift_score) * pr_het_der) + ((1-sift_score) * pr_het_anc * 0.5))

                total_hom_pr += pr_anc + pr_tv
                total_hom_pr_2 += pr_anc + pr_hom_der
                total_het_pr += pr_het_der + pr_het_anc
                total_sift_het += sift_hetLoad
                total_sift_hom += sift_homLoad
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

sift_out_file = out_file + "_sift_scores.txt"

with open(sift_out_file, "w") as file:
    (file.write(sample_name + "\t" + str(total_positions) + "\t"  + str(total_hom_pr) + "\t" + str(total_hom_pr_2) + "\t" +
        str(total_anc) + "\t" + str(total_tv) + "\t" + str(total_sift) + "\t" + str(total_sift/total_hom_pr) + "\t" + 
        str(total_sift_het/total_het) + "\t" + str(total_sift_hom/total_hom_pr_2) + "\t" + str((total_sift_het+total_sift_hom)/(total_het_pr+total_hom_pr_2)) + "\t" + str(h) + "\n"))
