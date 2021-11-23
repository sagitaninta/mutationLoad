import time
import sys
import csv
import numpy as np
import regex as re

#### Code ####

## Specify input file and extract sample name:
gpf_file = sys.argv[1]
score_file = sys.argv[2]
out_file = sys.argv[3]
h = float(sys.argv[4])
sample_name = re.split("\.", gpf_file)
sample_name = sample_name[0]

# Make python script work directly with gpf file:
def intersect_score(filename1, filename2):
    all_score=[]
    with open(filename1, "r") as f1, open(filename2, "r") as f2:
        gpf=f1.readlines()
        scr=f2.readlines()
        for i in gpf:
            gpf_col=i.split("\t")
            for j in scr:
                scr_col=j.split("\t")
                if gpf_col[0:3] == scr_col[0:3]:
                    its_score = scr_col[0:5] + gpf_col[3:13]
                    all_score.append(its_score)
    return(all_score)

start1=time.time()
bed_file=intersect_score(gpf_file,score_file)
end1=time.time()
print("Intersected gpf file and score file contaning",len(bed_file),"positions")
print("Time elapsed",end1-start1)

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
phylop_anc = 0.0 # Homozygous anc
phylop_tv = 0.0

# Heterozygous positions and total number of positions:
total_pos = 0.0
total_het = 0.0
total_anc = 0.0
total_tv = 0.0
total_phylop_pos = 0.0

# Load scores:
total_phylop_tv = 0.0
total_phylop_hom = 0.0
total_phylop_het = 0.0

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:

start2=time.time()

for line in bed_file:
    # Split line into scores, ancestral allele and posterior probabilities:
    post_prs = [float(i) for i in line[5:15]] # Split out posterior probabilities
    phylop = float(line[4]) # Extract phylop score
    anc = line[3] # Extract ancestral genotype
    
    # Calculate load assuming ancestral genotypes are homozygous:
    pr_transv = geno_prob(hom_genotypes, transversions, anc) # Add probability of homozygous transversions
    pr_hom_der = geno_prob(hom_genotypes, homDer, anc) # probability of homozygous derived (should be more than hom transversions
    pr_het_der = geno_prob(het_genotypes, hetDer, anc)
    pr_het_anc = geno_prob(het_genotypes, hetAnc, anc)
        
    phylop_tvLoad = phylop * pr_transv
    phylop_homLoad = phylop * pr_hom_der
    phylop_hetLoad = h * ((phylop * pr_het_der) + (phylop * pr_het_anc * 0.5))
    
    # Calculate probability of position being ancestral homozygous or tranversion homozygous:
    hom_anc = hom_genotypes[anc] # Co-ordinate of ancestral genotype in post_pr
    pr_anc = post_prs[hom_anc] # Probability of ancestral homozygous genotype
    pr_hom = pr_transv + pr_anc # Probability transv or anc homozygous
    pr_hom_2 = pr_hom_der + pr_anc 
    pr_het = pr_het_anc + pr_het_der
    
    # Calculate probability of position
    # Add load scores and probability of being homozygous to total:
    total_phylop_tv += phylop_tvLoad
    total_phylop_het += phylop_hetLoad
    total_phylop_hom += phylop_homLoad
    if phylop != 0:
        phylop_hom_pr += pr_hom  # Add pr if has a phylop score for position
        total_phylop_pos += 1
        phylop_anc += pr_anc
        phylop_tv += pr_transv
        total_hom_pr += pr_hom # Add pr hom to hom total
        total_hom_pr_2 += pr_hom_2
        total_het_pr += pr_het
        total_hom_der += pr_hom_der
        
    # Add 1 to total positions covered count and add heterozygous and transitions count:
    total_pos += 1.0
    total_het += np.sum(post_prs[1:4] + post_prs[5:7] + post_prs[8:9])
    total_anc += pr_anc
    total_tv += pr_transv
    
end2=time.time()

print("Done calculating, time elapsed",end2-start2,"\nNow writing output.")

## Write to file:

phylop_out_file = out_file + "_phylop_scores.txt"
summary_out_file = out_file + "_summary_cons_scores.txt"

with open(phylop_out_file, "w") as file:
    (file.write(sample_name + "\t" + str(total_phylop_pos) + "\t" +  str(phylop_hom_pr) + "\t" + str(phylop_anc) + "\t" 
        + str(total_phylop_tv) + "\t" + str(total_phylop_het) + "\t" + str(total_phylop_hom) + "\t" 
        + str(total_phylop_tv/phylop_hom_pr) + "\t" + str(total_phylop_het/total_het_pr) + "\t" + str(total_phylop_hom/total_hom_pr_2) + "\t" 
        + str((total_phylop_het+total_phylop_hom)/(total_het_pr+total_hom_pr_2)) + "\t" + str(h) + "\t" + str("phyloP") + "\n"))

with open(summary_out_file, "w") as file:
    (file.write(sample_name + "\t" + str(total_pos) + "\t" + str(total_hom_pr) + "\t" + str(total_hom_pr_2) + "\t" + str(total_hom_der) + "\t" +
        str(total_anc) + "\t" + str(total_tv) + "\t" + str(total_het) + "\n"))
