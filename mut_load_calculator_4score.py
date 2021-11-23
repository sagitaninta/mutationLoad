import time
import sys
import csv
import numpy as np
import regex as re

#### Code ####

## Specify input file and extract sample name:
gpf_file = sys.argv[1]
score_file = sys.argv [2]
sample_name = re.split("\.", gpf_file) # Split on dots in name
sample_name = sample_name[0] # Take name before first dot to get sample name
out_file = sys.argv[3]
h = float(sys.argv[4])

# Intersect gpf file and ancestral file
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
                    its_score = scr_col[0:8] + gpf_col[3:13]
                    all_score.append(its_score)
    return(all_score)

print("Now intersecting gpf file and score file")

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

print("Now start calculating mutation load from the SIFT scores.")

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:
start2=time.time()
for line in bed_file:
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
        sift_homLoad = ((1-sift[1]) * post_prs[4]) + ((1-sift[2]) * post_prs[7]) + ((1-sift[3]) * post_prs[9]) # sift score for CC, GG, TT
        sift_hetDerLoad = ((1-sift[1]) * (post_prs[5]+post_prs[6])) + ((1-sift[2]) * (post_prs[5]+post_prs[8])) + ((1-sift[3]) * (post_prs[6]+post_prs[8])) # sift score for CG,CT, CG,GT, CT,GT
        sift_hetAncLoad = 0.5 * (((1-sift[1]) * post_prs[1]) + ((1-sift[2]) * post_prs[2]) + ((1-sift[3]) * post_prs[3])) # sift score for AC, AG, AT
        sift_hetLoad = h * (sift_hetDerLoad + sift_hetAncLoad)
    elif anc == "C":
        sift_homLoad = ((1-sift[0]) * post_prs[0]) + ((1-sift[2]) * post_prs[7]) + ((1-sift[3]) * post_prs[9]) # sift score for AA, GG, TT
        sift_hetDerLoad = ((1-sift[0]) * (post_prs[2]+post_prs[3])) + ((1-sift[2]) * (post_prs[2]+post_prs[8])) + ((1-sift[3]) * (post_prs[3]+post_prs[8])) # sift score for AG,AT AG,GT, AT,GT
        sift_hetAncLoad = 0.5 * (((1-sift[0]) * post_prs[1]) + ((1-sift[2]) * post_prs[5]) + ((1-sift[3]) * post_prs[6])) # sift score for AC, CG, CT
        sift_hetLoad = h * (sift_hetDerLoad + sift_hetAncLoad)
    elif anc == "G":
        sift_homLoad = ((1-sift[0]) * post_prs[0]) + ((1-sift[1]) * post_prs[4]) + ((1-sift[3]) * post_prs[9]) # sift score for AA, CC, TT
        sift_hetDerLoad = ((1-sift[0]) * (post_prs[1]+post_prs[3])) + ((1-sift[1]) * (post_prs[1]+post_prs[6])) + ((1-sift[3]) * (post_prs[3]+post_prs[6])) # sift score for AC,AT AC,CT, AT,CT
        sift_hetAncLoad = 0.5 * (((1-sift[0]) * post_prs[2]) + ((1-sift[1]) * post_prs[5]) + ((1-sift[3]) * post_prs[8])) # sift score for AG, CG, GT
        sift_hetLoad = h * (sift_hetDerLoad + sift_hetAncLoad)
    elif anc == "T":
        sift_homLoad = ((1-sift[0]) * post_prs[0]) + ((1-sift[1]) * post_prs[4]) + ((1-sift[2]) * post_prs[7]) # sift score for AA, CC, GG
        sift_hetDerLoad = ((1-sift[0]) * (post_prs[1]+post_prs[2])) + ((1-sift[1]) * (post_prs[1]+post_prs[5])) + ((1-sift[2]) * (post_prs[2]+post_prs[5])) # sift score for A(C,G), C(A,G), G(A,C)
        sift_hetAncLoad = 0.5 * (((1-sift[0]) * post_prs[3]) + ((1-sift[1]) * post_prs[6]) + ((1-sift[2]) * post_prs[8])) # sift score for AT, CT, GT
        sift_hetLoad = h * (sift_hetDerLoad + sift_hetAncLoad)
        
    total_hom_pr += pr_anc + pr_tv
    total_hom_pr_2 += pr_anc + pr_hom_der
    total_het_pr += pr_het_der + pr_het_anc
    total_sift_het += sift_hetLoad
    total_sift_hom += sift_homLoad
    total_anc += pr_anc
    total_tv += pr_tv
    total_positions += 1.0
end2=time.time()

print("Done calculating, time elapsed",end2-start2,"\nNow writing output.")
print("Total prob of homozygous positions is",total_hom_pr_2,"and with transitions removed is",total_hom_pr)
print("Total prob of heterozygous positions is",total_het_pr)

sift_out_file = out_file + "_sift_scores.txt"

with open(sift_out_file, "w") as file:
    (file.write(sample_name + "\t" + str(total_positions) + "\t"  + str(total_hom_pr) + "\t" + str(total_anc) + "\t" 
        + str(total_sift) + "\t" + str(total_sift_het) + "\t" + str(total_sift_hom) + "\t" 
        + str(total_sift/total_hom_pr) + "\t" + str(total_sift_het/total_het_pr) + "\t" + str(total_sift_hom/total_hom_pr_2) + "\t" 
        + str((total_sift_het+total_sift_hom)/(total_het_pr+total_hom_pr_2)) + "\t" + str(h) + "\t" + str("SIFT") + "\n"))
