import time
import sys
import csv
import numpy as np

start = time.time()

###### Functions ########

def prior_calculator(heterozygosity):
	"Calculate hom and het priors from input heterozygosity value"
	#np.set_printoptions(linewidth=np.inf) #Set numpy print options if needed
	het = np.divide(heterozygosity, 6.0)
	q = np.subtract(1.0, heterozygosity)
	hom = np.divide(q, 4.0)
	priors = [hom, het, het, het, hom, het, het, hom, het, hom]
	return priors



def float_gls(line):
	"Extract gls and convert from strings to floats"
	gls = [float (i) for i in line[2:12]]
	return gls



def extract_best_genotype(gls):
	"Return genotype with likelihood ratio = 0.0 (best genotype)"
	genotypes = ["AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"]   #List genotypes alphabetically
	for value in range(0,10):
		if gls[value] == 0.0:
			best_geno = genotypes[value]     # Store genotype with highest gl
	return best_geno



def gl_calculator(gls, priors = [0.1]*10):               # Default = uniform prior
	"Calculate genotype likelihoods for all genotypes"
	#np.set_printoptions(linewidth=np.inf)
	values = np.power(10, gls)            # Raise 10 to power of gls
	values = np.multiply(values, priors)   # Multiply by prior
	sum_prob = np.sum(values)
	probabilities = np.divide(values, sum_prob) # Obtain list of probabilities by value/sum_prob
	probabilities = np.round(probabilities, decimals = 5)
	# If just want Pr of best genotype use:                #Note: Change here for extracting pr of best genotype or all genotypes 
	#best_pr = str(np.amax(probabilities))
	#return best_pr
	# Or if want Pr of all genotypes use:
	list_pr = probabilities.tolist()
	list_pr = [str(i) for i in list_pr]
	return list_pr


def append_line(file, line):
	with open(file, "a") as target_file:
		target_file.write(line)


###### Code  #######

# Specify input files from command line arguments:
glf_file = sys.argv[1]
output_file = sys.argv[2]

# Calculate heterozygosity/input heterozygosity value:
heterozygosity = 0.0014

# Write posterior probabilities to file:
with open(glf_file, "r") as file:
	glf = csv.reader(file, delimiter = "\t")
	priors = prior_calculator(heterozygosity)   # Comment out if using uniform prior
	for line in glf:
		gls = float_gls(line)
		if np.sum(gls) == 0:
			continue
		else:
			post_probabilities = gl_calculator(gls, priors) # Add in priors if not uniform; if uniform just provide gls (default priors = 0.1)
			bed_pos = str(int(line[1])-1) # Add in 0 based co-ordinate to convert to bed file
			out_line = [line[0], bed_pos, line[1]] + post_probabilities  # Note: Change gl_calculator function output to get probabilties for all genotypes or for just the bestgenotype
			out_line = "\t".join(out_line) # Separate values by tabs
			out_line = out_line + "\n" 
			append_line(output_file, out_line)  # Add line to out file

# Time taken to run script:
end = time.time()

print(end-start)

