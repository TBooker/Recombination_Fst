import sys, argparse
from collections import Counter
import numpy as np


def count_alleles(alleleString):
	As = alleleString.count('A')
	Cs = alleleString.count('C')
	Ts = alleleString.count('T')
	Gs = alleleString.count('G')
	return  np.array( [ As, Cs, Ts, Gs ] )
	
def get_allele_counts(VCF):
	for line in open(VCF):
		x = line.strip().split()
		if line.startswith('#'):
			print(line)
			if line.startswith('#CHROM'):
				n_samples = len( x[9:] )
			continue
		p1 = x[9:9+int(n_samples/2)] 
		p2 = x[9+int(n_samples/2):]

		p1_counts_p = p1.count('1')
		p1_counts_q = p1.count('0')

		p2_counts_p = p2.count('1')
		p2_counts_q = p2.count('0')
		
		
		pos = x[1]
		
		p1_cownts = [p1_counts_p, p1_counts_q, 0, 0]
		p2_cownts = [p2_counts_p, p2_counts_q, 0, 0]
		yield(pos,
				str(np.sum(p1_cownts)),
				str(np.sum(p2_cownts)),
				':'.join( map(str, p1_cownts ) ),
				':'.join( map(str, p2_cownts ) ))
	
def main():
	# Define the command line arguments
	parser = argparse.ArgumentParser(description="This script will run through the stdpopsim VCF data, reporting allele counts for variable sites")
	
	parser.add_argument("-v", "--VCF", 
		required = True,
		dest = "VCF",
		type =str, 
		help = "Specify the VCF file")

	parser.add_argument("-o", "--output", 
		required = True,
		dest = "output",
		type = str, 
		help = "Give the name for the output file")
		
	args = parser.parse_args()

	output = open(args.output, 'w')

	count =0
	
	for pos, p1_depth, p2_depth, p1_cownts, p2_cownts in  get_allele_counts(args.VCF):
		count +=1
		print(pos, p1_depth, p2_depth, p1_cownts, p2_cownts)
		output.write('\t'.join([ pos, p1_depth, p2_depth, p1_cownts, p2_cownts ]) + '\n' )
	output.close()
		
if __name__ == '__main__':
    main()
    
    
    
    
    
    
