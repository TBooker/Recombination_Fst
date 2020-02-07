import sys, argparse
from collections import Counter
import numpy as np

def count_alleles(alleleString):
	As = alleleString.count('A')
	Cs = alleleString.count('C')
	Ts = alleleString.count('T')
	Gs = alleleString.count('G')
	return  np.array( [ As, Cs, Ts, Gs ] )
	
def get_allele_counts(pop1, pop2):
	for p1_raw, p2_raw in zip( open(pop1), open(pop2) ):
		p1 = p1_raw.strip().split()
		p2 = p2_raw.strip().split()
## If there is no polymorphism, ignore the site and move to the next
		
		p1_counts = count_alleles(p1[3])
		p2_counts = count_alleles(p2[3])

## If there are 3 zeroes in the summed array of the two populations' allele count arrays
#		print(p1_counts + p2_counts, np.count_nonzero( p1_counts + p2_counts ) == 1)
		if np.count_nonzero( p1_counts + p2_counts ) == 1:
#			print(p1_counts + p2_counts, '!')
			yield(None, None, None, None,None)
			
		else:
#		zeros = np.count_nonzero(arr==0)
#		p1_counts  = ':'.join( map(str,
			yield(p1[0],
					str(np.sum(p1_counts)),
					str(np.sum(p2_counts)),
					':'.join( map(str, p1_counts ) ),
					':'.join( map(str, p2_counts ) ))

def main():
	# Define the command line arguments
	parser = argparse.ArgumentParser(description="This script will run through the Reinhardt et al data and get the allele counts, reporting allele counts for variable site")
	
	parser.add_argument("-p1", "--pop1", 
		required = True,
		dest = "pop1",
		type =str, 
		help = "Specify pileup file for population 1")

	parser.add_argument("-p2", "--pop2", 
		required = True,
		dest = "pop2",
		type =str, 
		help = "Specify pileup file for population 2")

	parser.add_argument("-o", "--output", 
		required = True,
		dest = "output",
		type = str, 
		help = "Give the name for the output file")
		
	args = parser.parse_args()

	output = open(args.output, 'w')

	count =0
	for pos, p1_depth, p2_depth, p1_cownts, p2_cownts in get_allele_counts(args.pop1, args.pop2):
		count +=1
		if p1_cownts == None:
			continue
#		print(pos, p1_depth, p2_depth, p1_cownts, p2_cownts)
		output.write('\t'.join([ pos, p1_depth, p2_depth, p1_cownts, p2_cownts ]) + '\n' )
#		print(count, '\n')
#		if count == 1000:
#			break
	output.close()

# We want

if __name__ == '__main__':
    main()
