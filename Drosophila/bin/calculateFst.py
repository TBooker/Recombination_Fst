import sys, argparse
import numpy as np
import pandas as pd

def isMultiAlellic( site ):
	alleles =  np.array( list( map(int, site[3].split(':'))) ) + np.array( list( map(int, site[3].split(':'))) ) 
	if list(alleles).count(0) == 1:
		return True
	else:
		return False

def alleleFreqs( alleleCounts ):

	majorAllele = alleleCounts.index(max( alleleCounts ) )

	return( max( alleleCounts ), sum(alleleCounts) -max( alleleCounts ), majorAllele )

def Fst_haploid( demeCounts, p1_major, p2_major):


	demeFreqs = [max(deme)/sum(deme) for deme in demeCounts]



	r = len(demeFreqs) # Number of demes
	sum_n_i_Sq = sum( [ sum(deme)**2 for deme in demeCounts] )

	sum_n_i = sum( [ sum(deme) for deme in demeCounts] )
	
	nc = (1/(r-1)) * ( sum([sum(deme) for deme in demeCounts]) - (sum_n_i_Sq/sum_n_i) )
	
	n_bar = sum_n_i/r 

	if p1_major == p2_major:
		p_a_squig = sum( demeFreqs ) / len( demeFreqs )
	elif p1_major != p2_major:
		temp = list( np.array(demeCounts[0]) + np.array(demeCounts[1]) )
		p_a_squig = max(temp)/ sum(temp)
#		print( p_a_squig)
		

	sumSq_p = [sum(c)*(p - p_a_squig)**2 for p,c in zip(demeFreqs, demeCounts) ]
	
	S_A_sq = (1/((r-1)*n_bar)) * sum(sumSq_p)

	T1 = S_A_sq - (1/(n_bar -1))*( p_a_squig*(1-p_a_squig) - (S_A_sq*(r-1)/r))

	T2 = ((nc -1)/(n_bar -1))*p_a_squig*(1-p_a_squig) + ( 1 + ((r-1)*(n_bar - nc))/(n_bar - 1))*S_A_sq/r
	return p_a_squig, T1, T2, T1/T2


def main():
	# Define the command line arguments
	parser = argparse.ArgumentParser(description="This script will run through the Reinhardt et al data and get the allele counts, reporting allele counts for variable site")
	
	parser.add_argument("--freqs", 
		required = True,
		dest = "freqs",
		type =str, 
		help = "the file with allele frequnecy information")

	parser.add_argument("--empiric", 
		required = False,
		dest = "empiric",
		action = 'store_true',
		help = "Is the data from the Drosophila data?")
		
	parser.add_argument("--MAF", 
		required = False,
		dest = "MAF",
		type = float,
		help = "Apply an MA filter to the data",
		default = 0)
		
	parser.add_argument("--windows", 
		required = False,
		dest = "windows",
		type = int,
		help = "What window size do you want to analyse?",
		default = 10000)
	args = parser.parse_args()

	chrom_dict = {'chr2L': 23011544, 'chr2R': 21146708, 'chr3L': 24543557, 'chr3R': 27905053}
	chromName = args.freqs.split('.')[-2]
	chrom = chrom_dict[ chromName ]
	
	resultList = []

	count = 0 

	for i in open(args.freqs):
## The mean and std.dev of coverage in the Drosophila data are 23 and 8, respectively  
		site = i.strip().split()
		
		if int(site[1]) > 23 + 8 or int(site[1]) < 23 - 8:
			continue 

		p1_alleles = list(map( int, site[3].split(':') ))
		p2_alleles = list(map( int, site[4].split(':') ))

		if isMultiAlellic( site ):
			continue

		count +=1 
		
		p1a1, p1a2, p1_major = alleleFreqs( p1_alleles )
		p2a1, p2a2, p2_major = alleleFreqs( p2_alleles )
		if p1_major != p2_major:
			p1a2, p1a1, p1_major = alleleFreqs( p1_alleles )

#		print(site)
#		print(p1_major, p2_major)
#		print(p2a1, p2a2)
		resultList.append( [int(site[0])]  + list( Fst_haploid([ [ p1a1, p1a2 ], [ p2a1, p2a2 ] ], p1_major, p2_major)) )

#		if count == 1000: break
	Fst = pd.DataFrame( resultList, columns = ['pos','p_bar','T1','T2','theta']) 
	Fst = Fst[1-Fst['p_bar'] > args.MAF]
	
	count = 0
	for w in range(0,chrom, args.windows):
		count +=1
		start = w
		end = w + args.windows - 1
		variants = Fst[ (Fst['pos'] >= start) & (Fst['pos'] <= end)]
		nvar = len(variants) 
		print(chromName, start, end, nvar, variants.T1.sum(), variants.T2.sum(),  variants.T1.sum()/ variants.T2.sum())

#		print(start, end)
#		if count == 100:
#			break
		
		
if __name__ == '__main__':
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
