import argparse
import pandas as pd
##  Here's a script to calculate the average recombination rate in a specified chunk of a genome, given the genetic map

def weighted_average(rec_rates):
## Calculate the weighted average of recombination rates across the analysis window
	print('weighted average')
	rec_rates['length'] = rec_rates.END - rec_rates.START + 1
	r_rate = sum((rec_rates.rate * rec_rates.length))/ sum( rec_rates.length )
	return(r_rate)
#	print(rec_rates)
#	print(r_rate)
	
def get_rec_rate(row, rec_map):
#	get the segments of the recombination map that overlap with the start and end of the genomic window that you are interested in
	if row.BIN_START == 0: 
		recStart = rec_map[(row.BIN_START+1 >= rec_map['START']) & (row.BIN_START < rec_map['END'])]
	else:
		recStart = rec_map[(row.BIN_START+1 >= rec_map['START']) & (row.BIN_START < rec_map['END'])]
		
	recEnd = rec_map[(row.BIN_END >= rec_map['START']) & (row.BIN_END <= rec_map['END'])]

# If there is no recEnd then your analysis window falls off the recombination map and the rate can't be calculated
	if len(recEnd ) == 0:
		return None

# If there a single recombination rate segment coveres the analysis window, just report that.

	if recStart.index == recEnd.index:
		r_rate = float(recStart['rate'])

# If the analysis window spans >1 recombination rate segment
	elif recEnd.index == recStart.index + 1:
		
# Get the recombination rate segments
		rec_rates =  rec_map.iloc[ recStart.index[0] :recEnd.index[0]+1 ].copy()
#		print( rec_rates)
# There is no need to calculate the weighted average if the different segments of the map have the same rate
		if len (set(rec_rates.rate) ) == 1:
			r_rate = list(rec_rates.rate)[0] 
		else:
#			length = int( row.BIN_END - row.BIN_START + 1 )
#			lower = recStart.END - row.BIN_START + 1
#			upper = row.BIN_END - recEnd.START + 1

			rec_rates.iloc[0, 0] = row.BIN_START
			rec_rates.iloc[-1, 2] = row.BIN_END

			r_rate = weighted_average(rec_rates)
#			rec_map.iloc[ recStart.index[0] :recEnd.index[0]+1 ]

#			print('\n', length, lower, upper)
			print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
	
	return r_rate


def main():
	
## Start by defining the command line arguments, each one should be self-explanatory
	parser = argparse.ArgumentParser(description="Adds the recombination rates from a specified recombination map to the output from VCFtools. In this case, I use the genetic maps spat out from msprime and stdpopsim")

	parser.add_argument("-i", "--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The Fst output file")
	parser.add_argument("--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "Give the output file some informative name")
	parser.add_argument("--g_map", 
		required = True,
		dest = "g_map",
		type = str, 
		help = "specify the recombination map")
		
	args = parser.parse_args()
	
	if args.input  == args.output:
		print('Use a different name for your input and your output')
		return
	Gmap = pd.read_csv(args.g_map, names = ['START','rate'], sep = ' ')
	Gmap['START'] = Gmap['START'] + 1
	Gmap['END'] = (Gmap['START'] -1).shift(-1)
	Fst = pd.read_csv(args.input, sep = ' ', names = ['CHROM','BIN_START','BIN_END','VARS','T1','T2','WEIGHTED_FST'])
	print(Fst)

	rates = []
	for index, row  in Fst.iterrows():
#		if index != 2255: continue
#		print(index)
		rec_rate = get_rec_rate(row, Gmap) 
		rates.append( rec_rate )
	Fst['rec'] = rates
	Fst.to_csv(args.output, index = False)

main()
