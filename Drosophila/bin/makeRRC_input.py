## Make the Reinhardt et al (2014) recombination intervals
import sys


## Remember to run the RRC perl script in the directory that contains the Comeron/ subdirectory 
##2L:7500..9400
def main():
	if len( sys.argv ) < 3:
		print('USAGE:\n\tpython FstFile.txt output.txt')
		sys.exit()

	chrom = sys.argv[1].split('/')[-1].split('.')[1].split('r')[1]

	output = open(sys.argv[2], 'w')

	for i in open(sys.argv[1]):
		x = i.strip().split()
		output.write(chrom + ':' + x[1] + '..' + x[2] + '\n')
	output.close() 
main()
