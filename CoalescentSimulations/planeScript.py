import msprime
import argparse
import pandas as pd
import numpy as np
import allel
from sklearn.utils import shuffle
from multiprocessing import Pool

# Here's a script to examine properties of the distribution of windowed Fst


## 20 DEME RING MODEL
def migration_simulation_donut( r_rate ):

    n_demes = 20
    m = 0.1/3.

##  Allocate the initial sample
##  Sample 20 individuals from two demes
    population_configurations_sample_1 = [
        msprime.PopulationConfiguration(sample_size=0),
        msprime.PopulationConfiguration(sample_size=20)] 
    population_configurations_sample_2 = [
        msprime.PopulationConfiguration(sample_size=0),
        msprime.PopulationConfiguration(sample_size=20)] 

    population_configurations_empty = [msprime.PopulationConfiguration(sample_size=0) for i in range(n_demes - len(population_configurations_sample_1) * 2 )] 
    population_configurations = population_configurations_sample_1 + population_configurations_empty + population_configurations_sample_2
    
    migration_matrix = []
    for i in range(n_demes):
        temp = [0] * n_demes
        if i ==0:
             temp[i+1] = m
             temp[n_demes-1] = m
        elif i == n_demes-1:
            temp[i-1] = m
            temp[0] = m
        else:
            temp[i-1] = m
            temp[i+1] = m
        migration_matrix.append( temp )
    
    # We pass these values to the simulate function, and ask it
    # to run the required number of replicates.
    forest = msprime.simulate(Ne=1000,
        population_configurations=population_configurations,
        length = 10000,
        migration_matrix=migration_matrix,
        recombination_rate = r_rate)
    return( forest )

## 20 DEME CLINE MODEL
def migration_simulation_cline( r_rate ):

    n_demes = 20
    m = 0.1/3.

##  Allocate the initial sample
##  Sample 20 individuals from two demes
    population_configurations_sample_1 = [
        msprime.PopulationConfiguration(sample_size=0),
        msprime.PopulationConfiguration(sample_size=20)] 
    population_configurations_sample_2 = [
        msprime.PopulationConfiguration(sample_size=0),
        msprime.PopulationConfiguration(sample_size=20)] 

    population_configurations_empty = [msprime.PopulationConfiguration(sample_size=0) for i in range(n_demes - len(population_configurations_sample_1) * 2 )] 
    population_configurations = population_configurations_sample_1 + population_configurations_empty + population_configurations_sample_2
    
# Now we set up the migration matrix. Since this is a symmetric
# island model, we have the same rate of migration between all
# pairs of subpopulations. Diagonal elements must be zero.
	
    migration_matrix = []
    for i in range(n_demes):
        temp = [0] * n_demes
        if i ==0:
             temp[i+1] = m
        elif i == n_demes-1:
            temp[i-1] = m
        else:
            temp[i-1] = m
            temp[i+1] = m
        migration_matrix.append( temp )
# We pass these values to the simulate function, and ask it
# to run the required number of replicates.
    forest = msprime.simulate(Ne=500,
        population_configurations=population_configurations,
        length = 10000,
        migration_matrix=migration_matrix,
        recombination_rate = r_rate)
    return( forest )
    
## 100 DEME ISLAND MODEL
def migration_simulation_Island( r_rate ):

    n_demes = 100
    m = 2*1/10000
	
    # Allocate the initial sample
    population_configurations_samples = [
        msprime.PopulationConfiguration(sample_size=20),
        msprime.PopulationConfiguration(sample_size=20)] 
    population_configurations_empty = [msprime.PopulationConfiguration(sample_size=0) for i in range(n_demes - len(population_configurations_samples) )] 
    population_configurations = population_configurations_samples + population_configurations_empty
    
   # Now we set up the migration matrix. Since this is a symmetric
   # island model, we have the same rate of migration between all
   # pairs of subpopulations. Diagonal elements must be zero.
	
    migration_matrix = []
    for i in range(n_demes):
        temp = [m] * n_demes
        temp[i] = 0
        migration_matrix.append( temp )
# We pass these values to the simulate function, and ask it
# to run the required number of replicates.
    forest = msprime.simulate(Ne=100,
        population_configurations=population_configurations,
        length = 10000,
        migration_matrix=migration_matrix,
        recombination_rate = r_rate)
    return( forest )

## 2 DEMES
def migration_simulation_2patch( r_rate ):
	d = 2
	m = 1/10000

# Allocate the initial sample (20 from each deme)
	population_configurations = [
		msprime.PopulationConfiguration(sample_size=20),
		msprime.PopulationConfiguration(sample_size=20)]
		
# Now we set up the migration matrix. 
	migration_matrix = [
		[0, m],
		[m, 0] ]
# We pass these values to the simulate function, and ask it
# to run the required number of replicates.
	forest = msprime.simulate(Ne=10000,
		population_configurations=population_configurations,
		length = 10000,
		migration_matrix=migration_matrix,
		recombination_rate = r_rate)
	return( forest )

	
## Now write a function to run the simulation many times for each 
## value of Fst
def RecombinationRepper( pooled_args ): #  provide r_rate, model_function, reps, samples 
	
	r_rate = pooled_args[0]
	model_function = pooled_args[1]
	reps = pooled_args[2]
	samples = pooled_args[3]

  
	mean_Fst_dists = []
	var_Fst_dists = []
	mean_SE_dists = []
	mean_SE_dists_shuf = []
	tree_counts = []
	mean_Dxy_dists = []
	mean_Tajima_dists = []
	mean_diversity_dists = []
	mean_H12_dists = []

	for t in range( reps ):
		print(t)

		new_tree  = migration_simulation_2patch( r_rate )
		# Add mutations to the tree
		count = 0
		for r in new_tree.trees():
			count += 1
#		print(t, count)
		tree_counts.append( count )
		
		new_tree_dist_Fst = []
		new_tree_dist_SE = []
		new_tree_dist_SE_shuf = []
		new_tree_dist_Dxy = []
		new_tree_dist_diversity = []
		new_tree_dist_Tajima = []
		new_tree_dist_H12 = []


		for i in range( samples ): ## Repeat 100 times
			# Add mutations to the tree
			muts = 0
			while muts == 0:
				mutated_tree = msprime.mutate(new_tree, 1.25e-7)
				muts = len( [ v for v  in mutated_tree.variants() ] )
			# Get the genotype matrix, ready for using sci-kit.allel
			msprime_genotype_matrix = mutated_tree.genotype_matrix()
			# Convert msprime's haplotype matrix into genotypes by randomly merging chromosomes
			haplotype_array = allel.HaplotypeArray( msprime_genotype_matrix )

			genotype_array = haplotype_array.to_genotypes(ploidy=2)

			shuffled_genotypes = shuffle(genotype_array, random_state=0)

			ac1 = haplotype_array.count_alleles(subpop=[s for s in range(0,100)])
			ac2 = haplotype_array.count_alleles(subpop=[s for s in range(100,200)])

## Calculate Tajima's D
			Tajimas_D  = allel.tajima_d(ac1) 

## Calculate Dxy
			dxy = sum( allel.mean_pairwise_difference_between(ac1, ac2) ) / 10000.

## Calculate Garud's H statistics for the population
	## Grab the haplotypes for 400SNPs from deme 1
			hapslice = haplotype_array[:400,0:100]
			H_vector  = allel.garud_h( hapslice ) 

## Calculate Diversity
			pi = sum(allel.mean_pairwise_difference(ac1) ) / 10000.
			
			subpopulations = [[p for p in range(0,50)], [z for z in range(50,100)]]
			mean_fst = allel.average_weir_cockerham_fst(genotype_array, blen = 100, subpops=subpopulations)
			mean_fst_shuf = allel.average_weir_cockerham_fst(shuffled_genotypes, blen = 100, subpops=subpopulations)
			
			new_tree_dist_Fst.append( mean_fst[0] )
			new_tree_dist_SE.append( mean_fst[1] )
			new_tree_dist_SE_shuf.append( mean_fst_shuf[1] )
			new_tree_dist_Tajima.append( Tajimas_D )
			new_tree_dist_Dxy.append( dxy )
			new_tree_dist_H12.append( H_vector[1] )
			new_tree_dist_diversity.append(  pi )

		mean_Fst_dists.append( np.mean( new_tree_dist_Fst ) )

		var_Fst_dists.append( np.sqrt( np.var( new_tree_dist_Fst ) ) )

		mean_SE_dists.append( np.mean( new_tree_dist_SE ) )

		mean_SE_dists_shuf.append( np.mean( new_tree_dist_SE_shuf ) )
		
		mean_Dxy_dists.append( np.mean( new_tree_dist_Dxy ) )

		mean_Tajima_dists.append( np.mean( new_tree_dist_Tajima ) )

		mean_H12_dists.append( np.mean( new_tree_dist_H12 ) )

		mean_diversity_dists.append( np.mean( new_tree_dist_diversity ) )

	return [r_rate, mean_Fst_dists, mean_SE_dists, mean_SE_dists_shuf, var_Fst_dists, tree_counts, mean_Dxy_dists, mean_Tajima_dists, mean_diversity_dists,mean_H12_dists ]

def main():
	
## Start by defining the command line arguments, each one should be self-explanatory
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("-m", "--model", 
		required = True,
		dest = "model",
		type =str, 
		help = "The name of the model that you want to simulate. [island, 2patch, cline, ring]")
	parser.add_argument("--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "Give the output file some informative name")
	parser.add_argument("--reps", 
		required = False,
		dest = "reps",
		type = int, 
		help = "Give the number of simulation replicates that you want to perform [1000]",
		default = 1000)
	parser.add_argument("--samples", 
		required = False,
		dest = "samples",
		type = int, 
		help = "How many times you want to add mutations to the trees that you simulate? [100]",
		default = 100)
	parser.add_argument("--nproc", 
		required = False,
		dest = "nproc",
		type = int, 
		help = "How many threads do you want to use? [1]",
		default = 1)
		
	args = parser.parse_args()

## Make sure the input model matches one of the ones there is a function for
	if args.model.upper() not in ['ISLAND', '2PATCH', 'CLINE', 'RING']:
		print('please give your model from one of the following:' +'\n' + '\n'.join(['ISLAND', '2PATCH', 'CLINE', 'RING'] +'\nAll three should give Fst of around 0.11'))
		return
	if args.model.upper() == 'ISLAND':
		sim_function = migration_simulation_Island
	elif args.model.upper() == '2PATCH':
		sim_function = migration_simulation_2patch
	elif args.model.upper() == 'CLINE':
		sim_function = migration_simulation_cline
	elif args.model.upper() == 'RING':
		sim_function = migration_simulation_donut
		
# Inititialise a pool that you can map arguments to
	p = Pool( args.nproc )
	r_list = []

# Store the arguments for the main Recombination function in a list
	pool_arguments = [[r, sim_function, args.reps, args.samples]  for r in [ 0 , 1e-8, 1e-7, 1e-6 ]]
	
# Here's the meat, we map the results to the pool, i.e. run RecombinationRepper in parallel over the specified number of threads
	simulation_results = p.map(RecombinationRepper, pool_arguments)


	for r in simulation_results:
#	r_rate, mean_Fst_dists, mean_SE_dists, mean_SE_dists_shuf, var_Fst_dists, tree_counts, d_xy, Tajima's D, pi, H12
		fst_dict = {'rec': r[0], 'weighted_Fst': r[1], 'mean_SE': r[2], 'mean_SE_shuf': r[3], 'variance': r[4], 'num_trees': r[5], 'd_xy' : r[6], 'Tajima_D' : r[7], 'diversity' : r[8], 'H12' : r[9]} 
		r_list.append( pd.DataFrame( fst_dict ) )


#	pd.concat(r_dict).to_csv('Fst_density.2deme.csv', index = False) 

#	pd.concat(r_dict).to_csv('Fst_density.cline.csv', index = False) 

	pd.concat(r_list).to_csv(args.output, index = False) 
#	pd.concat(r_dict).to_csv('Fst_density.donutModel.csv', index = False) 
	return
	
main()

