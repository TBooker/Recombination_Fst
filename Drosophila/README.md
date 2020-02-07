# Analysis of simulated and empirical *Drosophila melanogaster* datasets

This directory contains the scripts and analysis files that were used to generate the results shown in Figure 2

# Simulated data

I used *stdpopsim* to model the entire *D. melanogaster* genome incorporating the recombination rate variation estimated by Comeron et al (2012).

I used the script [run_2Pop_simulations.py](bin/run_2Pop_simulations.py) to run the genome simualtions. Using *stdpopsim* was really fun and I'll definitely make use of it in the future, look how easy it is to simulate chromosome 2L:

```python
import stdpopsim,sys

## I would like to simulate Drosophila melanogaster please
species = stdpopsim.get_species("DroMel")
## I have specified the desired chromosome arm at the command line, let's 
contig = species.get_contig("chr2L", genetic_map = "ComeronCrossover_dm6")

## For testing, it is good to model a lil chunk of chromosome
#contig = species.get_contig("chr2L", length_multiplier = 0.10)

## You can grab the genetic map out of the simulations using:
#for p, r in zip( contig.recombination_map.get_positions() , contig.recombination_map.get_rates() ):
#	print( p , r )

Ne = species.population_size/10

model = stdpopsim.IsolationWithMigration(2*Ne, Ne, Ne, Ne, 1.5/Ne, 1.5/Ne)

print("NA", "N1", "N2", "T", "M12", "M21")
print( 2*Ne, Ne, Ne, Ne, 1.5/Ne, 1.5/Ne)

## I want to simulate 20 samples from each population
samples = model.get_samples(20,20)

## I will simulate using msprime
engine = stdpopsim.get_default_engine()

print("running simulation")
ts = engine.simulate(model, contig, samples)

## Save the simulated data to a VCF
with open("drosophilaSimulated."chr2L".vcf", "w") as vcf_file:
	ts.write_vcf(vcf_file)


```

# Empirical data
