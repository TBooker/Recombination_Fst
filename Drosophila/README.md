# Analysis of simulated and empirical *Drosophila melanogaster* datasets

This directory contains the scripts and analysis files that were used to generate the results shown in Figure 2
![](../writeUp/threePanelFigure.png)

The script [bin/plotDrosophila_plusSimulations.R](bin/plotDrosophila_plusSimulations.R) contains the code to do all the statistical analyses and the plotting for this part of the paper. 

Below I demonstrate the code used to generate the simulated data and to analyse the empirical data.


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
with open("drosophilaSimulated.chr2L.vcf", "w") as vcf_file:
	ts.write_vcf(vcf_file)


```

I then calculated Fst in 10,000bp analysis windows from the simulated VCFs. To make the results comparable with the empirical *Drosophila* data, I used the haploid method for calculating Fst from Weir's book "Genetic Data Analysis" (1990; pp 145-148).


### Here's a pipeline to run and analsye the simulation data
I use the script (read_VCF.py)[bin/read_VCF.py] to generate an input file from the VCF that is the same format as the one I get from the empirical data (see below). 

Run the simulation script for each of the normal autosomes (sorry dot chromosome!). Here I make use of GNU parallel, but you don't have to parallelise it this way if you don't want to.

```
parallel "python3.6 bin/run_2Pop_simulations.py chr{}" ::: 2L 2R 3L 3R
```

Now, convert the VCF file for each autosome into an input format that we can get from the Reinhardt et al data...

```
parallel "python bin/read_VCF.py --VCF drosophilaSimulated.chr{}.vcf -o drosophilaSimulated.chr{}.txt" ::: 2L 2R 3L 3R
```
This produces a file for each chromosome with the number of alleles of each type for each biallelic SNP in the simulated population

Now we calculate Fst for each site and get the weighted average (using the ratio of averages approach outlined in Weir's textbook) in analysis windows of 10,000bp.

```
parallel "python bin/calculateFst.py --freqs drosophilaSimulated.chr{}.txt --MAF 0.05 > drosophilaSimulated.chr{}.fst" ::: 2L 2R 3L 3R
```

Now grab the recombination rates from the genetic map and add that to the analysis files.

```
parallel "python bin/add_recombination_rates.py --input drosophilaSimulated.chr{}.fst --output drosophilaSimulated.chr{}.fst.csv --g_map ../chr{}.map.txt" ::: 2L 2R 3L 3R
```
This script just looks up the recombination rates from the map given in **stdpopsim**. It let's you know if has had to calculate a weighted average recombination rate for a particular analysis window.

Now just compress the VCF files abd remove the intermediate files to save a little disc space
```
gzip *vcf
rm *txt
rm *fst
```
This process should result in the files:

	drosophilaSimulated.CHROM.vcf.gz
	drosophilaSimulated.CHROM.fst.csv
Except that there will be one per autosome.

# Empirical data
