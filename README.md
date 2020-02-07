# Recombination rate variation and population genetic summary statistics

A repository containing the scripts and analysis results behind the manuscript "Variation in recombination affects detection of Fst outliers under neutrality". The manuscript is currently in review and we will deposit it on BioRXiv in the coming days.

There are two main parts to this repository:

  1. The sampling distribution of Fst calculated in analysis windows of fixed physical size for coalescent simulations of an island model - these simulations went into Figure 1 in the MS. There are additional simulation results 

  2. Analysis of simulated and empirical *Drosophila melanogaster* datasets. These results went into Figure 2 of the paper. 

Within each of the directories, you'll find all the code to run the simualtions and data analysis that went into the paper, plus a few extra bits and pieces. 

## Motivation

This study was motivated by the following pattern:
![](Recombination_Fst/writeUp/RecRateVariation.png)

The figure shows how variable recombination rates are across the genomes of three animal species that have been subject to genome scan studies. All three organsisms exhibit wide variation in recombination rate, over at least an order of magnitude. What we wanted to know was, does such variation in recombination rate affect our the outcome of genome scans? 

We think that it does, and specifically that population genetic summary statistics calculated in regions of low recombination may have a different sampling distribution to statistics calculated in more highly recombining regions. 
