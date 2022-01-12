# The effect of recombination rate on genomic selection in simulation

This repository contains code underlying the conference contribution:

M. Johnsson, A.M. Johansson, T. Klingström and E. Rius-Vilarrasa (2022)
"The effect of recombination rate on genomic selection in simulation" 
Submitted to WCGALP 2022.


## Simulations with variable recombination rate

The main contents is a simulation using `AlphaSimR` with real published
linkage maps from cattle and chickens.

* `R/simulation_functions.R` -- helper functions for the AlphaSimR simulation

* `R/simulate_chromosome_length.R` -- simulation with variable chromosome
length but uniform recombination rate

* `R/simulate_karyotype.R` -- simulation with real chromsome lengths
but uniform recombination rate

* `R/simulate_real_map.R` -- simulation with real linkage maps


## Preparataion of maps

* `R/make_cattle_genome_table.R` and `R/make_chicken_genome_table.R` process
linkage maps and physical maps for simulation

* `R/plot_map_difference.R` illustrates the adjustment of simulated maps to
real maps


## Linkage equilibrium simulation

For comparison, we run a simulation with similar genetic parameters but
genotypes are drawn independently each generation.

* `R/simulation_functions_linkage_equilibrium.R` -- helper functions for
LE simulation

* `R/simulate_linkage_equilibrium.R` -- run LE simulation


## Bash scripts to run

The `scripts` directory contains bash scripts to run replicates of the
simulations with different parameters.


## Summaries of results

* `R/summarise_results.R` and `R/summarise_results_linkage_equilibrium.R` 
summarises simulation runs and puts statistics in the `outputs` directory.

* `R/plot_results.R` contains plotting code