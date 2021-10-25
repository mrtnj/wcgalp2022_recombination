#!/bin/bash

## Run simulation to generate founder genomes for the chicken

set -eu

if [ ! -d simulations/chicken_genome ]; then
  mkdir -p simulations/chicken_genome
fi


Rscript R/simulate_chicken_genome_founders.R

