#!/bin/bash

## Run simulation to generate founder genomes for cattle and chicken

set -eu

if [ ! -d simulations/chicken_genome ]; then
  mkdir -p simulations/chicken_genome
fi


Rscript R/simulate_karyotype_founders.R \
  annotation/chicken_genome_table.txt \
  simulations/chicken_genome/chicken_genome_founders.Rds


if [ ! -d simulations/cattle_genome ]; then
  mkdir -p simulations/cattle_genome
fi


Rscript R/simulate_karyotype_founders.R \
  annotation/cattle_genome_table.txt \
  simulations/cattle_genome/cattle_genome_founders.Rds
  