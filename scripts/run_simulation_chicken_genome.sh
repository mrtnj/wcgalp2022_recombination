#!/bin/bash

## Run simulations with approximation of chicken karyotype

set -eu

if [ ! -d simulations/chicken_genome ]; then
  mkdir -p simulations/chicken_genome
fi


## Additive


if [ ! -d simulations/chicken_genome/additive/ ]; then
  mkdir simulations/chicken_genome/additive/
fi

for TOT_QTL in 100 1000 10000; do

  for REP in {1..20}; do
    
    if [ ! -d simulations/chicken_genome/totqtl${TOT_QTL} ]; then
      mkdir simulations/chicken_genome/totqtl${TOT_QTL}
    fi

    Rscript R/simulate_chicken_genome.R \
      ${TOT_QTL} \
      0 \
      0 \
      simulations/chromosome_length/additive/totqtl${TOT_QTL}/populations_${REP}.Rds \
      simulations/chromosome_length/additive/totqtl${TOT_QTL}/results_${REP}.Rds
    
  done
  
done



## Dominance

if [ ! -d simulations/chicken_genome/dominance/ ]; then
  mkdir simulations/chicken_genome/dominance/
fi

for TOT_QTL in 100 1000 10000; do

  for REP in {1..20}; do
    
    if [ ! -d simulations/chicken_genome/totqtl${TOT_QTL} ]; then
      mkdir simulations/chicken_genome/totqtl${TOT_QTL}
    fi

    Rscript R/simulate_chicken_genome.R \
      ${TOT_QTL} \
      0.2 \
      0.1 \
      simulations/chromosome_length/additive/totqtl${TOT_QTL}/populations_${REP}.Rds \
      simulations/chromosome_length/additive/totqtl${TOT_QTL}/results_${REP}.Rds
    
  done
  
done
