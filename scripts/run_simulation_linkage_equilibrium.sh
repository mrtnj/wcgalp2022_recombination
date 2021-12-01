#!/bin/bash

## Run simulations with linkage equilibrium

set -eu


## Dominance

if [ ! -d simulations/linkage_equilibrium ]; then
  mkdir simulations/linkage_equilibrium
fi

if [ ! -d simulations/linkage_equilibrium/dominance/ ]; then
  mkdir simulations/linkage_equilibrium/dominance/
fi

for TOT_QTL in 100 1000 10000; do

  for REP in {1..20}; do
    
    if [ ! -d simulations/linkage_equilibrium/dominance/totqtl${TOT_QTL} ]; then
      mkdir simulations/linkage_equilibrium/dominance/totqtl${TOT_QTL}
    fi

    Rscript R/simulate_linkage_equilibrium.R \
      simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}/simparam_${REP}.Rds \
      simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}/genotypes_${REP}.Rds \
      1 \
      directional \
      simulations/linkage_equilibrium/dominance/totqtl${TOT_QTL}/results.Rds
    
  done
  
done

