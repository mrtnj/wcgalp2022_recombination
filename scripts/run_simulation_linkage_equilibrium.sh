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

for START in 1 10; do

  for TOT_QTL in 100 1000 10000; do

    for REP in {1..50}; do
    
      if [ ! -d simulations/linkage_equilibrium/dominance/totqtl${TOT_QTL}_start${START} ]; then
        mkdir simulations/linkage_equilibrium/dominance/totqtl${TOT_QTL}_start${START}
      fi

      Rscript R/simulate_linkage_equilibrium.R \
        simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}/simparam_${REP}.Rds \
        simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}/genotypes_${REP}.Rds \
        ${START} \
        directional \
        simulations/linkage_equilibrium/dominance/totqtl${TOT_QTL}_start${START}/results_${REP}.Rds \
        simulations/linkage_equilibrium/dominance/totqtl${TOT_QTL}_start${START}/generations_${REP}.Rds
        
      done
    
  done
  
done

