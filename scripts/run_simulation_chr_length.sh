#!/bin/bash

## Run simulations with varying chromosome lengths

set -eu

if [ ! -d simulations/chromosome_length ]; then
  mkdir -p simulations/chromosome_length
fi




## Dominance

if [ ! -d simulations/chromosome_length/dominance/ ]; then
  mkdir simulations/chromosome_length/dominance/
fi

for NQTL in 100 1000 10000; do

  for GEN_LENGTH in 0.5 1 2; do

    for REP in {1..50}; do
    
      if [ ! -d simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL} ]; then
        mkdir simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}
      fi
  
      Rscript R/simulate_chromosome_length.R \
        ${GEN_LENGTH} \
        ${NQTL} \
        0.2 \
        0.1 \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/genotypes_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/phenotypes_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/results_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/simparam_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/genotypes_gs_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/phenotypes_gs_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/results_gs_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/generation1_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/generation10_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/generation20_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/model_${REP}.Rds
    
    done
  
  done
  
done
