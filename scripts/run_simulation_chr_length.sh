#!/bin/bash

## Run simulations with varying chromosome lengths

set -eu

if [ ! -d simulations/chromosome_length ]; then
  mkdir -p simulations/chromosome_length
fi


## Additive


if [ ! -d simulations/chromosome_length/additive/ ]; then
  mkdir simulations/chromosome_length/additive/
fi

for NQTL in 10 100 1000; do

  for GEN_LENGTH in 0.5 1 2; do

    for REP in {1..20}; do
    
      if [ ! -d simulations/chromosome_length/chrlen${GEN_LENGTH}_nqtl${NQTL} ]; then
        mkdir simulations/chromosome_length/chrlen${GEN_LENGTH}_nqtl${NQTL}
      fi

      Rscript R/simulate_chromosome_length.R \
        ${GEN_LENGTH} \
        ${NQTL} \
        0 \
        0 \
        simulations/chromosome_length/additive/chrlen${GEN_LENGTH}_nqtl${NQTL}/populations_${REP}.Rds \
        simulations/chromosome_length/additive/chrlen${GEN_LENGTH}_nqtl${NQTL}/results_${REP}.Rds
    
    done
  
  done
  
done



## Dominance

if [ ! -d simulations/chromosome_length/dominance/ ]; then
  mkdir simulations/chromosome_length/dominance/
fi

for NQTL in 10 100 1000; do

  for GEN_LENGTH in 0.5 1 2; do

    for REP in {1..20}; do
    
      if [ ! -d simulations/chromosome_length/chrlen${GEN_LENGTH}_nqtl${NQTL} ]; then
        mkdir simulations/chromosome_length/chrlen${GEN_LENGTH}_nqtl${NQTL}
      fi
  
      Rscript R/simulate_chromosome_length.R \
        ${GEN_LENGTH} \
        ${NQTL} \
        0.2 \
        0.1 \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_nqtl${NQTL}/populations_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_nqtl${NQTL}/results_${REP}.Rds
    
    done
  
  done
  
done
