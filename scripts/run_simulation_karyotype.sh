#!/bin/bash

## Run simulations with approximation of real karyotypes

set -eu



## Chicken genome


## Dominance

if [ ! -d simulations/chicken_genome/karyotype/ ]; then
  mkdir simulations/chicken_genome/karyotype/
fi

if [ ! -d simulations/chicken_genome/karyotype/dominance/ ]; then
  mkdir simulations/chicken_genome/karyotype/dominance/
fi

for TOT_QTL in 100 1000 10000; do

  for REP in {1..20}; do
    
    if [ ! -d simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL} ]; then
      mkdir simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}
    fi

    Rscript R/simulate_karyotype.R \
      ${TOT_QTL} \
      0.2 \
      0.1 \
      annotation/chicken_genome_table.txt \
      simulations/chicken_genome/chicken_genome_founders.Rds \
      simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}/populations_${REP}.Rds \
      simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}/results_${REP}.Rds \
      simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}/simparam_${REP}.Rds \
      simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}/populations_gs_${REP}.Rds \
      simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}/results_gs_${REP}.Rds 
    
  done
  
done




############################################################################

## Cattle genome

## Dominance

if [ ! -d simulations/cattle_genome/karyotype/ ]; then
  mkdir simulations/cattle_genome/karyotype/
fi


if [ ! -d simulations/cattle_genome/karyotype/dominance/ ]; then
  mkdir simulations/cattle_genome/karyotype/dominance/
fi

for TOT_QTL in 100 1000 10000; do

  for REP in {1..20}; do
    
    if [ ! -d simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL} ]; then
      mkdir simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}
    fi

    Rscript R/simulate_karyotype.R \
      ${TOT_QTL} \
      0.2 \
      0.1 \
      annotation/cattle_genome_table.txt \
      simulations/cattle_genome/cattle_genome_founders.Rds \
      simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}/populations_${REP}.Rds \
      simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}/results_${REP}.Rds \
      simulations/cattle_genome/karyotype/dominance/totqtl${TOT_QTL}/simparam_${REP}.Rds \
      simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}/populations_gs_${REP}.Rds \
      simulations/chicken_genome/karyotype/dominance/totqtl${TOT_QTL}/results_gs_${REP}.Rds 
    
  done
  
done
