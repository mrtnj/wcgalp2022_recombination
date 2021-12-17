#!/bin/bash

## Run simulations with approximation of real karyotypes

set -eu



## Chicken 


## Dominance

if [ ! -d simulations/chicken_genome/real_map ]; then
  mkdir simulations/chicken_genome/real_map
fi

if [ ! -d simulations/chicken_genome/real_map/dominance/ ]; then
  mkdir simulations/chicken_genome/real_map/dominance/
fi

for TOT_QTL in 100 1000 10000; do

  for REP in {1..50}; do
    
    if [ ! -d simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL} ]; then
      mkdir simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}
    fi

    Rscript R/simulate_real_map.R \
      ${TOT_QTL} \
      0.2 \
      0.1 \
      annotation/chicken_genome_table.txt \
      simulations/chicken_genome/chicken_genome_founders.Rds \
      annotation/elferink2010_GRCg6a.txt \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/genotypes_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/phenotypes_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/results_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/simparam_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/genotypes_gs_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/phenotypes_gs_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/results_gs_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/generation1_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/generation10_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/generation20_${REP}.Rds \
      simulations/chicken_genome/real_map/dominance/totqtl${TOT_QTL}/model_${REP}.Rds 
    
  done
  
done




############################################################################

## Cattle genome

## Dominance

if [ ! -d simulations/cattle_genome/real_map ]; then
  mkdir simulations/cattle_genome/real_map
fi

if [ ! -d simulations/cattle_genome/real_map/dominance/ ]; then
  mkdir simulations/cattle_genome/real_map/dominance/
fi

for TOT_QTL in 100 1000 10000; do

  for REP in {1..50}; do
    
    if [ ! -d simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL} ]; then
      mkdir simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}
    fi

    Rscript R/simulate_real_map.R \
      ${TOT_QTL} \
      0.2 \
      0.1 \
      annotation/cattle_genome_table.txt \
      simulations/cattle_genome/cattle_genome_founders.Rds \
      annotation/ma2015_ARS-UCD1.2.txt \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/genotypes_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/phenotypes_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/results_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/simparam_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/genotypes_gs_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/phenotypes_gs_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/results_gs_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/generation1_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/generation10_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/generation20_${REP}.Rds \
      simulations/cattle_genome/real_map/dominance/totqtl${TOT_QTL}/model_${REP}.Rds
      
    
  done
  
done
