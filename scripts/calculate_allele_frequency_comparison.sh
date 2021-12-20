#!/bin/bash

## Run allele frequency comparison from simulation outpu8t

set -eu



if [ ! -d simulations/chromosome_length/dominance/allele_frequency_comparison ]; then
  mkdir simulations/chromosome_length/dominance/allele_frequency_comparison
fi



## Chromosome length

for NQTL in 10000; do

  for GEN_LENGTH in 0.5 1 2; do

    for REP in {1..50}; do
    
      if [ ! -d simulations/chromosome_length/dominance/allele_frequency_comparison/chrlen${GEN_LENGTH}_totqtl${NQTL} ]; then
        mkdir simulations/chromosome_length/dominance/allele_frequency_comparison/chrlen${GEN_LENGTH}_totqtl${NQTL}
      fi
  
      Rscript R/calculate_allele_frequency_comparison.R \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/simparam_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/genotypes_${REP}.Rds \
        simulations/chromosome_length/dominance/chrlen${GEN_LENGTH}_totqtl${NQTL}/phenotypes_${REP}.Rds \
        NULL \
        simulations/chromosome_length/dominance/allele_frequency_comparison/chrlen${GEN_LENGTH}_totqtl${NQTL}/correlations_asr_${REP}.txt \
        NULL

    done
  
  done
  
done


## Cattle genome

if [ ! -d simulations/cattle_genome/real_map/dominance/allele_frequency_comparison/ ]; then
  mkdir simulations/cattle_genome/real_map/dominance/allele_frequency_comparison/
fi


for NQTL in 10000; do

    for REP in {1..50}; do
    
      if [ ! -d simulations/cattle_genome/real_map/dominance/allele_frequency_comparison/totqtl${NQTL} ]; then
        mkdir simulations/cattle_genome/real_map/dominance/allele_frequency_comparison/totqtl${NQTL} 
      fi
  
      Rscript R/calculate_allele_frequency_comparison.R \
        simulations/cattle_genome/real_map/dominance/totqtl${NQTL}/simparam_${REP}.Rds \
        simulations/cattle_genome/real_map/dominance/totqtl${NQTL}/genotypes_${REP}.Rds \
        simulations/cattle_genome/real_map/dominance/totqtl${NQTL}/phenotypes_${REP}.Rds \
        simulations/linkage_equilibrium/dominance/totqtl${NQTL}_start1/generations_${REP}.Rds \
        simulations/cattle_genome/real_map/dominance/allele_frequency_comparison/totqtl${NQTL}/correlations_asr_${REP}.txt \
        simulations/cattle_genome/real_map/dominance/allele_frequency_comparison/totqtl${NQTL}/correlations_le_${REP}.txt
  
  done
  
done


## Chicken genome

if [ ! -d simulations/chicken_genome/real_map/dominance/allele_frequency_comparison/ ]; then
  mkdir simulations/chicken_genome/real_map/dominance/allele_frequency_comparison/
fi


for NQTL in 10000; do

    for REP in {1..50}; do
    
      if [ ! -d simulations/chicken_genome/real_map/dominance/allele_frequency_comparison/totqtl${NQTL} ]; then
        mkdir simulations/chicken_genome/real_map/dominance/allele_frequency_comparison/totqtl${NQTL} 
      fi
  
      Rscript R/calculate_allele_frequency_comparison.R \
        simulations/chicken_genome/real_map/dominance/totqtl${NQTL}/simparam_${REP}.Rds \
        simulations/chicken_genome/real_map/dominance/totqtl${NQTL}/genotypes_${REP}.Rds \
        simulations/chicken_genome/real_map/dominance/totqtl${NQTL}/phenotypes_${REP}.Rds \
        NULL \
        simulations/chicken_genome/real_map/dominance/allele_frequency_comparison/totqtl${NQTL}/correlations_asr_${REP}.txt \
        NULL
  
  done
  
done
