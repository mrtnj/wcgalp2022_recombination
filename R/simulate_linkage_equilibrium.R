
library(AlphaSimR)
library(purrr)
library(tibble)

source("R/simulation_functions.R")
source("R/simulation_functions_linkage_equlibrium.R")



args <- commandArgs(trailingOnly = TRUE)

simparam_file <- args[1]
genotypes_file <- args[2]
starting_generation <- as.numeric(args[3])
selection_type <- args[4]
out_file_results <- args[5]
out_file_generations <- args[6]


## Read starting point from ASR simulation

load(simparam_file)
genotypes <- readRDS(genotypes_file)


## Extract parameters

founder_geno <- genotypes[[starting_generation]]


a <- simparam$traits[[1]]@addEff
d <- simparam$traits[[1]]@domEff

Ve <- simparam$varE

n_ind <- 1000
n_loc <- ncol(founder_geno)


generations_le <-
  run_selection_linkage_equilibrium(founder_geno,
                                    a = a,
                                    d = d,
                                    Ve = Ve,
                                    n_ind = n_ind,
                                    n_gen = 20,
                                    prop = 0.15,
                                    selection_type = selection_type)


saveRDS(generations_le$results,
        file = out_file_results)

saveRDS(generations_le$generations,
        file = out_file_generations)