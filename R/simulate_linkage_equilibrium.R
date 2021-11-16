
library(AlphaSimR)
library(purrr)
library(tibble)

source("R/simulation_functions.R")
source("R/simulation_functions_linkage_equlibrium.R")


args <- commandArgs(trailingOnly = TRUE)

simparam_file <- args[1]
generations_file <- args[2]
starting_generation <- as.numeric(args[3])
selection_type <- args[4]
out_file <- args[5]


## Read starting point from ASR simulation

simparam <- readRDS(simparam_file)
generations <- readRDS(generations_file)

pop <- generations[starting_generation]


## Extract parameters

founder_geno <- pullQtlGeno(pop,
                            simParam = simparam)


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


saveRDS(generations_le,
        file = out_file)