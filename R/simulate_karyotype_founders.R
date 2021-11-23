
## Create founders for a simulation using real karyotypes


library(AlphaSimR)
library(readr)
library(tibble)
library(purrr)

source("R/simulation_functions.R")


args <- commandArgs(trailingOnly = TRUE)

genome_table_file <- args[1]
out_file <- args[2]

genome_table <- read_tsv(genome_table_file)


founder_chr <- pmap(list(bp = genome_table$length,
                         genLen = genome_table$genetic_length / 100,
                         Ne = 1000,
                         nInd = 1000,
                         nChr = 1),
                    runMacs2)

founders <- Reduce(cChr, founder_chr)


saveRDS(founders,
        file = out_file)

