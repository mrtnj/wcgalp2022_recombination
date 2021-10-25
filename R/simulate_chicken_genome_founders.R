
## Create founders for chicken karyotype simulations


library(AlphaSimR)
library(readr)
library(tibble)
library(purrr)

source("R/simulation_functions.R")



genome_table <- read_tsv("annotation/chicken_genome_table.txt")


founder_chr <- pmap(list(bp = genome_table$length,
                         genLen = genome_table$genetic_length,
                         Ne = 1000,
                         nInd = 1000,
                         nChr = 1),
                    runMacs2)

founders <- Reduce(cChr, founder_chr)


saveRDS(founders,
        file = "simulations/chicken_genome/chicken_genome_founders.Rds")

