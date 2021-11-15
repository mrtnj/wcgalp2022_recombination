
## Run simulation with approximation of actual karyotype


library(AlphaSimR)
library(readr)
library(tibble)
library(purrr)

source("R/simulation_functions.R")


args <- commandArgs(trailingOnly = TRUE)


total_qtl_number <- as.numeric(args[1])
mean_dd <- as.numeric(args[2])
var_dd <- as.numeric(args[3])
out_file_population <- args[4]
out_file_results <- args[5]


founders <- readRDS("simulations/chicken_genome/chicken_genome_founders.Rds")

genome_table <- read_tsv("annotation/chicken_genome_table.txt")

simulation <- make_simulation(founders,
                              n_qtl_per_chr = genome_table$n_qtl,
                              mean_dd = mean_dd,
                              var_dd = var_dd)

pop <- simulation$pop
simparam <- simulation$simparam

generations <- run_breeding(pop,
                            20,
                            simparam)

results <- get_stats(generations)


saveRDS(generations,
        file = out_file_population)

saveRDS(results,
        file = out_file_results)













