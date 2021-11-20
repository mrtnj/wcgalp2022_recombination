
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
genome_table_file <- args[4] ## e.g. "annotation/chicken_genome_table.txt"
founder_file <- args[5] ## e.g. "simulations/chicken_genome/chicken_genome_founders.Rds"
out_file_population <- args[6]
out_file_results <- args[7]
out_file_simparam <- args[8]
out_file_population_gs <- args[9]
out_file_results_gs <- args[10]


genome_table <- read_tsv(genome_table_file)

founders <- readRDS(founder_file)


simulation <- make_simulation(founders,
                              n_qtl_per_chr = round(genome_table$physical_length_fraction * total_qtl_number),
                              n_snps_per_chr = round(genome_table$physical_length_fraction * 50000),
                              mean_dd = mean_dd,
                              var_dd = var_dd)

pop <- simulation$pop
simparam <- simulation$simparam

## Phenotypic selection

generations <- run_breeding(pop,
                            40,
                            simparam)

results <- get_stats(generations)



## Genomic selection

training <- Reduce(c, generations[17:20])

model <- RRBLUP(pop = training,
                simParam = simparam)


generations_gs <- run_genomic_selection(generations[[20]],
                                        model,
                                        20,
                                        simparam)

results_gs <- get_stats(generations_gs)


saveRDS(generations,
        file = out_file_population)

saveRDS(results,
        file = out_file_results)

save(simparam,
     file = out_file_simparam)

saveRDS(generations_gs,
        file = out_file_population_gs)

saveRDS(results_gs,
        file = out_file_results_gs)

