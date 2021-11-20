## Simulate chromosome length variation

library(AlphaSimR)
library(tibble)
library(purrr)

source("R/simulation_functions.R")


args <- commandArgs(trailingOnly = TRUE)


gen_length <- as.numeric(args[1])
n_qtl_per_chr <- as.numeric(args[2])
mean_dd <- as.numeric(args[3])
var_dd <- as.numeric(args[4])
out_file_population <- args[5]
out_file_results <- args[6]
out_file_simparam <- args[7]
out_file_population_gs <- args[8]
out_file_results_gs <- args[9]


founders <- runMacs2(nInd = 1000, 
                     nChr = 10,
                     genLen = gen_length)

simulation <- make_simulation(founders,
                              n_qtl_per_chr = n_qtl_per_chr,
                              n_snps_per_chr = 5000,
                              mean_dd = mean_dd,
                              var_dd = var_dd)

## Phenotypic selection

pop <- simulation$pop
simparam <- simulation$simparam

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