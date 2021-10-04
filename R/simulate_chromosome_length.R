## Simulate chromosome length variation

library(AlphaSimR)
library(tibble)

source("R/simulation_functions.R")


args <- commandArgs(trailingOnly = TRUE)


gen_length <- as.numeric(args[1])
n_qtl_per_chr <- as.numeric(args[2])
mean_dd <- as.numeric(args[3])
var_dd <- as.numeric(args[4])
out_file_population <- args[5]
out_file_results <- args[6]


founders <- runMacs2(nInd = 1000,
                     nChr = 10,
                     genLen = gen_length)

simulation <- make_simulation(founders,
                              n_qtl_per_chr = n_qtl_per_chr,
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

