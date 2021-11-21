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
out_file_geno <- args[5]
out_file_pheno <- args[6]
out_file_results <- args[7]
out_file_simparam <- args[8]
out_file_geno_gs <- args[9]
out_file_pheno_gs <- args[10]
out_file_results_gs <- args[11]


founders <- runMacs2(nInd = 1000,
                     Ne = 1000,
                     nChr = 10,
                     genLen = gen_length)

n_snps_per_chr <- 5000

simulation <- make_simulation(founders,
                              n_qtl_per_chr = n_qtl_per_chr,
                              n_snps_per_chr = n_snps_per_chr,
                              mean_dd = mean_dd,
                              var_dd = var_dd)

## Phenotypic selection

pop <- simulation$pop
simparam <- simulation$simparam

generations <- run_breeding(pop,
                            20,
                            simparam)

results <- get_stats(generations, simparam)



## Genomic selection

training <- Reduce(c, generations[7:10])

model <- RRBLUP(pop = training,
                simParam = simparam)


generations_gs <- run_genomic_selection(generations[[10]],
                                        model,
                                        10,
                                        simparam)

results_gs <- get_stats(generations_gs, simparam)



## Save genotypes and phenotypes

genotypes <- map(generations,
                 pullQtlGeno,
                 simParam = simparam)

genotypes_gs <- map(generations_gs,
                    pullQtlGeno,
                    simParam = simparam)


phenotypes <- map(generations,
                  get_pheno_stats,
                  simparam = simparam)

phenotypes_gs <- map(generations_gs,
                     get_pheno_stats,
                     simparam = simparam)



saveRDS(genotypes,
        file = out_file_geno)

saveRDS(phenotypes,
        file = out_file_pheno)

saveRDS(results,
        file = out_file_results)

save(simparam,
     file = out_file_simparam)

saveRDS(genotypes_gs,
        file = out_file_geno_gs)

saveRDS(phenotypes_gs,
        file = out_file_pheno_gs)

saveRDS(results_gs,
        file = out_file_results_gs)