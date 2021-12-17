## Simulate chromosome length variation

library(AlphaSimR)
library(tibble)
library(purrr)

source("R/simulation_functions.R")


## Set up dumping of environment for debugging purposes

options(error = quote(dump.frames("dump_chr_length", TRUE, TRUE)))


args <- commandArgs(trailingOnly = TRUE)


gen_length <- as.numeric(args[1])
n_total_qtl <- as.numeric(args[2])
mean_dd <- as.numeric(args[3])
var_dd <- as.numeric(args[4])
out_file_geno <- args[5]
out_file_pheno <- args[6]
out_file_results <- args[7]
out_file_simparam <- args[8]
out_file_geno_gs <- args[9]
out_file_pheno_gs <- args[10]
out_file_results_gs <- args[11]
out_file_generation1 <- args[12]
out_file_generation10 <- args[13]
out_file_generation20 <- args[14]


founders <- runMacs2(nInd = 1000,
                     Ne = 1000,
                     nChr = 20,
                     genLen = gen_length)

n_snps_per_chr <- 5000

print("Make simulation")

simulation <- make_simulation(founders,
                              n_qtl_per_chr = n_total_qtl / 20,
                              n_snps_per_chr = n_snps_per_chr,
                              mean_dd = mean_dd,
                              var_dd = var_dd)

## Phenotypic selection

print("Phenotypic selection")

pop <- simulation$pop
simparam <- simulation$simparam

generations <- run_breeding(pop,
                            20,
                            simparam)

print("Get results")

results <- get_stats(generations, simparam)



## Genomic selection

print("Training genomic model")

training <- Reduce(c, generations[7:10])

model <- RRBLUP(pop = training,
                simParam = simparam)


print("Genomic selection")

generations_gs <- run_genomic_selection(generations[[10]],
                                        model,
                                        10,
                                        simparam)

print("Get results")

results_gs <- get_stats(generations_gs, simparam)



## Save genotypes and phenotypes

print("Save genotypes and phenotypes")

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

saveRDS(simparam,
     file = out_file_simparam)

saveRDS(genotypes_gs,
        file = out_file_geno_gs)

saveRDS(phenotypes_gs,
        file = out_file_pheno_gs)

saveRDS(results_gs,
        file = out_file_results_gs)

saveRDS(generations[[1]],
        file = out_file_generation1)

saveRDS(generation10[[10]],
        file = out_file_generation10)

saveRDS(generation20[[20]],
        file = out_file_generation20)