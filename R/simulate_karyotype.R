
## Run simulation with approximation of actual karyotype


library(AlphaSimR)
library(readr)
library(tibble)
library(purrr)

source("R/simulation_functions.R")


## Set up dumping of environment for debugging purposes

options(error = quote(dump.frames("dump_karyotype", TRUE, TRUE)))




args <- commandArgs(trailingOnly = TRUE)


total_qtl_number <- as.numeric(args[1])
mean_dd <- as.numeric(args[2])
var_dd <- as.numeric(args[3])
genome_table_file <- args[4] ## e.g. "annotation/chicken_genome_table.txt"
founder_file <- args[5] ## e.g. "simulations/chicken_genome/chicken_genome_founders.Rds"
out_file_geno <- args[6]
out_file_pheno <- args[7]
out_file_results <- args[8]
out_file_simparam <- args[9]
out_file_geno_gs <- args[10]
out_file_pheno_gs <- args[11]
out_file_results_gs <- args[12]


genome_table <- read_tsv(genome_table_file)

founders <- readRDS(founder_file)


print("Make simulation")

simulation <- make_simulation(founders,
                              n_qtl_per_chr = round(genome_table$physical_length_fraction * total_qtl_number),
                              n_snps_per_chr = round(genome_table$physical_length_fraction * 50000),
                              mean_dd = mean_dd,
                              var_dd = var_dd)

pop <- simulation$pop
simparam <- simulation$simparam


## Phenotypic selection

print("Phenotypic selection")

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