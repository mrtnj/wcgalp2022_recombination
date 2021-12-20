
## Compare allele frequency changes over the first generation of genomic selection

library(tibble)

source("R/simulation_functions.R")

args <- commandArgs(trailingOnly = TRUE)

simparam_file <- args[1] ## "simulations/cattle_genome/real_map/dominance/totqtl10000/simparam_1.Rds"
genotype_file <- args[2] ## "simulations/cattle_genome/real_map/dominance/totqtl10000/genotypes_1.Rds"
phenotype_file <- args[3] ## "simulations/cattle_genome/real_map/dominance/totqtl10000/phenotypes_gs_1.Rds"

le_generations_file <- args[4] ## "simulations/linkage_equilibrium/dominance/totqtl10000_start1/generations_1.Rds"

if (le_generations_file == "NULL") {
  le_generations_file <- NULL 
}

outfile_correlation_asr <- args[5]
outfile_correlation_le <- args[6]



## Load data from real map simulation 

genotypes <- readRDS(genotype_file)


## Load simparam to find QTL

simparam <- readRDS(simparam_file)


## Load results to get varP

pheno <- readRDS(phenotype_file)

var_p <- var(pheno[[1]]$pheno)



## Extract allele frequency change over first generation and over whole simulation

geno1 <- genotypes[[1]]
geno2 <- genotypes[[2]]
geno20 <- genotypes[[20]]

f_gen1 <- get_qtl_low_allele_freq(geno1, simparam)
f_gen2 <- get_qtl_low_allele_freq(geno2, simparam)
f_gen20 <- get_qtl_low_allele_freq(geno20, simparam)


## Calculate expected frequency change from theory

freq_change_asr2 <- compare_predicted_frequency_change(f_gen1,
                                                       f_gen2,
                                                       prop_male = 50/500,
                                                       prop_female = 100/500,
                                                       var_p = var_p,
                                                       simparam = simparam)


freq_change_asr20 <- compare_predicted_frequency_change(f_gen1,
                                                        f_gen20,
                                                        prop_male = 50/500,
                                                        prop_female = 100/500,
                                                        var_p = var_p,
                                                        simparam = simparam)



## Save data

correlation <- c(cor(freq_change_asr2$predicted_delta_f, freq_change_asr2$delta_f),
                 cor(freq_change_asr20$predicted_delta_f, freq_change_asr20$delta_f))

write(correlation,
      file = outfile_correlation_asr)





## If provided, also test the matching linkage equilibrium simulation

if (! is.null(le_generations_file)) {

  generations_le <- readRDS(le_generations_file)

  var_p_le <- var(generations_le[[1]]$phenotypes)
  
  f_gen1_le <- get_qtl_low_allele_freq_linkage_equilibrium(generations_le[[1]]$freq, simparam)
  f_gen2_le <- get_qtl_low_allele_freq_linkage_equilibrium(generations_le[[2]]$freq, simparam)
  f_gen20_le <- get_qtl_low_allele_freq_linkage_equilibrium(generations_le[[20]]$freq, simparam)
  
  freq_change_le2 <- compare_predicted_frequency_change(f_gen1_le,
                                                        f_gen2_le,
                                                        prop_male = 50/500,
                                                        prop_female = 100/500,
                                                        var_p = var_p_le,
                                                        simparam = simparam)
  
  freq_change_le20 <- compare_predicted_frequency_change(f_gen1_le,
                                                         f_gen20_le,
                                                         prop_male = 50/500,
                                                         prop_female = 100/500,
                                                         var_p = var_p_le,
                                                         simparam = simparam)
  
  correlation_le <- c(cor(freq_change_le2$predicted_delta_f, freq_change_le2$delta_f),
                      cor(freq_change_le20$predicted_delta_f, freq_change_le20$delta_f))
  
  write(correlation,
        file = outfile_correlation_le)
  
}


