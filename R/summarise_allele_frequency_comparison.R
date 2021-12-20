
## Read the results from allele frequency comparison and make a table

library(tibble)
library(purrr)

## Chr length

n_reps <- 50
gen_length <- c(0.5, 1, 2)
nqtl <- 10000
replicate <- 1:n_reps

file_prefix_chr <- paste("simulations/chromosome_length/dominance/allele_frequency_comparison",
                         "/chrlen",
                         gen_length,
                         "_totqtl",
                         nqtl,
                         "/",
                         sep = "")


filename_chr <- unlist(map(file_prefix_chr,
                           function (x) paste(x,
                                              "correlations_asr",
                                              replicate,
                                              ".txt",
                                              sep = "")))

files_chr <- tibble(chr_length = rep(gen_length, each = n_reps),
                    n_qtl = rep(nqtl, each = n_reps * length(gen_length)),
                    replicate = rep(replicate, length(gen_length)),
                    filename = filename_chr,
                    case = as.character(1:length(filename_chr)))


## Real map

genome <- c("chicken", "cattle")

file_prefix_real_map <- paste("simulations/",
                              genome,
                              "_genome/real_map/dominance/allele_frequency_comparison/totqtl",
                              nqtl,
                              "/",
                              sep = "")

filename_real_map <- unlist(map(file_prefix_real_map,
                                function (x) paste(x,
                                                   "correlations_asr",
                                                   replicate,
                                                   ".txt",
                                                   sep = "")))


files_real_map <- tibble(species = rep(genome, each = n_reps),
                         n_qtl = rep(nqtl, each = n_reps * length(genome)),
                         replicate = rep(replicate, length(genome)),
                         filename = filename_real_map,
                         case = as.character(1:length(filename_real_map)))




## Linkage equilibrium

genome_le <- "cattle"

file_prefix_le <- paste("simulations/",
                        genome_le,
                        "_genome/cattle_genome/dominance/allele_frequency_comparison/totqtl",
                        nqtl,
                        "/",
                        sep = "")

filename_le <- unlist(map(file_prefix_real_map,
                          function (x) paste(x,
                                             "correlations_le",
                                             replicate,
                                             ".txt",
                                             sep = "")))


files_le <- tibble(species = rep(genome_le, each = n_reps),
                   n_qtl = rep(nqtl, each = n_reps * length(genome_le)),
                   replicate = rep(replicate, length(genome_le)),
                   filename = filename_le,
                   case = as.character(1:length(filename_le)))



## Read the results

read_correlations <- function(filename) {
  correlations <- scan(filename)
  
  tibble(correlation_gen2 = correlations[1],
         correlations_gen20 = correlations[2])
  
}

results_chr <- read_results(files_chr)

results_real_map <- read_results(files_real_map)

results_le <- read_results(files_le)