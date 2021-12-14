
## Read results files and summarise them

library(dplyr)
library(purrr)
library(tibble)



## Chromosome length
n_reps <- 50

gen_length <- rep(rep(c(0.5, 1, 2), each = 3), 2) 
nqtl <- rep(rep(c(100, 1000, 10000), 3), 2)
replicate <- 1:n_reps

file_prefix_chr <- paste("simulations/chromosome_length/dominance",
                         "/chrlen",
                         gen_length,
                         "_totqtl",
                         nqtl,
                         "/",
                         sep = "")

filename_chr <- unlist(map(file_prefix_chr,
                           function (x) paste(x,
                                              "results_",
                                              replicate,
                                              ".Rds",
                                              sep = "")))

filename_gs_chr <- unlist(map(file_prefix_chr,
                              function (x) paste(x,
                                                 "results_gs_",
                                                 replicate,
                                                 ".Rds",
                                                 sep = "")))


files_chr <- tibble(chr_length = rep(gen_length, each = n_reps),
                    n_qtl = rep(nqtl, each = n_reps),
                    replicate = rep(replicate, length(gen_length)),
                    filename = filename_chr,
                    filename_gs = filename_gs_chr,
                    case = as.character(1:length(filename_chr)))




## Karyotype

genome <- rep(c("chicken", "cattle"), each = 3)
nqtl_karyotype <- rep(c(100, 1000, 10000), 2)


file_prefix_karyotype <- paste("simulations/",
                               genome,
                               "_genome/karyotype/dominance/totqtl",
                               nqtl_karyotype,
                               "/",
                               sep = "")

filename_karyotype <- unlist(map(file_prefix_karyotype,
                                 function (x) paste(x,
                                                    "results_",
                                                    replicate,
                                                    ".Rds",
                                                    sep = "")))

filename_gs_karyotype <- unlist(map(file_prefix_karyotype,
                                    function (x) paste(x,
                                                       "results_gs_",
                                                       replicate,
                                                       ".Rds",
                                                       sep = "")))


files_karyotype <- tibble(species = rep(genome, each = n_reps),
                          n_qtl = rep(nqtl_karyotype, each = n_reps),
                          replicate = rep(replicate, length(genome)),
                          filename = filename_karyotype,
                          filename_gs = filename_gs_karyotype,
                          case = as.character(1:length(filename_karyotype)))



## Real map

file_prefix_real_map <- paste("simulations/",
                              genome,
                              "_genome/real_map/dominance/totqtl",
                              nqtl_karyotype,
                              "/",
                              sep = "")

filename_real_map <- unlist(map(file_prefix_real_map,
                                function (x) paste(x,
                                                   "results_",
                                                   replicate,
                                                   ".Rds",
                                                   sep = "")))

filename_gs_real_map <- unlist(map(file_prefix_real_map,
                                   function (x) paste(x,
                                                      "results_gs_",
                                                      replicate,
                                                      ".Rds",
                                                      sep = "")))


files_real_map <- tibble(species = rep(genome, each = n_reps),
                n_qtl = rep(nqtl_karyotype, each = n_reps),
                replicate = rep(replicate, length(genome)),
                filename = filename_real_map,
                filename_gs = filename_gs_real_map,
                case = as.character(1:length(filename_real_map)))



read_results <- function(files) {

  results <- map_dfr(files$filename,
                     readRDS,
                     .id = "case")
  
  
  results_metadata <- inner_join(results, files)

  
  results_gs <- map_dfr(files$filename_gs,
                        readRDS,
                        .id = "case")
  
  
  results_gs_metadata <- inner_join(results_gs, files)
  
  list(results = results_metadata,
       results_gs = results_gs_metadata)
}  


results_chr <- read_results(files_chr)

results_karyotype <- read_results(files_karyotype)

results_real_map <- read_results(files_real_map)




## Summarisation

get_summary_stats_chr <- function(results) {
  summarise(group_by(results, gen, chr_length, n_qtl),
            average_gain = mean(mean_g),
            lower_gain = quantile(mean_g, 0.05),
            upper_gain = quantile(mean_g, 0.95),
            average_var_a = mean(var_a),
            lower_var_a = quantile(var_a, 0.05),
            upper_var_a = quantile(var_a, 0.95),
            average_accuracy = mean(accuracy, na.rm = TRUE),
            lower_accuracy = quantile(accuracy, 0.05, na.rm = TRUE),
            upper_accuracy = quantile(accuracy, 0.95, na.rm = TRUE))
}

get_summary_stats_karyotype <- function(results) {
  summarise(group_by(results, gen, species, n_qtl),
            average_gain = mean(mean_g),
            lower_gain = quantile(mean_g, 0.05),
            upper_gain = quantile(mean_g, 0.95),
            average_var_a = mean(var_a),
            lower_var_a = quantile(var_a, 0.05),
            upper_var_a = quantile(var_a, 0.95),
            average_accuracy = mean(accuracy, na.rm = TRUE),
            lower_accuracy = quantile(accuracy, 0.05, na.rm = TRUE),
            upper_accuracy = quantile(accuracy, 0.95, na.rm = TRUE))
}

summary_stats_chr <- get_summary_stats_chr(results_chr$results)
summary_stats_chr_gs <- get_summary_stats_chr(results_chr$results_gs)

summary_stats_karyotype <- get_summary_stats_karyotype(results_karyotype$results)
summary_stats_karyotype_gs <- get_summary_stats_karyotype(results_karyotype$results_gs)

summary_stats_real_map <- get_summary_stats_karyotype(results_real_map$results)
summary_stats_real_map_gs <- get_summary_stats_karyotype(results_real_map$results_gs)



summary_stats_chr$case <- paste(summary_stats_chr$chr_length, "M chromosomes")
summary_stats_chr$total_n_qtl <- summary_stats_chr$n_qtl

summary_stats_chr_gs$case <- paste(summary_stats_chr_gs$chr_length, "M chromosomes")
summary_stats_chr_gs$total_n_qtl <- summary_stats_chr_gs$n_qtl

summary_stats_karyotype$case <- paste(summary_stats_karyotype$species, "uniform recombination")
summary_stats_karyotype$total_n_qtl <- summary_stats_karyotype$n_qtl

summary_stats_karyotype_gs$case <- paste(summary_stats_karyotype_gs$species, "uniform recombination")
summary_stats_karyotype_gs$total_n_qtl <- summary_stats_karyotype_gs$n_qtl

summary_stats_real_map$case <- paste(summary_stats_real_map$species, "real map")
summary_stats_real_map$total_n_qtl <- summary_stats_real_map$n_qtl

summary_stats_real_map_gs$case <- paste(summary_stats_real_map_gs$species, "real map")
summary_stats_real_map_gs$total_n_qtl <- summary_stats_real_map_gs$n_qtl


cols <- c("gen", "average_gain", "lower_gain", "upper_gain",
          "average_var_a", "lower_var_a", "upper_var_a",
          "average_accuracy", "lower_accuracy", "upper_accuracy",
          "case", "total_n_qtl")

combined_phenotypic <- rbind(summary_stats_chr[, cols],
                             summary_stats_karyotype[, cols],
                             summary_stats_real_map[, cols])

combined_gs <- rbind(summary_stats_chr_gs[, cols],
                     summary_stats_karyotype_gs[, cols],
                     summary_stats_real_map_gs[, cols])


saveRDS(combined_phenotypic,
        file = "outputs/summary_stats_phenotypic.Rds")

saveRDS(combined_gs,
        file = "outputs/summary_stats_genomic.Rds")
