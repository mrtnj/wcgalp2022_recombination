
## Summarise the results of linkage equilibrium simulations

library(dplyr)
library(purrr)
library(tibble)


nqtl <- c(100, 1000, 10000)
replicate <- 1:50


file_prefix_start1 <- paste("simulations/linkage_equilibrium/dominance/",
                            "totqtl",
                            nqtl,
                            "_start1/",
                            sep = "")

filename_start1 <- unlist(map(file_prefix_start1,
                              function (x) paste(x,
                                                 "results_",
                                                 replicate,
                                                 ".Rds",
                                                 sep = "")))

files_start1 <- tibble(total_n_qtl = rep(nqtl, each = 50),
                       replicate = rep(replicate, 3),
                       filename = filename_start1,
                       case = as.character(1:length(filename_start1)))




read_results <- function(files) {
  
  results <- map_dfr(files$filename,
                     readRDS,
                     .id = "case")
  
  
  inner_join(results, files)
  
}


results_start1 <- read_results(files_start1)



get_summary_stats <- function(results) {
  summarise(group_by(results, gen, n_qtl),
            average_mean_g = mean(mean_g),
            sd_mean_g = sd(mean_g),
            lower_mean_g = quantile(mean_g, 0.05),
            upper_mean_g = quantile(mean_g, 0.95),
            average_var_a = mean(var_a),
            sd_var_a = sd(var_a),
            lower_var_a = quantile(var_a, 0.05),
            upper_var_a = quantile(var_a, 0.95))
}

stats_start1 <- get_summary_stats(results_start1)



saveRDS(stats_start1,
        file = "outputs/summary_stats_linkage_equilibrium_start1.Rds")
