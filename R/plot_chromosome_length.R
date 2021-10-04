
library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(tibble)


gen_length <- c(0.5, 1, 1.5, 2)
replicate <- 1:20

file_prefix <- paste("simulations/chromosome_length/gen_length_",
                     gen_length,
                     "/results_", sep = "")

filename <- unlist(map(file_prefix,
                       function (x) paste(x,
                                          replicate,
                                          ".Rds",
                                          sep = "")))


files <- tibble(chr_length = rep(gen_length, each = 20),
                replicate = rep(replicate, length(gen_length)),
                filename = filename,
                case = as.character(1:length(filename)))


results <- map_dfr(files$filename,
                   readRDS,
                   .id = "case")


results_metadata <- inner_join(results, files)




summary_stats <- summarise(group_by(results_metadata, gen, chr_length),
          average = mean(mean_g),
          lower = quantile(mean_g, 0.05),
          upper = quantile(mean_g, 0.95))


plot_lines <- qplot(x = gen, y = mean_g, data = results_metadata, colour = chr_length, geom = "line", group = paste(replicate, chr_length))


plot_stats <- ggplot() +
  geom_pointrange(aes(x = gen, y = average, ymin = lower, ymax = upper, colour = factor(chr_length)),
                  data = summary_stats,
                  position = position_dodge(width = 1))

