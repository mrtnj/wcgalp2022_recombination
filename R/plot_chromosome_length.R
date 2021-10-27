
library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(tibble)


architecture <- rep(c("additive", "dominance"), each = 9)
gen_length <- rep(rep(c(0.5, 1, 2), each = 3), 2) 
nqtl <- rep(rep(c(10, 100, 1000), 3), 2)
replicate <- 1:20

file_prefix <- paste("simulations/chromosome_length/",
                     architecture,
                     "/chrlen",
                     chr_length,
                     "_nqtl",
                     nqtl,
                     "/",
                     sep = "")

filename <- unlist(map(file_prefix,
                       function (x) paste(x,
                                          "results_",
                                          replicate,
                                          ".Rds",
                                          sep = "")))


files <- tibble(chr_length = rep(gen_length, each = 20),
                n_qtl = rep(nqtl, each = 20),
                architecture = rep(architecture, each = 20),
                replicate = rep(replicate, length(gen_length)),
                filename = filename,
                case = as.character(1:length(filename)))


results <- map_dfr(files$filename,
                   readRDS,
                   .id = "case")


results_metadata <- inner_join(results, files)




summary_stats <- summarise(group_by(results_metadata, gen, chr_length, n_qtl, architecture),
          average = mean(mean_g),
          lower = quantile(mean_g, 0.05),
          upper = quantile(mean_g, 0.95))


plot_lines <- qplot(x = gen, y = mean_g,
                    data = results_metadata,
                    colour = chr_length,
                    geom = "line", group = paste(replicate, chr_length)) +
  facet_wrap(architecture ~ n_qtl)


plot_stats <- ggplot() +
  geom_pointrange(aes(x = gen, y = average, ymin = lower, ymax = upper, colour = factor(chr_length)),
                  data = summary_stats,
                  position = position_dodge(width = 1)) +
  facet_wrap(architecture ~ n_qtl)


plot_stats_lines <- ggplot() +
  geom_line(aes(x = gen, y = average, colour = factor(chr_length)),
            data = summary_stats) +
  geom_ribbon(aes(x = gen, y = average, ymin = lower, ymax = upper, fill = factor(chr_length)),
              data = summary_stats, alpha = 1/4) +
  facet_wrap(architecture ~ n_qtl)

