
library(dplyr)
library(ggplot2)

combined_phenotypic <- readRDS("outputs/summary_stats_phenotypic.Rds")
combined_genomic <- readRDS("outputs/summary_stats_genomic.Rds")




plot_stats_lines <- ggplot() +
  geom_line(aes(x = gen, y = average_gain, colour = factor(case)),
            data = combined_phenotypic) +
  geom_ribbon(aes(x = gen, y = average_gain, ymin = lower_gain, ymax = upper_gain, fill = factor(case)),
              data = combined_phenotypic, alpha = 1/4) +
  facet_wrap( ~ total_n_qtl)



combined_phenotypic_subset <- filter(combined_phenotypic,
                                     gen %in% c(5, 10, 15, 20),
                                     !(case %in% c("chicken uniform recombination",
                                                   "cattle uniform recombination")))


colours <- c("0.5 M chromosomes" = "grey0",
             "1 M chromosomes" = "grey30",
             "2 M chromosomes" = "grey60",
             "cattle real map" = "blue",
             "chicken real map" = "red")


plot_stats_bar <- ggplot() +
  geom_pointrange(aes(y = average_gain, x = gen,
                      ymin = lower_gain, ymax = upper_gain, colour = factor(case)),
                  data = combined_phenotypic_subset,
                  position = position_dodge(width = 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = colours) +
  facet_wrap( ~ total_n_qtl) +
  xlim(0, 22) +
  xlab("Generation") +
  ylab("Genetic mean")




combined_genomic_subset <- filter(combined_genomic,
                                  !(case %in% c("chicken uniform recombination",
                                                "cattle uniform recombination")))


plot_stats_gs_decline <- ggplot() +
  geom_line(aes(y = average_accuracy, x = gen,
                      colour = factor(case)),
                  data = combined_genomic_subset,
            size = 1) +
  scale_x_continuous(breaks = 1:10) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = colours) +
  facet_wrap( ~ total_n_qtl) +
  xlab("Generation") +
  ylab("Accuracy") +
  ylim(0, 1)

