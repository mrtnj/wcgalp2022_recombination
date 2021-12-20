
library(dplyr)
library(ggplot2)
library(patchwork)

combined_phenotypic <- readRDS("outputs/summary_stats_phenotypic.Rds")
combined_phenotypic$total_n_qtl <- paste(combined_phenotypic$total_n_qtl, "loci")
combined_genomic <- readRDS("outputs/summary_stats_genomic.Rds")
combined_genomic$total_n_qtl <- paste(combined_genomic$total_n_qtl, "loci")

linkage_equilibrium <- readRDS("outputs/summary_stats_linkage_equilibrium_start1.Rds")
linkage_equilibrium$total_n_qtl <- paste(linkage_equilibrium$total_n_qtl, "loci")

plot_stats_lines <- ggplot() +
  geom_line(aes(x = gen, y = average_mean_g, colour = factor(case)),
            data = combined_phenotypic) +
  geom_ribbon(aes(x = gen, y = average_mean_g, ymin = lower_mean_g, ymax = upper_mean_g, fill = factor(case)),
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



plot_stats_line <- ggplot() +
  geom_line(aes(y = average_mean_g, x = gen,
                colour = factor(case)),
            data = combined_phenotypic) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank()) +
  ##scale_colour_manual(values = colours) +
  facet_wrap( ~ total_n_qtl) +
  xlim(0, 22) +
  xlab("Generation") +
  ylab("Genetic mean")

plot_average_bar <- ggplot() +
  geom_pointrange(aes(y = average_mean_g, x = gen,
                      ymin = lower_mean_g,
                      ymax = upper_mean_g,
                      colour = factor(case)),
                  data = combined_phenotypic_subset,
                  position = position_dodge(width = 3)) +
  geom_line(aes(y = average_mean_g,
                 x = gen),
            linetype = 2,
            data = linkage_equilibrium) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = colours) +
  facet_wrap( ~ total_n_qtl) +
  xlim(0, 22) +
  xlab("Generation") +
  ylab("Genetic mean")


pdf("figures/phenotypic_selection_gain.pdf",
    height = 5,
    width = 10)
print(plot_average_bar)
dev.off()




plot_stats_var_a <- ggplot() +
  geom_pointrange(aes(y = average_var_a, x = gen,
                      ymin = lower_var_a,
                      ymax = upper_var_a,
                      colour = factor(case)),
                  data = combined_phenotypic_subset,
                  position = position_dodge(width = 3)) +
  geom_line(aes(y = average_var_a,
                x = gen),
            linetype = 2,
            data = linkage_equilibrium) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = colours) +
  facet_wrap( ~ total_n_qtl) +
  xlim(0, 22) +
  xlab("Generation") +
  ylab("Additive genetic variance")



pdf("figures/phenotypic_selection_additive_genetic_variance.pdf",
    height = 3,
    width = 10)
print(plot_stats_var_a)
dev.off()


## Combined genetic mean and additive variance figure

plot_combined1 <- plot_average_bar / plot_stats_var_a +
  plot_layout(guides = "collect") & theme(legend.position = "top")

pdf("figures/fig1_mean_var_a.pdf",
    width = 8.3,
    height = 8.3)
print(plot_combined1)
dev.off()





combined_genomic_subset <- filter(combined_genomic,
                                  !(case %in% c("chicken uniform recombination",
                                                "cattle uniform recombination")))


plot_stats_gs_decline <- ggplot() +
  geom_line(aes(y = average_accuracy, x = gen + 9,
                      colour = factor(case)),
                  data = combined_genomic_subset,
            size = 0.75) +
  scale_x_continuous(breaks = (1:10) + 9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = colours) +
  facet_wrap( ~ total_n_qtl) +
  xlab("Generation") +
  ylab("Accuracy") +
  ylim(0, 1)

pdf("figures/genomic_selection_accuracy.pdf",
    height = 3,
    width = 10)
print(plot_stats_gs_decline)
dev.off()


pdf("figures/fig2_genomic_selection_accuracy.pdf",
    width = 8.3,
    height = 8.3/2)
print(plot_stats_gs_decline & theme(legend.position = "top"))
dev.off()


plot_mean_g_gs <- ggplot() +
  geom_line(aes(y = average_mean_g, x = gen,
                colour = factor(case)),
            data = combined_genomic_subset,
            size = 0.75) +
  scale_x_continuous(breaks = 1:10) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = colours) +
  facet_wrap( ~ total_n_qtl) +
  xlab("Generation") +
  ylab("Genetic mean")




## Comparison between real map and karyotype only


colours_comparison <- c("blue", "lightblue", "red", "pink")

plot_real_map_accuracy <- ggplot() +
  geom_line(aes(y = average_accuracy, x = gen,
                colour = factor(case)),
            data = filter(combined_genomic,
                          case %in% c("chicken real map",
                                      "cattle real map",
                                      "chicken uniform recombination",
                                      "cattle uniform recombination")),
            size = 1) +
  scale_x_continuous(breaks = 1:10) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = colours_comparison) +
  facet_wrap( ~ total_n_qtl) +
  xlab("Generation") +
  ylab("Accuracy") +
  ylim(0, 1)



plot_real_map_mean_g <- ggplot() +
  geom_line(aes(y = average_mean_g, x = gen,
                colour = factor(case)),
            data = filter(combined_genomic,
                          case %in% c("chicken real map",
                                      "cattle real map",
                                      "chicken uniform recombination",
                                      "cattle uniform recombination")),
            size = 1) +
  scale_x_continuous(breaks = 1:10) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = colours_comparison) +
  facet_wrap( ~ total_n_qtl) +
  xlab("Generation") +
  ylab("Genetic mean")



pdf("figures/real_map_accuracy_comparison.pdf",
    height = 5,
    width = 10)
print(plot_real_map_accuracy)
dev.off()





## Values at the end

genetic_mean20 <- summarise(group_by(filter(combined_phenotypic, gen == 20 & total_n_qtl == "10000 loci"), case),
                            mean20 = mean(average_mean_g))

baseline_mean20 <- genetic_mean20$mean20[genetic_mean20$case == "0.5 M chromosomes"]

genetic_mean20$diff <- genetic_mean20$mean20 - baseline_mean20

genetic_mean20$percent_diff <- genetic_mean20$diff/baseline_mean20 * 100



accuracy2 <- summarise(group_by(filter(combined_genomic, gen == 2 & total_n_qtl == "10000 loci"), case),
                       acc2 = mean(average_accuracy))

baseline_accuracy <- accuracy2$acc2[accuracy2$case == "0.5 M chromosomes"]

accuracy2$diff <- accuracy2$acc2 - baseline_accuracy

accuracy2$percent_diff <- accuracy2$diff/baseline_accuracy * 100



