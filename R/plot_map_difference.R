

library(AlphaSimR)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)

source("R/simulation_functions.R")

chicken_genome_table <- read_tsv("annotation/chicken_genome_table.txt")
cattle_genome_table <- read_tsv("annotation/cattle_genome_table.txt")


## Prepare linkage maps

real_map_chicken <- read_tsv("annotation/elferink2010_GRCg6a.txt")
real_map_cattle <- read_tsv("annotation/ma2015_ARS-UCD1.2.txt")

real_map_chicken <- filter(real_map_chicken, chr %in% chicken_genome_table$chr)
real_map_chicken_split <- split(real_map_chicken, real_map_chicken$chr)
real_map_chicken_split <- real_map_chicken_split[match(chicken_genome_table$chr, names(real_map_chicken_split))]

real_map_cattle <- filter(real_map_cattle, chr %in% cattle_genome_table$chr)
real_map_cattle_split <- split(real_map_cattle, real_map_cattle$chr)
real_map_cattle_split <- real_map_cattle_split[match(cattle_genome_table$chr, names(real_map_cattle_split))]



founders_chicken <- readRDS("simulations/chicken_genome/chicken_genome_founders.Rds")
founders_cattle <- readRDS("simulations/cattle_genome/cattle_genome_founders.Rds")




simulation_chicken <- make_simulation(founders_chicken,
                                      n_qtl_per_chr = 100,
                                      n_snps_per_chr = 100,
                                      mean_dd = 0,
                                      var_dd = 0)

simparam_chicken <- simulation_chicken$simparam


simulation_cattle <- make_simulation(founders_cattle,
                                     n_qtl_per_chr = 100,
                                     n_snps_per_chr = 100,
                                     mean_dd = 0,
                                     var_dd = 0)

simparam_cattle <- simulation_cattle$simparam


## Create adjusted map and switch it out
new_map_chicken <- make_adjusted_map(simparam_chicken$genMap,
                                     real_map_chicken_split)

new_map_cattle <- make_adjusted_map(simparam_cattle$genMap,
                                     real_map_cattle_split)


names(new_map_chicken) <- names(real_map_chicken_split)
names(new_map_cattle) <- names(real_map_cattle_split)

flatten_map <- function(map) {
  
  map_dfr(map,
          function(chr)
            data.frame(pos = chr),
          .id = "chr_name")
  
}


chicken_df <- flatten_map(new_map_chicken)
chicken_old_df <- flatten_map(simparam_chicken$genMap)
chicken_df$old_pos <- chicken_old_df$pos


qplot(x = old_pos, y = pos, data = chicken_df) + facet_wrap(~ chr_name, scale = "free")

cattle_df <- flatten_map(new_map_cattle)
cattle_old_df <- flatten_map(simparam_cattle$genMap)
cattle_df$old_pos <- cattle_old_df$pos


qplot(x = old_pos, y = pos, data = cattle_df) + facet_wrap(~ chr_name, scale = "free")



qplot(x = old_pos, y = pos, data = cattle_df, geom = "line", group = chr_name)

qplot(x = old_pos, y = pos, data = chicken_df, geom = "line", group = chr_name)



### Focus in on cattle chr 5 and chicken chr8

case_studies <- rbind(transform(filter(chicken_df, chr_name == "chr8"),
                                case = "chicken chr8"),
                      transform(filter(cattle_df, chr_name == "chr5"),
                                case = "cattle chr5"))


plot_case_studies <- ggplot() +
  geom_point(aes(x = old_pos, y = pos, colour = case),
             data = case_studies,
             size = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.3)) +
  scale_colour_manual(values = c("blue", "red")) +
  ggtitle("Comparison of uniform and real recombination map") +
  xlab("Position on uniform map (M)") +
  ylab("Position on real map (M)")



pdf("figures/map_comparison.pdf")
print(plot_case_studies)
dev.off()
