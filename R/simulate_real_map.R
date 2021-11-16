
## Run simulation with approximation of actual karyotype and genetic map


library(AlphaSimR)
library(readr)
library(tibble)
library(purrr)

source("R/simulation_functions.R")


args <- commandArgs(trailingOnly = TRUE)


total_qtl_number <- as.numeric(args[1])
mean_dd <- as.numeric(args[2])
var_dd <- as.numeric(args[3])
genome_table_file <- args[4] ## e.g. "annotation/chicken_genome_table.txt"
founder_file <- args[5] ## e.g. "simulations/chicken_genome/chicken_genome_founders.Rds"
map_file <- args[6] ## e.g. "annotation/elferink2010_GRCg6a.txt"
out_file_population <- args[7]
out_file_results <- args[8]
out_file_simparam <- args[9]


genome_table <- read_tsv(genome_table_file)


## Prepare linkage map

real_map <- read_tsv("annotation/elferink2010_GRCg6a.txt")

real_map <- filter(real_map, chr %in% genome_table$chr)
real_map_chr <- split(real_map, real_map$chr)
real_map_chr <- real_map_chr[match(genome_table$chr, names(real_map_chr))]


founders <- readRDS(founder_file)


simulation <- make_simulation(founders,
                              n_qtl_per_chr = round(genome_table$physical_length_fraction * total_qtl_number),
                              mean_dd = mean_dd,
                              var_dd = var_dd)

pop <- simulation$pop
simparam <- simulation$simparam


## Create adjusted map and switch it out
new_map <- make_adjusted_map(simparam$genMap,
                             real_map_chr)

simparam$switchGenMap(genmap = new_map)


generations <- run_breeding(pop,
                            20,
                            simparam)

results <- get_stats(generations)


saveRDS(generations,
        file = out_file_population)

saveRDS(results,
        file = out_file_results)

save(simparam,
     file = out_file_simparam)












