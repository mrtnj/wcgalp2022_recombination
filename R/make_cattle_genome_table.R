
## Create table of chromosome sizes in the chicken genome

library(dplyr)
library(readr)
library(tibble)

## Physical chromosome lengths

assembly_info <- read_tsv("annotation/GCF_002263795.1_ARS-UCD1.2_assembly_report.txt",
                          comment = "#",
                          col_names = FALSE)

chr_length <- tibble(chr = assembly_info$X10[1:32],
                     length = assembly_info$X9[1:32])



## Genetic chromosome length

genetic_map <- read_tsv("annotation/ma2015_ARS-UCD1.2.txt")


genetic_chr_length <- summarise(group_by(genetic_map, chr),
                                genetic_length = max(position_cM))


## Combine

chr_length_combined <- inner_join(chr_length, genetic_chr_length)




total_number_of_qtl <- 10 * 1000

chr_length_combined$physical_length_fraction <- chr_length_combined$length /
  sum(chr_length_combined$length)

chr_length_combined$n_qtl <- round(total_number_of_qtl * chr_length_combined$physical_length_fraction)


write.table(chr_length_combined,
            file ="annotation/cattle_genome_table.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)
