
## Create table of chromosome sizes in the chicken genome

library(dplyr)
library(readr)
library(tibble)

## Physical chromosome lengths

assembly_info <- read_tsv("annotation/GCF_000002315.6_GRCg6a_assembly_report.txt",
                          comment = "#",
                          col_names = FALSE)

chr_length <- tibble(chr = assembly_info$X10[1:32],
                     length = assembly_info$X9[1:32])



## Genetic chromosome length

genetic_map <- read_tsv("annotation/Supplemental_Table_S1_Groenen.txt",
                        col_types = "ccnnnnc")
colnames(genetic_map) <- c("chr", "marker", "female_cM", "male_cM",
                           "average_cM", "bp", "X")

genetic_map$chr <- paste("chr", genetic_map$chr, sep = "")

genetic_chr_length <- summarise(group_by(genetic_map, chr),
                                genetic_length = max(average_cM))


## Combine

chr_length_combined <- full_join(chr_length, genetic_chr_length)

chr_length_combined <- filter(chr_length_combined, !is.na(length))


## Set genetic length of microchromosomes missing from genetic map to the
## average of a few other microchromosomes

generic_micro_rate <- mean(filter(chr_length_combined, chr %in% paste("chr", 17:28, sep = ""))$genetic_length) /
  mean(filter(chr_length_combined, chr %in% paste("chr", 17:28, sep = ""))$length)

chr_length_combined$genetic_length[is.na(chr_length_combined$genetic_length)] <- 
  chr_length_combined$length[is.na(chr_length_combined$genetic_length)] *  
  generic_micro_rate


total_number_of_qtl <- 10 * 1000

chr_length_combined$physical_length_fraction <- chr_length_combined$length /
  sum(chr_length_combined$length)

chr_length_combined$n_qtl <- round(total_number_of_qtl * chr_length_combined$physical_length_fraction)


write.table(chr_length_combined,
            file ="annotation/chicken_genome_table.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)
