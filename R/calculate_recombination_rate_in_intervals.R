
## Get recombiantion rate for between-marker intervals

library(readr)
library(purrr)

cattle_map <- read_tsv("annotation/ma2015_ARS-UCD1.2.txt")
chicken_map <- read_tsv("annotation/elferink2010_GRCg6a.txt")

cattle_map_chr <- split(cattle_map, cattle_map$chr)
chicken_map_chr <- split(chicken_map, chicken_map$chr)


get_interval_rate <- function(chr) {
  
  n_mar <- nrow(chr)
  chr$delta_physical <- c(NA, (chr$position_bp[2:n_mar] - chr$position_bp[-n_mar])/1e6)
  chr$delta_genetic <- c(NA, chr$position_cM[2:n_mar] - chr$position_cM[-n_mar])
  
  chr$r <- chr$delta_genetic/chr$delta_physical
  chr 
}


cattle_map_rate <- map_dfr(cattle_map_chr, get_interval_rate)
chicken_map_rate <- map_dfr(chicken_map_chr, get_interval_rate)

quantile(cattle_map_rate$r, c(0.25, 0.75), na.rm = TRUE)
quantile(chicken_map_rate$r, c(0.25, 0.75), na.rm = TRUE)