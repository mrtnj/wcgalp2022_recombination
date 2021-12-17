
library(assertthat)

make_simulation <- function(founders,
                            n_qtl_per_chr,
                            n_snps_per_chr,
                            mean_dd,
                            var_dd) {
  
  simparam <- SimParam$new(founders)
  simparam$setSexes("yes_sys")
  
  simparam$addTraitAD(nQtlPerChr = n_qtl_per_chr,
                      meanDD = mean_dd,
                      varDD = var_dd)
  simparam$addSnpChip(nSnpPerChr = n_snps_per_chr)
  
  simparam$setVarE(h2 = 0.3)
  
  pop <- newPop(founders,
                simParam = simparam)
  
  list(pop = pop,
       simparam = simparam)
}


run_breeding <- function(pop,
                         n_gen,
                         simparam) {
  
  generations <- vector(length = n_gen,
                        mode = "list")
  
  generations[[1]] <- pop
  
  for (gen_ix in 2:n_gen) {
    
    males <- generations[[gen_ix - 1]][generations[[gen_ix - 1]]@sex == "M"]
    females <- generations[[gen_ix - 1]][generations[[gen_ix - 1]]@sex == "F"]
    
    sires <- selectInd(pop = males,
                       nInd = 50,
                       simParam = simparam)
    dams <- selectInd(pop = females,
                      nInd = 100,
                      simParam = simparam)
    
    generations[[gen_ix]] <- randCross2(males = sires,
                                        females = dams,
                                        nProgeny = 10,
                                        nCrosses = 100,
                                        simParam = simparam)
    
  }
  
  
  generations  
  
}

run_genomic_selection <- function(pop,
                                  model,
                                  n_gen,
                                  simparam) {
  
  generations <- vector(length = n_gen,
                        mode = "list")
  
  generations[[1]] <- setEBV(pop,
                             solution = model,
                             simParam = simparam)
  
  for (gen_ix in 2:n_gen) {
    
    males <- generations[[gen_ix - 1]][generations[[gen_ix - 1]]@sex == "M"]
    females <- generations[[gen_ix - 1]][generations[[gen_ix - 1]]@sex == "F"]
    
    sires <- selectInd(pop = males,
                       use = "ebv",
                       nInd = 50,
                       simParam = simparam)
    dams <- selectInd(pop = females,
                      use = "ebv",
                      nInd = 100,
                      simParam = simparam)
    
    generations[[gen_ix]] <- randCross2(males = sires,
                                        females = dams,
                                        nProgeny = 10,
                                        nCrosses = 100,
                                        simParam = simparam)
    
    generations[[gen_ix]] <- setEBV(generations[[gen_ix]],
                                    solution = model,
                                    simParam = simparam)
    
  }
  
  
  generations  
  
}


## Peform stabilising selection based on a Gaussian fitness function
## with optimum at zero

select_stabilising <- function(pop,
                               prop) {
  
  fitness <- exp(-pop@pheno[,1]^2)
  
  rank <- order(fitness, decreasing = TRUE)
  
  n_to_select <- round(pop@nInd) * 0.15
  
  parent_ix <- rank[1:n_to_select]
  
  pop[parent_ix]
}

run_stabilising_selection <- function(pop,
                                      n_gen,
                                      simparam,
                                      prop) {
  
  generations <- vector(length = n_gen,
                        mode = "list")
  
  generations[[1]] <- pop
  
  for (gen_ix in 2:n_gen) {
    
    males <- generations[[gen_ix - 1]][generations[[gen_ix - 1]]@sex == "M"]
    females <- generations[[gen_ix - 1]][generations[[gen_ix - 1]]@sex == "F"]
    
    sires <- select_stabilising(pop = males,
                                prop = prop)
    dams <- select_stabilising(pop = females,
                               prop = prop)
    
    generations[[gen_ix]] <- randCross2(males = sires,
                                        females = dams,
                                        nProgeny = 10,
                                        nCrosses = 100,
                                        simParam = simparam)
    
  }
  
  generations  
  

}




get_gs_accuracy <- function(pop, simparam) {
  
  if (ncol(pop@ebv) > 0) {
    accuracy <- cor(bv(pop, simParam = simparam), ebv(pop))
  } else {
    accuracy <- NA_real_
  }
  
  accuracy
}


get_stats <- function(generations, simparam) {
  
  accuracy <- map_dbl(generations,
                      get_gs_accuracy,
                      simparam = simparam)
  
  tibble(gen = 1:length(generations),
         mean_g = unlist(lapply(generations, meanG)),
         var_g = unlist(lapply(generations, varG)),
         var_a = unlist(lapply(generations, varA, simParam = simparam)),
         accuracy = accuracy)
}



get_pheno_stats <- function(pop, simparam) {
  
  ebv <- rep(NA_real_, pop@nInd)
  
  if (ncol(ebv(pop) > 0)) {
    ebv <- ebv(pop)[, 1] 
  }
  
  tibble(pheno = pop@pheno[, 1],
         gv = gv(pop)[, 1],
         bv = bv(pop, simparam)[,1],
         ebv = ebv) 
  
}





##################################################################
## Linkage disequilibrium functions


## Get LD from haplotypes from pair of markers

get_ld <- function(haploA, haploB) {
  
  haploAB <- paste(haploA, haploB)
  
  pA <- sum(haploA)/length(haploA)
  pB <- sum(haploB)/length(haploB)
  
  pAB <- sum(haploAB == "1 1") / length(haploAB)
  
  D <- pAB - pA * pB
  
  r2 <- D^2 / (pA * (1 - pA) * pB * (1 - pB))
  
  r2
}


## Get LD from whole matrix of haplotypes

get_ld_matrix <- function(haplo) {
  
  fixed <- colSums(haplo) == 0 | colSums(haplo) == nrow(haplo)
  
  haplo <- haplo[,!fixed]
  
  n_qtl <- ncol(haplo)
  
  ld <- matrix(NA_real_,
               nrow = n_qtl,
               ncol = n_qtl)
  
  for (row_ix in 1:n_qtl) {
    for (col_ix in 1:n_qtl) {
      ld[row_ix, col_ix] <- get_ld(haplo[, row_ix], haplo[, col_ix])
    }
  }
  
  ld
}


## Get LD between QTN on a chromosome

get_chr_ld <- function(pop, simparam) {
  
  n_chr <- simparam$nChr
  
  chr_ld <- vector(length = n_chr,
                   mode = "list")
  
  for (chr_ix in 1:n_chr) {
    
    haplo <- pullQtlHaplo(pop,
                          chr = chr_ix,
                          simParam = simparam)
    
    chr_ld[[chr_ix]] <- get_ld_matrix(haplo)
    
  }
  chr_ld 
}




##################################################################

## Selection intensity functions


## Get the frequency of the allele with the lower genotypic value

get_qtl_low_allele_freq <- function(geno, simparam) {

  f <- colSums(geno)/2/nrow(geno)
  
  low_allele <- ifelse(simparam$traits[[1]]@addEff < 0,
                       1,
                       0)
  
  f_low <- ifelse(low_allele == 1,
                  f,
                  1 - f)
  
  f_low
}

## Get frequency of allele with the lower genotypic value from
## linkage equilibrium simulation

get_qtl_low_allele_freq_linkage_equilibrium <- function(f, simparam) {
  
  low_allele <- ifelse(simparam$traits[[1]]@addEff < 0,
                       1,
                       0)
  
  f_low <- ifelse(low_allele == 1,
                  f,
                  1 - f)
  
  f_low
}

## Calculate selection intensity i

calculate_selection_intensity <- function(prop) {
  dnorm(qnorm(1 - prop))/prop
}

assert_that(round(calculate_selection_intensity(0.2), 3) == 1.40)
assert_that(round(calculate_selection_intensity(0.1), 3) == 1.755)
assert_that(round(calculate_selection_intensity(0.02), 3) == 2.421)


## Calculate selection coefficient based in i, a and phenotypic sd

calculate_selection_coefficient <- function(prop_male,
                                            prop_female,
                                            a,
                                            pheno_sd) {
  i_male <- calculate_selection_intensity(prop_male)
  i_female <- calculate_selection_intensity(prop_female)
  i <- 0.5 * (i_male + i_female)
  
  i * 2 * a / pheno_sd
}




## Get allele frequency change, real and predicted, between two populations

compare_predicted_frequency_change <- function(f_gen1,
                                               f_gen2,
                                               prop_male,
                                               prop_female,
                                               var_p,
                                               simparam) {
  
  delta_f <- f_gen2 - f_gen1
  
  ## Since frequency is always for the lower allele, a is always positive
  a <- abs(simparam$traits[[1]]@addEff)
  
  pheno_sd <- sqrt(var_p)
  
  s <- calculate_selection_coefficient(prop_male,
                                       prop_female,
                                       a,
                                       pheno_sd)
  
  pred_delta_q <- -s * f_gen1^2 * (1 - f_gen1)
  
  tibble(locus = 1:length(delta_f),
         delta_f = delta_f,
         predicted_delta_f = pred_delta_q)
}



##################################################################

## Functions for using real genetic maps to adjust simulated maps



## Linear interpolation of genetic marker position

interpolate_cM <- function(marker_pos,
                           bp_before,
                           bp_after,
                           cM_before,
                           cM_after) {
  
  delta_bp <- bp_after - bp_before
  delta_cM <- cM_after - cM_before
  
  cM_per_bp <- delta_cM/delta_bp
  
  cM_before +  cM_per_bp * (marker_pos - bp_before)
}


## Adjust a simulated genetic map based on real genetic map; both
## real and simulated map need to be split by chromosome into a list

make_adjusted_map <- function(sim_map,
                              real_map) {
  
  new_map <- vector(length = length(sim_map),
                    mode = "list")
  
  n_chr <- length(new_map)
  
  for (chr_ix in 1:n_chr) {
    new_map[[chr_ix]] <- matrix(get_cM_positions_chr(sim_map[[chr_ix]],
                                                     real_map[[chr_ix]]),
                                ncol = 1) / 100
  }
  
  new_map
}


## Get positions in cM based on a real map for simulated map, which is taken
## to show the relative position in basepairs

get_cM_positions_chr <- function(sim_map_chr,
                                 real_map_chr) {
  
  chr_length_bp <- max(real_map_chr$position_bp)
  chr_length_morgan <- max(real_map_chr$position_cM) / 100
  n_markers <- nrow(sim_map_chr)
  
  cM <- numeric(n_markers)
  
  for (marker_ix in 1:n_markers) {
    
    marker_pos <- chr_length_bp * sim_map_chr[marker_ix] / chr_length_morgan
    
    ix_before <- max(which(real_map_chr$position_bp < marker_pos))
    ix_after <- min(which(real_map_chr$position_bp > marker_pos))
    
    ## First marker on chromosome
    if (ix_after == 1) {
      cM_position <- interpolate_cM(marker_pos,
                                    bp_before = 0,
                                    bp_after = real_map_chr$position_bp[ix_after],
                                    cM_before = 0,
                                    cM_after = real_map_chr$position_cM[ix_after])
    } else {
      ## Most cases
      cM_position <- interpolate_cM(marker_pos,
                                    bp_before = real_map_chr$position_bp[ix_before],
                                    bp_after = real_map_chr$position_bp[ix_after],
                                    cM_before = real_map_chr$position_cM[ix_before],
                                    cM_after = real_map_chr$position_cM[ix_after])
    }
    
    cM[marker_ix] <- cM_position
    
  }
  
  cM
}


