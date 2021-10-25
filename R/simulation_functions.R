

make_simulation <- function(founders,
                            n_qtl_per_chr,
                            mean_dd,
                            var_dd) {
  
  simparam <- SimParam$new(founders)
  simparam$setSexes("yes_sys")
  
  simparam$addTraitAD(nQtlPerChr = n_qtl_per_chr,
                      meanDD = mean_dd,
                      varDD = var_dd)
  
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



get_stats <- function(result) {
  tibble(gen = 1:length(result),
         mean_g = unlist(lapply(result, meanG)),
         var_g = unlist(lapply(result, varG)))
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

get_qtl_low_allele_freq <- function(pop, simparam) {
  geno <- pullQtlGeno(pop,
                      simParam = simparam)
  f <- colSums(geno)/2/nrow(geno)
  
  low_allele <- ifelse(simparam$traits[[1]]@addEff < 0,
                       1,
                       0)
  
  f_low <- ifelse(low_allele == 1,
                  f,
                  1 - f)
}

## Calculate selection intensity i

calculate_selection_intensity <- function(prop) {
  dnorm(qnorm(1 - prop))/prop
}


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

compare_predicted_frequency_change <- function(pop1,
                                               pop2,
                                               prop_male,
                                               prop_female,
                                               simparam) {
  
  f_gen1 <- get_qtl_low_allele_freq(pop1, simparam)
  f_gen2 <- get_qtl_low_allele_freq(pop2, simparam)
  
  delta_f <- f_gen2 - f_gen1
  
  ## Since frequency is always for the lower allele, a is always positive
  a <- abs(simparam$traits[[1]]@addEff)
  
  pheno_sd <- as.numeric(sqrt(varP(generations[[gen1]])))
  
  s <- calculate_selection_coefficient(prop_male,
                                       prop_female,
                                       a,
                                       pheno_sd)
  
  pred_delta_q <- -s * f_gen1^2 * (1 - f_gen1)
  
  tibble(locus = 1:length(delta_f),
         delta_f = delta_f,
         predicted_delta_f = pred_delta_q)
}
