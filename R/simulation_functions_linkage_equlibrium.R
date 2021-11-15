
library(purrr)


## LE simulation


## Draw genotypes for LE simulation from frequencies

draw_le_genotypes <- function(freq, n_ind) {

  geno <- map_dfc(freq,
                  function(p) rbinom(n_ind, size = 2, prob = p))

  geno
}


geno <- draw_le_genotypes(founder_p,
                          n_ind)

## Get genetic values

calculate_genetic_values <- function(geno, a, d) {

  geno_additive <- geno - 1
  
  genetic_value <- numeric(n_ind)
  
  for (locus_ix in 1:n_loc) {
    
    locus_genetic_value <- geno_additive[, locus_ix] * a[locus_ix] +
      as.numeric(geno_additive[, locus_ix] == 0) * d[locus_ix]
    
    genetic_value <- genetic_value + locus_genetic_value
    
  }
  
  genetic_value
}


## Add environmental variation to get phenotypic values

add_environmental_noise <- function(genetic_value,
                                    Ve) {
  
  pheno <- genetic_value + rnorm(n_ind, 0, sqrt(Ve))
  
  pheno

}


## Selection function for the linkage equilibrium model; selects
## individuals and returns new allele frequencies

selected_linkage_equilibrium <- function(geno, pheno, prop) {
  n_ind <- length(pheno)
  assert_that(nrow(geno) == n_ind)
  
  selected_ix <- order(pheno,
                       decreasing = TRUE)[1:(n_ind * prop)]
  
  selected_geno <- geno[selected_ix,]
  
  selected_freq <- colSums(selected_geno)/nrow(selected_geno)/2
  
  selected_freq
}


sim_linkage_equilibrium <- function(founder_geno, a, d, Ve, n_ind, n_gen, prop) {

  geno <- vector(length = n_gen,
                 mode = "list")
  
  genetic_values <- vector(length = n_gen,
                           mode = "list")
  
  phenotypes <- vector(length = n_gen,
                       mode = "list")
  
  selected_freq <- vector(length = n_gen,
                          mode = "list")
  
  geno[[1]] <- founder_geno
  genetic_values[[1]] <- calculate_genetic_values(geno = founder_geno, 
                                                  a = a,
                                                  d = d)
  founder_mean <- mean(genetic_values[[1]])
  genetic_values[[1]] <- genetic_values[[1]] - founder_mean
  
  phenotypes[[1]] <- add_environmental_noise(genetic_values[[1]], Ve)
  selected_freq[[1]] <- selected_linkage_equilibrium(geno = founder_geno,
                                                     pheno = phenotypes[[1]],
                                                     prop = prop)
  
  
  
  for (gen_ix in 2:n_gen) {
    
    ## Draw genotypes
    geno[[gen_ix]] <- draw_le_genotypes(freq = selected_freq[[gen_ix - 1]], 
                                        n_ind = n_ind)
    
    ## Genetic values
    genetic_values[[gen_ix]] <- calculate_genetic_values(geno = geno[[gen_ix]], 
                                                         a = a,
                                                         d = d) - founder_mean
    
    ## Phenotypes
    phenotypes[[gen_ix]] <- add_environmental_noise(genetic_value = genetic_values[[gen_ix]],
                                                    Ve = Ve)
    ## Selection
    selected_freq[[gen_ix]] <- selected_linkage_equilibrium(geno = geno[[gen_ix]],
                                                            pheno = phenotypes[[gen_ix]],
                                                            prop = prop)
    
  }
  
  list(freq,
       geno,
       genetic_values,
       phenotypes,
       selected_ix)

}


## Test run

founder_geno <- pullQtlGeno(pop,
                            simParam = simparam)

founder_p <- colSums(founder_geno)/nrow(founder_geno)/2

a <- simparam$traits[[1]]@addEff
d <- simparam$traits[[1]]@domEff

Ve <- simparam$varE

n_ind <- 1000
n_loc <- length(founder_p)


s <- sim_linkage_equilibrium(founder_geno, a, d, Ve, n_ind = 1000, n_gen = 20, prop = 0.15)

meang <- map_dbl(s[[3]], mean)

plot(meang)
