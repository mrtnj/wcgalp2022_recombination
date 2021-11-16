
library(assertthat)
library(purrr)


## LE simulation


## Draw genotypes for LE simulation from frequencies

draw_le_genotypes <- function(freq, n_ind) {

  geno <- map_dfc(freq,
                  function(p) rbinom(n_ind, size = 2, prob = p))

  geno
}



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


selected_stabilising_linkage_equilibrium <- function(geno, pheno, prop) {
  n_ind <- length(pheno)
  assert_that(nrow(geno) == n_ind)
  
  fitness <- exp(-pheno^2)
  
  selected_ix <- order(fitness,
                       decreasing = TRUE)[1:(n_ind * prop)]
  
  selected_geno <- geno[selected_ix,]
  
  selected_freq <- colSums(selected_geno)/nrow(selected_geno)/2
  
  selected_freq
}



run_selection_linkage_equilibrium <- function(founder_geno,
                                              a,
                                              d,
                                              Ve,
                                              n_ind,
                                              n_gen,
                                              selection_type,
                                              prop) {
  generations <- vector(length = n_gen,
                        mode = "list")
  
  ## Type of selection
  
  if (selection_type == "stabilising") {
    selection_function <- selected_stabilising_linkage_equilibrium
  } else if (selection_type == "directional") {
    selection_function <- selected_linkage_equilibrium 
  }
  
  ## Set up founder generation
  founder_p <- colSums(founder_geno)/2/nrow(founder_geno)
  
  drawn_founder_geno <- draw_le_genotypes(freq = founder_p,
                                          n_ind = n_ind)
  
  founder_genetic_values <- calculate_genetic_values(geno = drawn_founder_geno,
                                                     a = a,
                                                     d = d)
  founder_mean <- mean(founder_genetic_values)
  
  founder_phenotypes <- add_environmental_noise(founder_genetic_values - founder_mean, Ve)
  
  founder_selected_freq <- selection_function(geno = founder_geno,
                                              pheno = founder_phenotypes,
                                              prop = prop)
  
  generations[[1]] <- list(geno = drawn_founder_geno,
                           genetic_values = founder_genetic_values - founder_mean,
                           phenotypes = founder_phenotypes,
                           selected_freq = founder_selected_freq)
  
  for (gen_ix in 2:n_gen) {
    
    ## Draw genotypes
    geno <- draw_le_genotypes(freq = generations[[gen_ix - 1]]$selected_freq, 
                              n_ind = n_ind)
    
    ## Genetic values
    genetic_values <- calculate_genetic_values(geno = geno, 
                                               a = a,
                                               d = d) - founder_mean
    
    ## Phenotypes
    phenotypes <- add_environmental_noise(genetic_value = genetic_values,
                                          Ve = Ve)
    ## Selection
    selected_freq <- selection_function(geno = geno,
                                        pheno = phenotypes,
                                        prop = prop)
    
    generations[[gen_ix]] <- list(geno = geno,
                                  genetic_values = genetic_values,
                                  phenotypes = phenotypes,
                                  selected_freq = selected_freq)
    
  }
  
  generations
  
}


