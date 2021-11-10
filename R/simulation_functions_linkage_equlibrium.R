## LE simulation


founder_geno <- pullQtlGeno(pop,
                            simParam = simparam)

founder_p <- colSums(founder_geno)/nrow(founder_geno)/2

a <- simparam$traits[[1]]@addEff
d <- simparam$traits[[1]]@domEff

Ve <- simparam$varE

n_ind <- 1000
n_loc <- length(founder_p)


geno <- map_dfc(founder_p,
                function(p) rbinom(n_ind, size = 2, prob = p))


## Get genetic values

geno_additive <- geno - 1

genetic_value <- numeric(n_ind)

for (locus_ix in 1:n_loc) {
  
  locus_genetic_value <- geno_additive[, locus_ix] * a[locus_ix] +
    as.numeric(geno_additive[, locus_ix] == 0) * d[locus_ix]
  
  genetic_value <- genetic_value + locus_genetic_value
  
}

genetic_value <- genetic_value - mean(genetic_value)


## Get phenotypes and select

## varG/varP too small compared to ASR?

pheno <- genetic_value + rnorm(n_ind, 0, Ve)

prop <- 0.15

selected_ix <- order(pheno, decreasing = TRUE)[1:(n_ind * prop)]


parent_geno <- geno[selected_ix, ]

p <- colSums(parent_geno)/n_ind/2
