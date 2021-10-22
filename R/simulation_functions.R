

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