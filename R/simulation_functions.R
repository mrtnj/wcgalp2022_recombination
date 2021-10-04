

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
