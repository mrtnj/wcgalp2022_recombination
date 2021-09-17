

make_simulation <- function(founders) {
  simparam <- SimParam$new(founders)
  simparam$setSexes("yes_sys")
  
  simparam$addTraitA(nQtlPerChr = 1000)
  
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