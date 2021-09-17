## Simulate chromosome length variation

library(AlphaSimR)

source("R/simulation_functions.R")

  
founders_short <- runMacs2(nInd = 1000,
                           nChr = 10)

founders_long <- runMacs2(nInd = 1000,
                          nChr = 10,
                          genLen = 2)

simulation_short <- make_simulation(founders_short)
simulation_long <- make_simulation(founders_long)

pop_short <- simulation_short$pop
simparam_short <- simulation_short$simparam

pop_long <- simulation_long$pop
simparam_long <- simulation_long$simparam

result_short <- run_breeding(pop_short,
                             10,
                             simparam_short)

result_long <- run_breeding(pop_long,
                            10,
                            simparam_long)
