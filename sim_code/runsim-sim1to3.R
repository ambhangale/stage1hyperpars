## Aditi M. Bhangale
## Last updated: 4 April 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model
# runsim for simulations 1 and 2

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/sim_code")
# getwd()

source("functions-sim1to3.R")

library(doSNOW)

## conditions
# 8 (analType/prior) x 4 (n; group size) x 2 (G; number of groups)
## create separate runsim functions for ogsat and s1sat 

# common conditions for ogsat() and s1sat()
n <- c(6, 8, 10, 20)
G <- c(10, 25)
MCSampID <- 1:5

# specific conditions for s1sat()
priorType <- c("default", "thoughtful", "prophetic", "ANOVA", "FIML")

# all conditions
og_grid <- expand.grid(MCSampID = MCSampID, n = n, G = G)
og_grid$row_num <- 1:nrow(og_grid)

s1_grid <- rbind(expand.grid(MCSampID = MCSampID, n = n, G = G,
                       priorType = priorType, precision = 0.1),
                 expand.grid(MCSampID = MCSampID, n = n, G = G,
                             priorType = "prophetic", precision = c(0.05, 0.2)))
s1_grid$row_num <- 1:nrow(s1_grid)


# prepare parallel processing
# nClus <- 32
nClus <- parallel::detectCores() - 1
cl <- makeCluster(nClus)
registerDoSNOW(cl)

# run simulation
ogResult <- foreach(row_num = 1:nrow(og_grid),
                    .packages = c("mnormt", "parallel", "portableParallelSeeds",
                                  "TripleR", "srm", "car", "lavaan.srm")) %dopar% {
  
  out <- try(ogsat(MCSampID = og_grid[row_num, ]$MCSampID, og_grid[row_num, ]$n, 
               og_grid[row_num, ]$G), silent = T)
  if(inherits(out, "try-error")) out <- NULL
  
  return(out)
} #TODO try this again after removing TripleR and lavaan.srm because you don't need these for og

s1Result <- foreach(row_num = 1:nrow(s1_grid),
                    .packages = c("mnormt", "parallel", "portableParallelSeeds",
                                  "TripleR", "srm", "car", "lavaan.srm")) %dopar% {
  
  out <- try(s1sat(MCSampID = s1_grid[row_num, ]$MCSampID, 
               n = s1_grid[row_num, ]$n, G = s1_grid[row_num, ]$G,
               priorType = s1_grid[row_num, ]$priorType, 
               precision = s1_grid[row_num, ]$precision), silent = T)
  if(inherits(out, "try-error")) out <- NULL
  
  return(out)
}

# close cluster
stopCluster(cl)

# save results
saveRDS(ogResult, "results_1S-FIML.rds")
saveRDS(s1Result, "results_MCMC.rds")
