# Aditi M. Bhangale
# Last updated: 11 April 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the
# multivariate social relations model
# Simulations 1 to 3

# make runsim and shell files for simulation

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/sim_code")
# getwd()

source("functions_sim1to3.R")

# simConds for analTypes with precisions
simCondsA <- rbind(expand.grid(nSamps = 1000, n = c(6,8,10,20), G = c(10,25),
              analType = "prophetic", precision = c(0.05, 0.1, 0.2), sim = "sim1"), # sim1
  expand.grid(nSamps = 1000, n = c(6,8,10,20), G = c(10,25),
              analType = c("thoughtful", "ANOVA", "FIML"), precision = 0.1, sim = "sim2")) # sim2

simCondsA$wallTime <- c("06:00:00", "10:00:00", "14:00:00", "1-02:00:00", # prophetic_0.05_G10
                          "08:00:00", "12:00:00", "16:00:00", "1-12:00:00", # prophetic_0.05_G25
                          "08:00:00", "12:00:00", "14:00:00", "20:00:00", # prophetic_0.1_G10
                          "09:00:00", "13:00:00", "16:00:00", "3-23:59:59", # prophetic_0.1_G25
                          "10:00:00", "13:00:00", "16:00:00", "1-01:00:00", # prophetic_0.2_G10
                          "13:00:00", "16:00:00", "19:00:00", "1-23:59:59", # prophetic_0.2_G25
                          "07:00:00", "09:00:00", "11:00:00", "22:00:00", # thoughtful_0.1_G10
                          "08:00:00", "12:00:00", "15:00:00", "1-23:59:59", # thoughtful_0.1_G25
                          "12:00:00", "14:00:00", "18:00:00", "23:00:00", # ANOVA_0.1_G10
                          "12:00:00", "16:00:00", "20:00:00", "3-23:59:59", # ANOVA_0.1_G25
                          "2-08:00:00", "2-12:00:00", "2-22:00:00", "3-01:00:00", # FIML_0.1_G10
                          "2-08:00:00", "2-12:00:00", "3-02:00:00", "3-20:00:00" # FIML_0.1_G25
)

# simConds for analTypes without precisions
simCondsB <- expand.grid(nSamps = 1000, n = c(6,8,10,20), G = c(10,25),
            analType = c("default", "FIML1S"), sim = "sim1")
  
simCondsB$wallTime <- c("12:00:00", "14:00:00", "1-12:00:00", "2-12:00:00", # default_G10
                "20:00:00", "1-01:00:00", "1-18:00:00", "2-18:00:00", # default_G25
                "03:00:00", "12:00:00", "16:00:00", "20:00:00", # FIML1S_G10
                "05:00:00", "14:00:00", "18:00:00", "22:00:00" # FIML1S_G25
                )

# generate files  
for(i in 1:nrow(simCondsA)) {
  makeRunsim(nSamps = simCondsA[i, ]$nSamps, n = simCondsA[i, ]$n, 
             G = simCondsA[i, ]$G, analType = simCondsA[i, ]$analType, 
             precision = simCondsA[i, ]$precision, sim = simCondsA[i, ]$sim)
  makeShSnellius(n = simCondsA[i, ]$n, 
                 G = simCondsA[i, ]$G, analType = simCondsA[i, ]$analType, 
                 precision = simCondsA[i, ]$precision, sim = simCondsA[i, ]$sim,
                 wallTime = simCondsA[i, ]$wallTime)
}

for(i in 1:nrow(simCondsB)) {
  makeRunsim(nSamps = simCondsB[i, ]$nSamps, n = simCondsB[i, ]$n, 
             G = simCondsB[i, ]$G, analType = simCondsB[i, ]$analType, 
             precision = simCondsB[i, ]$precision, sim = simCondsB[i, ]$sim)
  makeShSnellius(n = simCondsB[i, ]$n, 
                 G = simCondsB[i, ]$G, analType = simCondsB[i, ]$analType, 
                 precision = simCondsB[i, ]$precision, sim = simCondsB[i, ]$sim,
                 wallTime = simCondsB[i, ]$wallTime)
}

