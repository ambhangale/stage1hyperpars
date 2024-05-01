# Aditi M. Bhangale
# Last updated: 1 May 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the
# multivariate social relations model
# Simulations 1 to 4

# Snellius status

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/sim_code")
# getwd()

###### SIMULATIONS 1 AND 2

# simConds for analTypes with precisions
simCondsA <- rbind(expand.grid(nSamps = 1000, n = c(6,8,10,20), G = c(10,25),
                               analType = "prophetic", precision = c(0.05, 0.1, 0.2), sim = "sim1"), # sim1
                   expand.grid(nSamps = 1000, n = c(6,8,10,20), G = c(10,25),
                               analType = c("thoughtful", "ANOVA", "FIML"), precision = 0.1, sim = "sim2")) # sim2

simCondsA$wallTime <- c("06:00:00", "10:00:00", "14:00:00", "5-00:00:00", # prophetic_0.05_G10
                        "08:00:00", "12:00:00", "16:00:00", "5-00:00:00", # prophetic_0.05_G25
                        "08:00:00", "12:00:00", "14:00:00", "5-00:00:00", # prophetic_0.1_G10
                        "09:00:00", "13:00:00", "16:00:00", "5-00:00:00", # prophetic_0.1_G25
                        "10:00:00", "13:00:00", "16:00:00", "5-00:00:00", # prophetic_0.2_G10
                        "13:00:00", "16:00:00", "19:00:00", "5-00:00:00", # prophetic_0.2_G25
                        "07:00:00", "09:00:00", "11:00:00", "5-00:00:00", # thoughtful_0.1_G10
                        "08:00:00", "12:00:00", "22:00:00", "5-00:00:00", # thoughtful_0.1_G25
                        "12:00:00", "14:00:00", "21:00:00", "5-00:00:00", # ANOVA_0.1_G10
                        "12:00:00", "16:00:00", "23:00:00", "5-00:00:00", # ANOVA_0.1_G25
                        "2-08:00:00", "2-12:00:00", "2-22:00:00", "5-00:00:00", # FIML_0.1_G10
                        "2-08:00:00", "5-00:00:00", "5-00:00:00", "5-00:00:00" # FIML_0.1_G25
)

# NA = not applicable (code not run yet)
# IQ = in queue
# S_T = sent to Terrence
# R_A = running on Aditi's account
# R_T = running on Terrence's account
# CD_A = completed on Aditi's account
# CD_T = completed on Terrence's account
# F_A = failed on A's Snellius
# R_WA = running on A's windows
# R_WT = running on T's windows
# CD_WA = completed on A's windows
# CD_WT = completed on T's windows

simCondsA$snelliusStatus <- c("CD_A", "CD_A", "CD_A", "CD_T", # prophetic_0.05_G10
                              "CD_A", "CD_A", "CD_A", "CD_T", # prophetic_0.05_G25
                              "CD_A", "CD_A", "CD_A", "CD_T", # prophetic_0.1_G10
                              "CD_A", "CD_A", "CD_A", "CD_T", # prophetic_0.1_G25
                              "CD_A", "CD_A", "CD_A", "CD_T", # prophetic_0.2_G10
                              "CD_A", "CD_A", "CD_A", "CD_T", # prophetic_0.2_G25
                              "CD_A", "CD_A", "CD_A", "CD_A", # thoughtful_0.1_G10
                              "CD_A", "CD_A", "CD_A", "CD_A", # thoughtful_0.1_G25
                              "CD_A", "CD_A", "CD_A", "CD_A", # ANOVA_0.1_G10
                              "CD_A", "CD_A", "CD_A", "CD_A", # ANOVA_0.1_G25
                              "CD_A", "CD_A", "CD_A", "CD_A", # FIML_0.1_G10
                              "CD_A", "CD_T", "CD_A", "NA" # FIML_0.1_G25
)

simCondsA$savedResult <- c("Y", "Y", "Y", "Y", # prophetic_0.05_G10
                           "Y", "Y", "Y", "Y", # prophetic_0.05_G25
                           "Y", "Y", "Y", "Y", # prophetic_0.1_G10
                           "Y", "Y", "Y", "Y", # prophetic_0.1_G25
                           "Y", "Y", "Y", "Y", # prophetic_0.2_G10
                           "Y", "Y", "Y", "Y", # prophetic_0.2_G25
                           "Y", "Y", "Y", "Y", # thoughtful_0.1_G10
                           "Y", "Y", "Y", "Y", # thoughtful_0.1_G25
                           "Y", "Y", "Y", "Y", # ANOVA_0.1_G10
                           "Y", "Y", "Y", "Y", # ANOVA_0.1_G25
                           "Y", "Y", "Y", "Y", # FIML_0.1_G10
                           "Y", "Y", "Y", "NA" # FIML_0.1_G25
)

# simConds for analTypes without precisions
simCondsB <- expand.grid(nSamps = 1000, n = c(6,8,10,20), G = c(10,25),
                         analType = c("default", "FIML1S"), sim = "sim1")

simCondsB$wallTime <- c("12:00:00", "14:00:00", "1-12:00:00", "2-12:00:00", # default_G10
                        "20:00:00", "1-01:00:00", "1-18:00:00", "2-18:00:00", # default_G25
                        "03:00:00", "1-08:00:00", "1-08:00:00", "5-00:00:00", # FIML1S_G10
                        "05:00:00", "2-18:00:00", "3-20:00:00", "5-00:00:00" # FIML1S_G25
)

simCondsB$snelliusStatus <- c("CD_A", "CD_A", "CD_A", "CD_A", # default_G10
                                "CD_A", "CD_A", "CD_A", "CD_A", # default_G25
                              "CD_A", "CD_A", "CD_A", "CD_WT", # FIML1S_G10
                              "CD_A", "CD_A", "CD_A", "CD_WA" # FIML1S_G25
)

simCondsB$savedResult <- c("Y", "Y", "Y", "Y", # default_G10
                           "Y", "Y", "Y", "Y", # default_G25
                           "Y", "Y", "Y", "Y", # FIML1S_G10
                           "Y", "Y", "Y", "Y" # FIML1S_G25
)

# simCondsA[simCondsA$snelliusStatus == "NA",]
# simCondsB[simCondsB$snelliusStatus == "NA",]

# rbind(data.frame(nSamps = 1000, n = c(6, 8, 10, 20, 8, 10, 20), 
#                  G = c(rep(10,4), rep(25, 3)), analType = "FIML", precision = 0.1), 
#       data.frame(nSamps = 1000, n = c(10, 20, 8, 10, 20), 
#                  G = c(rep(10,2), rep(25, 3)), analType = "ANOVA", precision = 0.1))


###### SIMULATION 3
simConds_sim3 <- expand.grid(nSamps = 1000, n = c(6,8), G = 50,
                             analType = c("default", "thoughtful"), sim = "sim3")

simConds_sim3$wallTime <- c("5-00:00:00", "5-00:00:00",
                            "5-00:00:00", "5-00:00:00")

simConds_sim3$snelliusStatus <- c("R_A", "R_A",
                                  "R_A", "R_A")

simConds_sim3$savedResult <- c("NA", "NA",
                               "NA", "NA")
