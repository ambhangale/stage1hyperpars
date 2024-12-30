## Aditi M. Bhangale
## Last updated: 30 December 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model

# Data processing

# rm(list = ls())

##README MCSampID 198 in FIML1S 6-10 is not in the output. This is probably because
## that sample did not converge

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/results_and_postprocessing")
# getwd()

# function 1: extract and load results----
extract_results <- function(analType = NULL, precision = NULL, sim) {
  # files <- grep("results", grep(sim, dir(), value = T), value = T)
  files <- if (!missing(analType) && !missing(precision)) {
    grep("results", grep(analType, grep(precision, grep(sim, dir(), value = T), 
                                        value = T), value = T), value = T)
  } else if (missing(precision) && !missing(analType)) {
    grep("results", grep(analType, grep(sim, dir(), value = T), value = T), value = T)
  } else {
    grep("results", grep(sim, dir(), value = T), value = T)
  }
  
  cov_list <- list()
  cor_list <- list()
  SD_list <- list()
  mPSRF_list <- list()
  all_list <- list()
  
  for (i in 1:length(files)) {
    file_results <- readRDS(files[i])
  
    # file_cov_list <- lapply(1:length(file_results), function(MCSampID) file_results[[MCSampID]]$cov)
    file_cov_result <- do.call("rbind", lapply(1:length(file_results), 
                                               function(MCSampID) file_results[[MCSampID]]$cov))
    # file_cor_list <- lapply(1:length(file_results), function(MCSampID) file_results[[MCSampID]]$cor)
    file_cor_result <- do.call("rbind", lapply(1:length(file_results), 
                                               function(MCSampID) file_results[[MCSampID]]$cor))
    # file_SD_list <- lapply(1:length(file_results), function(MCSampID) file_results[[MCSampID]]$SD)
    file_SD_result <- do.call("rbind", lapply(1:length(file_results), 
                                              function(MCSampID) file_results[[MCSampID]]$SD))
    
    if (isFALSE(grepl("FIML1S", files[i]))) {
      # file_mPSRF_list <- lapply(1:length(file_results), function(MCSampID) file_results[[MCSampID]]$mPSRF)
      file_mPSRF_result <- data.frame(do.call("rbind", lapply(1:length(file_results), 
                                                              function(MCSampID) file_results[[MCSampID]]$mPSRF)))
      
    } else {
      file_mPSRF_result <- NULL
    }
    cov_list[[i]] <- file_cov_result
    cor_list[[i]] <- file_cor_result
    SD_list[[i]] <- file_SD_result
    if (!is.null(file_mPSRF_result)) mPSRF_list[[i]] <- file_mPSRF_result
  }
  
  cov_result <- do.call("rbind", cov_list)
  cor_result <- do.call("rbind", cor_list)
  SD_result <- do.call("rbind", SD_list)
  mPSRF_result <- do.call("rbind", mPSRF_list)
  all_list <- list(cov_result = cov_result, cor_result = cor_result, 
                   SD_result = SD_result, mPSRF_result = mPSRF_result)
  
}

# foo <- extract_results(sim = "sim1")

#----

# function 2: convert results into long format (single sample)----
make_long_single <- function(dat, estType) {
  resType <- c("cov_result", "cor_result", "SD_result")
  common_cols <- c("par_names", "MCSampID", "n", "G", "condition", "level", "analType", "iter",
                   "RunTime")
  newcolNames_cov <- c(common_cols, "pop.cov", "cov.est", "cov.SE", "cov.low", "cov.up", 
                       "cov.MCMC_neff", "cov.MCMC_Rhat", "estType")
  newcolNames_cor <- c(common_cols, "pop.cor", "cor.est", "cor.SE", "cor.low", "cor.up", 
                       "prior1", "prior2", "cor.MCMC_neff", "cor.MCMC_Rhat", "estType")
  newcolNames_SD <- c(common_cols, "pop.SD", "SD.est", "SD.SE", "SD.low", "SD.up", 
                      "prior1", "prior2", "SD.MCMC_neff", "SD.MCMC_Rhat", "estType")
  
  all_list <- list()
  
  if (estType == "EAP") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_EAP <- dat[[i]][, c(common_cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", 
                                "Ecov.up", "Ecov.n_eff", "Ecov.Rhat")]
        if ("1S-FIML" %in% unique(cov_EAP$analType)) {
          cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_EAP$estType <- "EAP" # add estType
        colnames(cov_EAP) <- newcolNames_cov
        
      } else if (i == "cor_result") {
        cor_EAP <- dat[[i]][, c(common_cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low",
                                "Ecor.up", "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        if ("1S-FIML" %in% unique(cor_EAP$analType)) {
          cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_EAP$estType <- "EAP" # add estType
        colnames(cor_EAP) <- newcolNames_cor
        
      } else if (i == "SD_result") {
        sd_EAP <- dat[[i]][, c(common_cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                               "prior1", "prior2", "Esd.n_eff", "Esd.Rhat")]
        if ("1S-FIML" %in% unique(sd_EAP$analType)) { 
          sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_EAP$estType <- "EAP" # add estType
        colnames(sd_EAP) <- newcolNames_SD
        
      }
    }
    all_list <- list(cov_result = cov_EAP, cor_result = cor_EAP, 
                     SD_result = sd_EAP, mPSRF_result = dat$mPSRF_result)
    
  } else if (estType == "MAP") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_MAP <- dat[[i]][, c(common_cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up", 
                                "Ecov.n_eff", "Ecov.Rhat")]
        cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column, then rearrange columns
        cov_MAP <- cov_MAP[, c(common_cols, "pop.cov", "Mcov", "Mcov.SE", 
                               "Mcov.low", "Mcov.up", "Ecov.n_eff", "Ecov.Rhat")]
        if ("1S-FIML" %in% unique(cov_MAP$analType)) {
          cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_MAP$estType <- "MAP" # add estType
        colnames(cov_MAP) <- newcolNames_cov
        
      } else if (i == "cor_result") {
        cor_MAP <- dat[[i]][, c(common_cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
        cor_MAP <- cor_MAP[, c(common_cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                               "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        if ("1S-FIML" %in% unique(cor_MAP$analType)) {
          cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_MAP$estType <- "MAP" # add estType
        colnames(cor_MAP) <- newcolNames_cor
        
      } else if (i == "SD_result") {
        sd_MAP <- dat[[i]][, c(common_cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", 
                               "prior2", "Esd.n_eff", "Esd.Rhat")]
        sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
        sd_MAP <- sd_MAP[, c(common_cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                             "prior1", "prior2", "Esd.n_eff", "Esd.Rhat")]
        if ("1S-FIML" %in% unique(sd_MAP$analType)) {
          sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_MAP$estType <- "MAP" # add estType
        colnames(sd_MAP) <- newcolNames_SD
      } 
    }
    all_list <- list(cov_result = cov_MAP, cor_result = cor_MAP, 
                     SD_result = sd_MAP, mPSRF_result = dat$mPSRF_result)
    
  } else if (estType == "allMCMC") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_EAP <- dat[[i]][, c(common_cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", 
                                "Ecov.up", "Ecov.n_eff", "Ecov.Rhat")]
        if ("1S-FIML" %in% unique(cov_EAP$analType)) {
          cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_EAP$estType <- "EAP" # add estType
        colnames(cov_EAP) <- newcolNames_cov
        
        cov_MAP <- dat[[i]][, c(common_cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up", 
                                "Ecov.n_eff", "Ecov.Rhat")]
        cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column, then rearrange columns
        cov_MAP <- cov_MAP[, c(common_cols, "pop.cov", "Mcov", "Mcov.SE", 
                               "Mcov.low", "Mcov.up", "Ecov.n_eff", "Ecov.Rhat")]
        if ("1S-FIML" %in% unique(cov_MAP$analType)) {
          cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_MAP$estType <- "MAP" # add estType
        colnames(cov_MAP) <- newcolNames_cov
        
        cov_result <- rbind(cov_EAP, cov_MAP)
        
      } else if (i == "cor_result") {
        cor_EAP <- dat[[i]][, c(common_cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low", "Ecor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        if ("1S-FIML" %in% unique(cor_EAP$analType)) {
          cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_EAP$estType <- "EAP" # add estType
        colnames(cor_EAP) <- newcolNames_cor
        
        cor_MAP <- dat[[i]][, c(common_cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
        cor_MAP <- cor_MAP[, c(common_cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                               "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        if ("1S-FIML" %in% unique(cor_MAP$analType)) {
          cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_MAP$estType <- "MAP" # add estType
        colnames(cor_MAP) <- newcolNames_cor
        
        cor_result <- rbind(cor_EAP, cor_MAP)
        
      } else if (i == "SD_result") {
        sd_EAP <- dat[[i]][, c(common_cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                               "prior1", "prior2", "Esd.n_eff", "Esd.Rhat")]
        if ("1S-FIML" %in% unique(sd_EAP$analType)) { 
          sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_EAP$estType <- "EAP" # add estType
        colnames(sd_EAP) <- newcolNames_SD
        
        sd_MAP <- dat[[i]][, c(common_cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", 
                               "prior2", "Esd.n_eff", "Esd.Rhat")]
        sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
        sd_MAP <- sd_MAP[, c(common_cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                             "prior1", "prior2", "Esd.n_eff", "Esd.Rhat")]
        if ("1S-FIML" %in% unique(sd_MAP$analType)) {
          sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_MAP$estType <- "MAP" # add estType
        colnames(sd_MAP) <- newcolNames_SD
        
        SD_result <- rbind(sd_EAP, sd_MAP)
      }
    }
    all_list <- list(cov_result = cov_result, cor_result = cor_result, 
                     SD_result = SD_result, mPSRF_result = dat$mPSRF_result)
    
  } else if (estType == "FIML1S") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_ogest <- dat[[i]][, c(common_cols, "pop.cov", "ogcov", "ogcov.SE", 
                                  "ogcov.low", "ogcov.up")]
        cov_ogest <- cov_ogest[cov_ogest$analType == "1S-FIML", ] # remove redundant rows
        cov_ogest$cov.MCMC_neff <- NA; cov_ogest$cov.MCMC_Rhat <- NA
        cov_ogest$estType <- "ogFIML" # add estType
        colnames(cov_ogest) <- newcolNames_cov
        
      } else if (i == "cor_result") {
        cor_ogest <- dat[[i]][, c(common_cols, "pop.cor", "ogcor", "ogcor.SE", 
                                  "ogcor.low", "ogcor.up")]
        cor_ogest$prior1 <- NA; cor_ogest$prior2 <- NA 
        cor_ogest$cor.MCMC_neff <- NA; cor_ogest$cor.MCMC_Rhat <- NA # create dummy prior columns
        cor_ogest <- cor_ogest[cor_ogest$analType == "1S-FIML", ] # remove redundant rows
        cor_ogest$estType <- "ogFIML" # add estType
        colnames(cor_ogest) <- newcolNames_cor
        
      } else if (i == "SD_result") {
        sd_ogest <- dat[[i]][, c(common_cols, "pop.SD", "ogsd", "ogsd.SE", 
                                 "ogsd.low", "ogsd.up")]
        sd_ogest$prior1 <- NA; sd_ogest$prior2 <- NA 
        sd_ogest$SD.MCMC_neff <- NA; sd_ogest$SD.MCMC_Rhat <- NA # create dummy prior columns
        sd_ogest <- sd_ogest[sd_ogest$analType == "1S-FIML", ] # remove redundant rows
        sd_ogest$estType <- "ogFIML" # add estType
        colnames(sd_ogest) <- newcolNames_SD
      }
    }
    all_list <- list(cov_result = cov_ogest, cor_result = cor_ogest, 
                     SD_result = sd_ogest, mPSRF_result = NULL)
    
  } else if (estType == "all") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_EAP <- dat[[i]][, c(common_cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", 
                                "Ecov.up", "Ecov.n_eff", "Ecov.Rhat")]
        if ("1S-FIML" %in% unique(cov_EAP$analType)) {
          cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_EAP$estType <- "EAP" # add estType
        colnames(cov_EAP) <- newcolNames_cov
        
        cov_MAP <- dat[[i]][, c(common_cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up", 
                                "Ecov.n_eff", "Ecov.Rhat")]
        cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column, then rearrange columns
        cov_MAP <- cov_MAP[, c(common_cols, "pop.cov", "Mcov", "Mcov.SE", 
                               "Mcov.low", "Mcov.up", "Ecov.n_eff", "Ecov.Rhat")]
        if ("1S-FIML" %in% unique(cov_MAP$analType)) {
          cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_MAP$estType <- "MAP" # add estType
        colnames(cov_MAP) <- newcolNames_cov
        
        cov_ogest <- dat[[i]][, c(common_cols, "pop.cov", "ogcov", "ogcov.SE", 
                                  "ogcov.low", "ogcov.up")]
        cov_ogest <- cov_ogest[cov_ogest$analType == "1S-FIML", ] # remove redundant rows
        cov_ogest$cov.MCMC_neff <- NA; cov_ogest$cov.MCMC_Rhat <- NA
        cov_ogest$estType <- "ogFIML" # add estType
        colnames(cov_ogest) <- newcolNames_cov
        
        cov_result <- rbind(cov_EAP, cov_MAP, cov_ogest)
        
      } else if (i == "cor_result") {
        cor_EAP <- dat[[i]][, c(common_cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low", "Ecor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        if ("1S-FIML" %in% unique(cor_EAP$analType)) {
          cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_EAP$estType <- "EAP" # add estType
        colnames(cor_EAP) <- newcolNames_cor
        
        cor_MAP <- dat[[i]][, c(common_cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
        cor_MAP <- cor_MAP[, c(common_cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                               "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat")]
        if ("1S-FIML" %in% unique(cor_MAP$analType)) {
          cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_MAP$estType <- "MAP" # add estType
        colnames(cor_MAP) <- newcolNames_cor
        
        cor_ogest <- dat[[i]][, c(common_cols, "pop.cor", "ogcor", "ogcor.SE", 
                                  "ogcor.low", "ogcor.up")]
        cor_ogest$prior1 <- NA; cor_ogest$prior2 <- NA 
        cor_ogest$cor.MCMC_neff <- NA; cor_ogest$cor.MCMC_Rhat <- NA # create dummy prior columns
        cor_ogest <- cor_ogest[cor_ogest$analType == "1S-FIML", ] # remove redundant rows
        cor_ogest$estType <- "ogFIML" # add estType
        colnames(cor_ogest) <- newcolNames_cor
        
        cor_result <- rbind(cor_EAP, cor_MAP, cor_ogest)
        
      } else if (i == "SD_result") {
        sd_EAP <- dat[[i]][, c(common_cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                               "prior1", "prior2", "Esd.n_eff", "Esd.Rhat")]
        if ("1S-FIML" %in% unique(sd_EAP$analType)) { 
          sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_EAP$estType <- "EAP" # add estType
        colnames(sd_EAP) <- newcolNames_SD
        
        sd_MAP <- dat[[i]][, c(common_cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", 
                               "prior2", "Esd.n_eff", "Esd.Rhat")]
        sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
        sd_MAP <- sd_MAP[, c(common_cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                             "prior1", "prior2", "Esd.n_eff", "Esd.Rhat")]
        if ("1S-FIML" %in% unique(sd_MAP$analType)) {
          sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_MAP$estType <- "MAP" # add estType
        colnames(sd_MAP) <- newcolNames_SD
        
        sd_ogest <- dat[[i]][, c(common_cols, "pop.SD", "ogsd", "ogsd.SE", 
                                 "ogsd.low", "ogsd.up")]
        sd_ogest$prior1 <- NA; sd_ogest$prior2 <- NA 
        sd_ogest$SD.MCMC_neff <- NA; sd_ogest$SD.MCMC_Rhat <- NA # create dummy prior columns
        sd_ogest <- sd_ogest[sd_ogest$analType == "1S-FIML", ] # remove redundant rows
        sd_ogest$estType <- "ogFIML" # add estType
        colnames(sd_ogest) <- newcolNames_SD
        
        SD_result <- rbind(sd_EAP, sd_MAP, sd_ogest)
        
      }
    }
    all_list <- list(cov_result = cov_result, cor_result = cor_result, 
                     SD_result = SD_result, mPSRF_result = dat$mPSRF_result)
  }
  
  for (i in resType) {
    if (i == "cov_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f1@A", "f2@A~~f2@A", "f3@A~~f3@A", # ego var
                                                   "f1@P~~f1@P", "f2@P~~f2@P", "f3@P~~f3@P", # alter var
                                                   "f1@A~~f2@A", "f1@A~~f3@A", "f2@A~~f3@A", # e-e cov
                                                   "f1@P~~f2@P", "f1@P~~f3@P", "f2@P~~f3@P", # a-a cov
                                                   "f1@A~~f1@P", "f1@A~~f2@P", "f1@A~~f3@P", # e-a cov
                                                   "f2@A~~f2@P", "f2@A~~f3@P", "f3@A~~f3@P",
                                                   "f1@P~~f2@A", "f1@P~~f3@A", "f2@P~~f3@A", # a-e cov
                                                   "f1@AP~~f1@AP", "f2@AP~~f2@AP", "f3@AP~~f3@AP", # rel var
                                                   "f1@AP~~f1@PA", "f2@AP~~f2@PA", "f3@AP~~f3@PA", # dyadic cov
                                                   "f1@AP~~f2@AP", "f1@AP~~f3@AP", "f2@AP~~f3@AP", # intra cov
                                                   "f1@AP~~f2@PA", "f1@AP~~f3@PA", "f2@AP~~f3@PA" # inter cov
                                        ))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-10", "8-10", "10-10", "20-10", 
                                                   "6-25", "8-25", "10-25", "20-25"))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
    } else if (i == "cor_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f2@A", "f1@A~~f3@A", "f2@A~~f3@A", # e-e cov
                                                   "f1@P~~f2@P", "f1@P~~f3@P", "f2@P~~f3@P", # a-a cov
                                                   "f1@A~~f1@P", "f1@A~~f2@P", "f1@A~~f3@P", # e-a cov
                                                   "f2@A~~f2@P", "f2@A~~f3@P", "f3@A~~f3@P",
                                                   "f1@P~~f2@A", "f1@P~~f3@A", "f2@P~~f3@A", # a-e cov
                                                   "f1@AP~~f1@PA", "f2@AP~~f2@PA", "f3@AP~~f3@PA", # dyadic cov
                                                   "f1@AP~~f2@AP", "f1@AP~~f3@AP", "f2@AP~~f3@AP", # intra cov
                                                   "f1@AP~~f2@PA", "f1@AP~~f3@PA", "f2@AP~~f3@PA" # inter cov
                                        ))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-10", "8-10", "10-10", "20-10", 
                                                   "6-25", "8-25", "10-25", "20-25"))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
    } else if (i == "SD_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f1@A", "f2@A~~f2@A", "f3@A~~f3@A", # ego var
                                                   "f1@P~~f1@P", "f2@P~~f2@P", "f3@P~~f3@P", # alter var
                                                   "f1@AP~~f1@AP", "f2@AP~~f2@AP", "f3@AP~~f3@AP" # rel var
                                        ))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-10", "8-10", "10-10", "20-10", 
                                                   "6-25", "8-25", "10-25", "20-25"))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
    } else if (i == "mPSRF_result") {
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-10", "8-10", "10-10", "20-10", 
                                                   "6-25", "8-25", "10-25", "20-25"))
    }
  }
  
  all_list
}

#### READ ME: made a separate function for sim3

# estType = c("EAP", "MAP", "allMCMC", "FIML1S", "all")

# make_long_single(dat = foo, estType = "EAP") -> long1
# make_long_single(dat = foo, estType = "MAP") -> long2
# make_long_single(dat = foo, estType = "allMCMC") -> long3
# make_long_single(dat = foo, estType = "FIML1S") -> long4
# make_long_single(dat = foo, estType = "all") -> long5

#----

# function 3: convert results into long format (multiple samples at once)----

make_long <- function(...) {
  wide_list <- list(...)
  estTypes <- names(wide_list)
  long_list <- lapply(1:length(estTypes), function (i) make_long_single(wide_list[[i]], 
                                                                        estType = estTypes[i]))
  cov_long_list <- lapply(1:length(long_list), function(ll) long_list[[ll]]$cov_result)
  cov_long <- do.call("rbind", cov_long_list)
  cor_long_list <- lapply(1:length(long_list), function(mm) long_list[[mm]]$cor_result)
  cor_long <- do.call("rbind", cor_long_list)
  SD_long_list <- lapply(1:length(long_list), function(nn) long_list[[nn]]$SD_result)
  SD_long <- do.call("rbind", SD_long_list)
  mPSRF_list <- lapply(1:length(long_list), function(pp) long_list[[pp]]$mPSRF_result)
  mPSRF_long <- do.call("rbind", mPSRF_list)
  
  all_list <- list(cov_result = cov_long, cor_result = cor_long, 
                                 SD_result = SD_long, mPSRF_result = mPSRF_long)
  all_list
  }

# tryme <- make_long(allMCMC = extract_results(analType = "default", sim = "sim1"), 
#           allMCMC = extract_results(analType = "prophetic", precision = 0.1, 
#                                  sim = "sim1"), 
#           allMCMC = extract_results(sim = "sim2"), 
#           FIML1S = extract_results(analType = "FIML1S", sim = "sim1"))

#----

# function 4: convert results into long format for simulation 3 (extra $mPSRF column)----
make_long_sim3 <- function(dat, estType) {
  resType <- c("cov_result", "cor_result", "SD_result")
  common_cols <- c("par_names", "MCSampID", "n", "G", "condition", "level", "analType", "iter",
                   "RunTime")
  newcolNames_cov <- c(common_cols, "pop.cov", "cov.est", "cov.SE", "cov.low", "cov.up", 
                       "cov.MCMC_neff", "cov.MCMC_Rhat", "mPSRF", "estType")
  newcolNames_cor <- c(common_cols, "pop.cor", "cor.est", "cor.SE", "cor.low", "cor.up", 
                       "prior1", "prior2", "cor.MCMC_neff", "cor.MCMC_Rhat", "mPSRF", "estType")
  newcolNames_SD <- c(common_cols, "pop.SD", "SD.est", "SD.SE", "SD.low", "SD.up", 
                      "prior1", "prior2", "SD.MCMC_neff", "SD.MCMC_Rhat", "mPSRF", "estType")
  
  all_list <- list()
  
  if (estType == "EAP") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_EAP <- dat[[i]][, c(common_cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", 
                                "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cov_EAP$analType)) {
          cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_EAP$estType <- "EAP" # add estType
        colnames(cov_EAP) <- newcolNames_cov
        
      } else if (i == "cor_result") {
        cor_EAP <- dat[[i]][, c(common_cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low",
                                "Ecor.up", "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cor_EAP$analType)) {
          cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_EAP$estType <- "EAP" # add estType
        colnames(cor_EAP) <- newcolNames_cor
        
      } else if (i == "SD_result") {
        sd_EAP <- dat[[i]][, c(common_cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                               "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(sd_EAP$analType)) { 
          sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_EAP$estType <- "EAP" # add estType
        colnames(sd_EAP) <- newcolNames_SD
        
      }
    }
    all_list <- list(cov_result = cov_EAP, cor_result = cor_EAP, 
                     SD_result = sd_EAP, mPSRF_result = dat$mPSRF_result)
    
  } else if (estType == "MAP") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_MAP <- dat[[i]][, c(common_cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up", 
                                "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column, then rearrange columns
        cov_MAP <- cov_MAP[, c(common_cols, "pop.cov", "Mcov", "Mcov.SE", 
                               "Mcov.low", "Mcov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cov_MAP$analType)) {
          cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_MAP$estType <- "MAP" # add estType
        colnames(cov_MAP) <- newcolNames_cov
        
      } else if (i == "cor_result") {
        cor_MAP <- dat[[i]][, c(common_cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
        cor_MAP <- cor_MAP[, c(common_cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                               "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cor_MAP$analType)) {
          cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_MAP$estType <- "MAP" # add estType
        colnames(cor_MAP) <- newcolNames_cor
        
      } else if (i == "SD_result") {
        sd_MAP <- dat[[i]][, c(common_cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", 
                               "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
        sd_MAP <- sd_MAP[, c(common_cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                             "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(sd_MAP$analType)) {
          sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_MAP$estType <- "MAP" # add estType
        colnames(sd_MAP) <- newcolNames_SD
      } 
    }
    all_list <- list(cov_result = cov_MAP, cor_result = cor_MAP, 
                     SD_result = sd_MAP, mPSRF_result = dat$mPSRF_result)
    
  } else if (estType == "allMCMC") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_EAP <- dat[[i]][, c(common_cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", 
                                "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cov_EAP$analType)) {
          cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_EAP$estType <- "EAP" # add estType
        colnames(cov_EAP) <- newcolNames_cov
        
        cov_MAP <- dat[[i]][, c(common_cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up", 
                                "Ecov.n_eff", "Ecov.Rhat")]
        cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column, then rearrange columns
        cov_MAP <- cov_MAP[, c(common_cols, "pop.cov", "Mcov", "Mcov.SE", 
                               "Mcov.low", "Mcov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cov_MAP$analType)) {
          cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_MAP$estType <- "MAP" # add estType
        colnames(cov_MAP) <- newcolNames_cov
        
        cov_result <- rbind(cov_EAP, cov_MAP)
        
      } else if (i == "cor_result") {
        cor_EAP <- dat[[i]][, c(common_cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low", "Ecor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cor_EAP$analType)) {
          cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_EAP$estType <- "EAP" # add estType
        colnames(cor_EAP) <- newcolNames_cor
        
        cor_MAP <- dat[[i]][, c(common_cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
        cor_MAP <- cor_MAP[, c(common_cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                               "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cor_MAP$analType)) {
          cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_MAP$estType <- "MAP" # add estType
        colnames(cor_MAP) <- newcolNames_cor
        
        cor_result <- rbind(cor_EAP, cor_MAP)
        
      } else if (i == "SD_result") {
        sd_EAP <- dat[[i]][, c(common_cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                               "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(sd_EAP$analType)) { 
          sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_EAP$estType <- "EAP" # add estType
        colnames(sd_EAP) <- newcolNames_SD
        
        sd_MAP <- dat[[i]][, c(common_cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", 
                               "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
        sd_MAP <- sd_MAP[, c(common_cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                             "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(sd_MAP$analType)) {
          sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_MAP$estType <- "MAP" # add estType
        colnames(sd_MAP) <- newcolNames_SD
        
        SD_result <- rbind(sd_EAP, sd_MAP)
      }
    }
    all_list <- list(cov_result = cov_result, cor_result = cor_result, 
                     SD_result = SD_result, mPSRF_result = dat$mPSRF_result)
    
  } else if (estType == "FIML1S") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_ogest <- dat[[i]][, c(common_cols, "pop.cov", "ogcov", "ogcov.SE", 
                                  "ogcov.low", "ogcov.up")]
        cov_ogest <- cov_ogest[cov_ogest$analType == "1S-FIML", ] # remove redundant rows
        cov_ogest$cov.MCMC_neff <- NA; cov_ogest$cov.MCMC_Rhat <- NA; cov_ogest$mPSRF <- NA
        cov_ogest$estType <- "ogFIML" # add estType
        colnames(cov_ogest) <- newcolNames_cov
        
      } else if (i == "cor_result") {
        cor_ogest <- dat[[i]][, c(common_cols, "pop.cor", "ogcor", "ogcor.SE", 
                                  "ogcor.low", "ogcor.up")]
        cor_ogest$prior1 <- NA; cor_ogest$prior2 <- NA 
        cor_ogest$cor.MCMC_neff <- NA; cor_ogest$cor.MCMC_Rhat <- NA; cor_ogest$mPSRF <- NA # create dummy prior columns
        cor_ogest <- cor_ogest[cor_ogest$analType == "1S-FIML", ] # remove redundant rows
        cor_ogest$estType <- "ogFIML" # add estType
        colnames(cor_ogest) <- newcolNames_cor
        
      } else if (i == "SD_result") {
        sd_ogest <- dat[[i]][, c(common_cols, "pop.SD", "ogsd", "ogsd.SE", 
                                 "ogsd.low", "ogsd.up")]
        sd_ogest$prior1 <- NA; sd_ogest$prior2 <- NA 
        sd_ogest$SD.MCMC_neff <- NA; sd_ogest$SD.MCMC_Rhat <- NA; sd_ogest$mPSRF <- NA # create dummy prior columns
        sd_ogest <- sd_ogest[sd_ogest$analType == "1S-FIML", ] # remove redundant rows
        sd_ogest$estType <- "ogFIML" # add estType
        colnames(sd_ogest) <- newcolNames_SD
      }
    }
    all_list <- list(cov_result = cov_ogest, cor_result = cor_ogest, 
                     SD_result = sd_ogest, mPSRF_result = NULL)
    
  } else if (estType == "all") {
    for (i in resType) {
      if (i == "cov_result") {
        cov_EAP <- dat[[i]][, c(common_cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", 
                                "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cov_EAP$analType)) {
          cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_EAP$estType <- "EAP" # add estType
        colnames(cov_EAP) <- newcolNames_cov
        
        cov_MAP <- dat[[i]][, c(common_cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up", 
                                "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column, then rearrange columns
        cov_MAP <- cov_MAP[, c(common_cols, "pop.cov", "Mcov", "Mcov.SE", 
                               "Mcov.low", "Mcov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cov_MAP$analType)) {
          cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cov_MAP$estType <- "MAP" # add estType
        colnames(cov_MAP) <- newcolNames_cov
        
        cov_ogest <- dat[[i]][, c(common_cols, "pop.cov", "ogcov", "ogcov.SE", 
                                  "ogcov.low", "ogcov.up")]
        cov_ogest <- cov_ogest[cov_ogest$analType == "1S-FIML", ] # remove redundant rows
        cov_ogest$cov.MCMC_neff <- NA; cov_ogest$cov.MCMC_Rhat <- NA; cov_ogest$mPSRF <- NA
        cov_ogest$estType <- "ogFIML" # add estType
        colnames(cov_ogest) <- newcolNames_cov
        
        cov_result <- rbind(cov_EAP, cov_MAP, cov_ogest)
        
      } else if (i == "cor_result") {
        cor_EAP <- dat[[i]][, c(common_cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low", "Ecor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cor_EAP$analType)) {
          cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_EAP$estType <- "EAP" # add estType
        colnames(cor_EAP) <- newcolNames_cor
        
        cor_MAP <- dat[[i]][, c(common_cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", 
                                "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
        cor_MAP <- cor_MAP[, c(common_cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                               "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(cor_MAP$analType)) {
          cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        cor_MAP$estType <- "MAP" # add estType
        colnames(cor_MAP) <- newcolNames_cor
        
        cor_ogest <- dat[[i]][, c(common_cols, "pop.cor", "ogcor", "ogcor.SE", 
                                  "ogcor.low", "ogcor.up")]
        cor_ogest$prior1 <- NA; cor_ogest$prior2 <- NA 
        cor_ogest$cor.MCMC_neff <- NA; cor_ogest$cor.MCMC_Rhat <- NA; cor_ogest$mPSRF <- NA # create dummy prior columns
        cor_ogest <- cor_ogest[cor_ogest$analType == "1S-FIML", ] # remove redundant rows
        cor_ogest$estType <- "ogFIML" # add estType
        colnames(cor_ogest) <- newcolNames_cor
        
        cor_result <- rbind(cor_EAP, cor_MAP, cor_ogest)
        
      } else if (i == "SD_result") {
        sd_EAP <- dat[[i]][, c(common_cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                               "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(sd_EAP$analType)) { 
          sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_EAP$estType <- "EAP" # add estType
        colnames(sd_EAP) <- newcolNames_SD
        
        sd_MAP <- dat[[i]][, c(common_cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", 
                               "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
        sd_MAP <- sd_MAP[, c(common_cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                             "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
        if ("1S-FIML" %in% unique(sd_MAP$analType)) {
          sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML", ] # remove redundant rows
        }
        sd_MAP$estType <- "MAP" # add estType
        colnames(sd_MAP) <- newcolNames_SD
        
        sd_ogest <- dat[[i]][, c(common_cols, "pop.SD", "ogsd", "ogsd.SE", 
                                 "ogsd.low", "ogsd.up")]
        sd_ogest$prior1 <- NA; sd_ogest$prior2 <- NA 
        sd_ogest$SD.MCMC_neff <- NA; sd_ogest$SD.MCMC_Rhat <- NA; sd_ogest$mPSRF <- NA # create dummy prior columns
        sd_ogest <- sd_ogest[sd_ogest$analType == "1S-FIML", ] # remove redundant rows
        sd_ogest$estType <- "ogFIML" # add estType
        colnames(sd_ogest) <- newcolNames_SD
        
        SD_result <- rbind(sd_EAP, sd_MAP, sd_ogest)
        
      }
    }
    all_list <- list(cov_result = cov_result, cor_result = cor_result, 
                     SD_result = SD_result, mPSRF_result = dat$mPSRF_result)
  }
  
  for (i in resType) {
    if (i == "cov_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f1@A", "f2@A~~f2@A", "f3@A~~f3@A", # ego var
                                                   "f1@P~~f1@P", "f2@P~~f2@P", "f3@P~~f3@P", # alter var
                                                   "f1@A~~f2@A", "f1@A~~f3@A", "f2@A~~f3@A", # e-e cov
                                                   "f1@P~~f2@P", "f1@P~~f3@P", "f2@P~~f3@P", # a-a cov
                                                   "f1@A~~f1@P", "f1@A~~f2@P", "f1@A~~f3@P", # e-a cov
                                                   "f2@A~~f2@P", "f2@A~~f3@P", "f3@A~~f3@P",
                                                   "f1@P~~f2@A", "f1@P~~f3@A", "f2@P~~f3@A", # a-e cov
                                                   "f1@AP~~f1@AP", "f2@AP~~f2@AP", "f3@AP~~f3@AP", # rel var
                                                   "f1@AP~~f1@PA", "f2@AP~~f2@PA", "f3@AP~~f3@PA", # dyadic cov
                                                   "f1@AP~~f2@AP", "f1@AP~~f3@AP", "f2@AP~~f3@AP", # intra cov
                                                   "f1@AP~~f2@PA", "f1@AP~~f3@PA", "f2@AP~~f3@PA" # inter cov
                                        ))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-50", "8-50", "6-100"))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
    } else if (i == "cor_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f2@A", "f1@A~~f3@A", "f2@A~~f3@A", # e-e cov
                                                   "f1@P~~f2@P", "f1@P~~f3@P", "f2@P~~f3@P", # a-a cov
                                                   "f1@A~~f1@P", "f1@A~~f2@P", "f1@A~~f3@P", # e-a cov
                                                   "f2@A~~f2@P", "f2@A~~f3@P", "f3@A~~f3@P",
                                                   "f1@P~~f2@A", "f1@P~~f3@A", "f2@P~~f3@A", # a-e cov
                                                   "f1@AP~~f1@PA", "f2@AP~~f2@PA", "f3@AP~~f3@PA", # dyadic cov
                                                   "f1@AP~~f2@AP", "f1@AP~~f3@AP", "f2@AP~~f3@AP", # intra cov
                                                   "f1@AP~~f2@PA", "f1@AP~~f3@PA", "f2@AP~~f3@PA" # inter cov
                                        ))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-50", "8-50", "6-100"))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
    } else if (i == "SD_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f1@A", "f2@A~~f2@A", "f3@A~~f3@A", # ego var
                                                   "f1@P~~f1@P", "f2@P~~f2@P", "f3@P~~f3@P", # alter var
                                                   "f1@AP~~f1@AP", "f2@AP~~f2@AP", "f3@AP~~f3@AP" # rel var
                                        ))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-50", "8-50", "6-100"))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
    } else if (i == "mPSRF_result") {
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-50", "8-50", "6-100"))
    }
  }
  
  all_list
}

#-----

# function 5: convert results into long format for simulation 4 (analType names are different than sim3)----
make_long_sim4 <- function(dat) { 
  resType <- c("cov_result", "cor_result", "SD_result")
  common_cols <- c("par_names", "MCSampID", "n", "G", "condition", "level", "analType", "iter",
                   "RunTime")
  newcolNames_cov <- c(common_cols, "pop.cov", "cov.est", "cov.SE", "cov.low", "cov.up", 
                       "cov.MCMC_neff", "cov.MCMC_Rhat", "mPSRF", "estType")
  newcolNames_cor <- c(common_cols, "pop.cor", "cor.est", "cor.SE", "cor.low", "cor.up", 
                       "prior1", "prior2", "cor.MCMC_neff", "cor.MCMC_Rhat", "mPSRF", "estType")
  newcolNames_SD <- c(common_cols, "pop.SD", "SD.est", "SD.SE", "SD.low", "SD.up", 
                      "prior1", "prior2", "SD.MCMC_neff", "SD.MCMC_Rhat", "mPSRF", "estType")
  
  all_list <- list()
  
  for (i in resType) {
    if (i == "cov_result") {
      cov_EAP <- dat[[i]][, c(common_cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", 
                              "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
      if ("1S-FIML-smallvar" %in% unique(cov_EAP$analType)) {
        cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML-smallvar", ] # remove redundant rows
      }
      cov_EAP$estType <- "EAP" # add estType
      colnames(cov_EAP) <- newcolNames_cov
      
      cov_MAP <- dat[[i]][, c(common_cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up", 
                              "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
      cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column, then rearrange columns
      cov_MAP <- cov_MAP[, c(common_cols, "pop.cov", "Mcov", "Mcov.SE", 
                             "Mcov.low", "Mcov.up", "Ecov.n_eff", "Ecov.Rhat", "mPSRF")]
      if ("1S-FIML-smallvar" %in% unique(cov_MAP$analType)) {
        cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML-smallvar", ] # remove redundant rows
      }
      cov_MAP$estType <- "MAP" # add estType
      colnames(cov_MAP) <- newcolNames_cov
      
      cov_ogest <- dat[[i]][, c(common_cols, "pop.cov", "ogcov", "ogcov.SE", 
                                "ogcov.low", "ogcov.up")]
      cov_ogest <- cov_ogest[cov_ogest$analType == "1S-FIML-smallvar", ] # remove redundant rows
      cov_ogest$cov.MCMC_neff <- NA; cov_ogest$cov.MCMC_Rhat <- NA; cov_ogest$mPSRF <- NA
      cov_ogest$estType <- "ogFIML" # add estType
      colnames(cov_ogest) <- newcolNames_cov
      
      cov_result <- rbind(cov_EAP, cov_MAP, cov_ogest)
      
    } else if (i == "cor_result") {
      cor_EAP <- dat[[i]][, c(common_cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low", "Ecor.up", 
                              "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
      if ("1S-FIML-smallvar" %in% unique(cor_EAP$analType)) {
        cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML-smallvar", ] # remove redundant rows
      }
      cor_EAP$estType <- "EAP" # add estType
      colnames(cor_EAP) <- newcolNames_cor
      
      cor_MAP <- dat[[i]][, c(common_cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", 
                              "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
      cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
      cor_MAP <- cor_MAP[, c(common_cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                             "prior1", "prior2", "Ecor.n_eff", "Ecor.Rhat", "mPSRF")]
      if ("1S-FIML-smallvar" %in% unique(cor_MAP$analType)) {
        cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML-smallvar", ] # remove redundant rows
      }
      cor_MAP$estType <- "MAP" # add estType
      colnames(cor_MAP) <- newcolNames_cor
      
      cor_ogest <- dat[[i]][, c(common_cols, "pop.cor", "ogcor", "ogcor.SE", 
                                "ogcor.low", "ogcor.up")]
      cor_ogest$prior1 <- NA; cor_ogest$prior2 <- NA 
      cor_ogest$cor.MCMC_neff <- NA; cor_ogest$cor.MCMC_Rhat <- NA; cor_ogest$mPSRF <- NA # create dummy prior columns
      cor_ogest <- cor_ogest[cor_ogest$analType == "1S-FIML-smallvar", ] # remove redundant rows
      cor_ogest$estType <- "ogFIML" # add estType
      colnames(cor_ogest) <- newcolNames_cor
      
      cor_result <- rbind(cor_EAP, cor_MAP, cor_ogest)
      
    } else if (i == "SD_result") {
      sd_EAP <- dat[[i]][, c(common_cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                             "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
      if ("1S-FIML-smallvar" %in% unique(sd_EAP$analType)) { 
        sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML-smallvar", ] # remove redundant rows
      }
      sd_EAP$estType <- "EAP" # add estType
      colnames(sd_EAP) <- newcolNames_SD
      
      sd_MAP <- dat[[i]][, c(common_cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", 
                             "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
      sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
      sd_MAP <- sd_MAP[, c(common_cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                           "prior1", "prior2", "Esd.n_eff", "Esd.Rhat", "mPSRF")]
      if ("1S-FIML-smallvar" %in% unique(sd_MAP$analType)) {
        sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML-smallvar", ] # remove redundant rows
      }
      sd_MAP$estType <- "MAP" # add estType
      colnames(sd_MAP) <- newcolNames_SD
      
      sd_ogest <- dat[[i]][, c(common_cols, "pop.SD", "ogsd", "ogsd.SE", 
                               "ogsd.low", "ogsd.up")]
      sd_ogest$prior1 <- NA; sd_ogest$prior2 <- NA 
      sd_ogest$SD.MCMC_neff <- NA; sd_ogest$SD.MCMC_Rhat <- NA; sd_ogest$mPSRF <- NA # create dummy prior columns
      sd_ogest <- sd_ogest[sd_ogest$analType == "1S-FIML-smallvar", ] # remove redundant rows
      sd_ogest$estType <- "ogFIML" # add estType
      colnames(sd_ogest) <- newcolNames_SD
      
      SD_result <- rbind(sd_EAP, sd_MAP, sd_ogest)
      
    }
  }
  all_list <- list(cov_result = cov_result, cor_result = cor_result, 
                   SD_result = SD_result, mPSRF_result = dat$mPSRF_result)
  
  for (i in resType) {
    if (i == "cov_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f1@A", "f2@A~~f2@A", "f3@A~~f3@A", # ego var
                                                   "f1@P~~f1@P", "f2@P~~f2@P", "f3@P~~f3@P", # alter var
                                                   "f1@A~~f2@A", "f1@A~~f3@A", "f2@A~~f3@A", # e-e cov
                                                   "f1@P~~f2@P", "f1@P~~f3@P", "f2@P~~f3@P", # a-a cov
                                                   "f1@A~~f1@P", "f1@A~~f2@P", "f1@A~~f3@P", # e-a cov
                                                   "f2@A~~f2@P", "f2@A~~f3@P", "f3@A~~f3@P",
                                                   "f1@P~~f2@A", "f1@P~~f3@A", "f2@P~~f3@A", # a-e cov
                                                   "f1@AP~~f1@AP", "f2@AP~~f2@AP", "f3@AP~~f3@AP", # rel var
                                                   "f1@AP~~f1@PA", "f2@AP~~f2@PA", "f3@AP~~f3@PA", # dyadic cov
                                                   "f1@AP~~f2@AP", "f1@AP~~f3@AP", "f2@AP~~f3@AP", # intra cov
                                                   "f1@AP~~f2@PA", "f1@AP~~f3@PA", "f2@AP~~f3@PA" # inter cov
                                        ))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                        levels = c("6-10", "10-10"))
      
    } else if (i == "cor_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f2@A", "f1@A~~f3@A", "f2@A~~f3@A", # e-e cov
                                                   "f1@P~~f2@P", "f1@P~~f3@P", "f2@P~~f3@P", # a-a cov
                                                   "f1@A~~f1@P", "f1@A~~f2@P", "f1@A~~f3@P", # e-a cov
                                                   "f2@A~~f2@P", "f2@A~~f3@P", "f3@A~~f3@P",
                                                   "f1@P~~f2@A", "f1@P~~f3@A", "f2@P~~f3@A", # a-e cov
                                                   "f1@AP~~f1@PA", "f2@AP~~f2@PA", "f3@AP~~f3@PA", # dyadic cov
                                                   "f1@AP~~f2@AP", "f1@AP~~f3@AP", "f2@AP~~f3@AP", # intra cov
                                                   "f1@AP~~f2@PA", "f1@AP~~f3@PA", "f2@AP~~f3@PA" # inter cov
                                        ))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                       levels = c("6-10", "10-10"))
      
    } else if (i == "SD_result") {
      all_list[[i]]$par_names <- factor(all_list[[i]]$par_names, 
                                        levels = c("f1@A~~f1@A", "f2@A~~f2@A", "f3@A~~f3@A", # ego var
                                                   "f1@P~~f1@P", "f2@P~~f2@P", "f3@P~~f3@P", # alter var
                                                   "f1@AP~~f1@AP", "f2@AP~~f2@AP", "f3@AP~~f3@AP" # rel var
                                        ))
      all_list[[i]]$level <- factor(all_list[[i]]$level, levels = c("dyad", "case"),
                                    labels = c("Dyad level", "Case level"))
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                       levels = c("6-10", "10-10"))
      
    } else if (i == "mPSRF_result") {
      all_list[[i]]$condition <- factor(all_list[[i]]$condition, 
                                       levels = c("6-10", "10-10"))
    }
  }
  
  all_list
} # extract all estimates
#----

# sim1result <- make_long_single(extract_results(sim = "sim1"), estType = "all")
# saveRDS(sim1result, file = "sim1result.rds")
# 
# sim2result <- make_long(EAP = extract_results(analType = "prophetic", precision = 0.1, sim = "sim1"),
#                         FIML1S = extract_results(analType = "FIML1S", sim = "sim1"),
#                         EAP = extract_results(sim = "sim2"))
# saveRDS(sim2result, file = "sim2result.rds")
# 
# sim3result <- make_long_sim3(dat = extract_results(sim = "sim3"), estType = "EAP")
# saveRDS(sim3result, file = "sim3result.rds")
# 
# sim4result <- make_long_sim4(dat = extract_results(sim = "sim4"))
# saveRDS(sim4result, file = "sim4result.rds")
# 
# BMAtestresult <- make_long_single(extract_results(sim = "BMAtest"), estType = "EAP")
# saveRDS(BMAtestresult, file = "BMAresult.rds")

