## Aditi M. Bhangale
## Last updated: 24 April 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model

# Data processing

# rm(list = ls())

##README MCSampID 198 in FIML1S 6-10 is not in the output. This is probably because
## that sample did not converge

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/results_and_postprocessing")
# getwd()

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

# dat <- extract_results(analType = "default", sim = "sim1")
# dat <- extract_results(analType = "FIML1S", sim = "sim1")

make_long_single <- function(dat) {
  resType <- c("cov_result", "cor_result", "SD_result")
  cols <- c("par_names", "MCSampID", "n", "G", "condition", "level", "analType", "iter",
            "RunTime")
  MCMC_analTypes <- c("MCMC-default-NA", "MCMC-prophetic-0.05", "MCMC-prophetic-0.1",
                      "MCMC-prophetic-0.2", "MCMC-ANOVA-0.1", "MCMC-FIML-0.1",
                      "MCMC-thoughtful-0.1") #TODO start from here on 24 April 2024
  
  all_list <- list()
  
  for (i in resType) {
    if (i == "cov_result") {
      cov_EAP <- dat[[i]][, c(cols, "pop.cov", "Ecov", "Ecov.SE", "Ecov.low", "Ecov.up")]
      if ("1S-FIML" %in% unique(cov_EAP$analType)) cov_EAP <- cov_EAP[!cov_EAP$analType == "1S-FIML", ] # remove redundant rows
      cov_EAP$estType <- "EAP" # add estType
      
      cov_MAP <- dat[[i]][, c(cols, "pop.cov", "Mcov", "Mcov.low", "Mcov.up")]
      cov_MAP$Mcov.SE <- NA # add dummy Mcov.SE column
      cov_MAP <- cov_MAP[, c(cols, "pop.cov", "Mcov", "Mcov.SE", "Mcov.low", "Mcov.up")]
      if ("1S-FIML" %in% unique(cov_MAP$analType)) cov_MAP <- cov_MAP[!cov_MAP$analType == "1S-FIML", ] # remove redundant rows
      cov_MAP$estType <- "MAP" # add estType
      
      if ("1S-FIML" %in% unique(dat[[i]]$analType)) {
        cov_ogest <- dat[[i]][, c(cols, "pop.cov", "ogcov", "ogcov.SE", "ogcov.low", "ogcov.up")]
        cov_ogest <- cov_ogest[cov_ogest$analType == "1S-FIML", ] # remove redundant rows
        cov_ogest$estType <- "ogFIML" # add estType
      } else {
        cov_ogest <- NULL
      }
      
      # new column names
      newcolNames <- c(cols, "pop.cov", "cov.est", "cov.SE", "cov.low", "cov.up", "estType")
      colnames(cov_EAP) <- newcolNames
      colnames(cov_MAP) <- newcolNames
      if (!is.null(cov_ogest)) colnames(cov_ogest) <- newcolNames
      
      cov_result <- rbind(cov_EAP, cov_MAP, cov_ogest)
      
    } else if (i == "cor_result") {
      cor_EAP <- dat[[i]][, c(cols, "pop.cor", "Ecor", "Ecor.SE", "Ecor.low", "Ecor.up", 
                              "prior1", "prior2")]
      if ("1S-FIML" %in% unique(cor_EAP$analType)) cor_EAP <- cor_EAP[!cor_EAP$analType == "1S-FIML", ] # remove redundant rows
      cor_EAP$estType <- "EAP" # add estType
      
      cor_MAP <- dat[[i]][, c(cols, "pop.cor", "Mcor", "Mcor.low", "Mcor.up", "prior1", "prior2")]
      cor_MAP$Mcor.SE <- NA # add dummy Mcor.SE column
      cor_MAP <- cor_MAP[, c(cols, "pop.cor", "Mcor", "Mcor.SE", "Mcor.low", "Mcor.up", 
                             "prior1", "prior2")]
      if ("1S-FIML" %in% unique(cor_MAP$analType)) cor_MAP <- cor_MAP[!cor_MAP$analType == "1S-FIML", ] # remove redundant rows
      cor_MAP$estType <- "MAP" # add estType
      
      if ("1S-FIML" %in% unique(dat[[i]]$analType)) {
        cor_ogest <- dat[[i]][, c(cols, "pop.cor", "ogcor", "ogcor.SE", "ogcor.low", "ogcor.up")]
        cor_ogest$prior1 <- cor_ogest$prior2 <- NA # create dummy prior columns
        cor_ogest <- cor_ogest[cor_ogest$analType == "1S-FIML", ] # remove redundant rows
        cor_ogest$estType <- "ogFIML" # add estType
      } else {
        cor_ogest <- NULL
      }
      
      # new column names
      newcolNames <- c(cols, "pop.cor", "cor.est", "cor.SE", "cor.low", "cor.up", 
                       "prior1", "prior2", "estType")
      colnames(cor_EAP) <- newcolNames
      colnames(cor_MAP) <- newcolNames
      if (!is.null(cor_ogest)) colnames(cor_ogest) <- newcolNames
      
      cor_result <- rbind(cor_EAP, cor_MAP, cor_ogest)
      
    } else if (i == "SD_result") {
      sd_EAP <- dat[[i]][, c(cols, "pop.SD", "Esd", "Esd.SE", "Esd.low", "Esd.up", 
                             "prior1", "prior2")]
      if ("1S-FIML" %in% unique(sd_EAP$analType)) sd_EAP <- sd_EAP[!sd_EAP$analType == "1S-FIML", ] # remove redundant rows
      sd_EAP$estType <- "EAP" # add estType
      
      sd_MAP <- dat[[i]][, c(cols, "pop.SD", "Msd", "Msd.low", "Msd.up", "prior1", "prior2")]
      sd_MAP$Msd.SE <- NA # add dummy Msd.SE column
      sd_MAP <- sd_MAP[, c(cols, "pop.SD", "Msd", "Msd.SE", "Msd.low", "Msd.up", 
                           "prior1", "prior2")]
      if ("1S-FIML" %in% unique(sd_MAP$analType)) sd_MAP <- sd_MAP[!sd_MAP$analType == "1S-FIML", ] # remove redundant rows
      sd_MAP$estType <- "MAP" # add estType
      
      if ("1S-FIML" %in% unique(dat[[i]]$analType)) {
        sd_ogest <- dat[[i]][, c(cols, "pop.SD", "ogsd", "ogsd.SE", "ogsd.low", "ogsd.up")]
        sd_ogest$prior1 <- sd_ogest$prior2 <- NA # create dummy prior columns
        sd_ogest <- sd_ogest[sd_ogest$analType == "1S-FIML", ] # remove redundant rows
        sd_ogest$estType <- "ogFIML" # add estType
      } else {
        sd_ogest <- NULL
      }
      
      # new column names
      newcolNames <- c(cols, "pop.SD", "SD.est", "SD.SE", "SD.low", "SD.up", 
                       "prior1", "prior2", "estType")
      colnames(sd_EAP) <- newcolNames
      colnames(sd_MAP) <- newcolNames
      if (!is.null(sd_ogest)) colnames(sd_ogest) <- newcolNames
      
      SD_result <- rbind(sd_EAP, sd_MAP, sd_ogest)
    }
  }
  
  all_list <- list(cov_result = cov_result, cor_result = cor_result, 
                   SD_result = SD_result, mPSRF_result = dat$mPSRF_result)
  
  all_list
}

## ugh, this is super tricky because what if something doesn't have EAP and MAP and only has
# FIML1S?

bar <- make_long_single(foo)

make_long <- function(dat1 = NULL, dat2 = NULL, dat3 = NULL, dat4 = NULL) {
  if(!is.missing(dat1)) long1 <- make_long_single(dat1)
  if(!is.missing(dat2)) long2 <- make_long_single(dat2)
  if(!is.missing(dat3)) long3 <- make_long_single(dat3)
  if(!is.missing(dat4)) long4 <- make_long_single(dat4)
  
  cov_result <- rbind(long1$cov_result)
}

