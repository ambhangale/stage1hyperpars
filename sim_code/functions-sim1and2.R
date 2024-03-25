## Aditi M. Bhangale
## Last updated: 25 March 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model
# Simulations 1 and 2

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/sim_code")
# getwd()

# MCSampID = 1; n = 5; G = 3
# rr.data <- genGroups(MCSampID = MCSampID, n = n, G = G)
# rr.vars <- c("V1", "V2", "V3")
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"
# precision < 0.1
# targetCorr <- 0.3

# function 0: generate level-specific (co)variance matrices----

getSigma <- function(return_mats = TRUE){
  
  ## person effects
  
  Ui_names <- c("V1@A", "V2@A", "V3@A", "V1@P", "V2@P", "V3@P") # person-level indicator names
  Uf_names <- c("f1@A", "f1@P") # person-level factor names
  
  # person-level lambda
  LAM_U <- matrix(0, 6, 2)
  LAM_U[1,1] <- 1
  LAM_U[2,1] <- 1.2
  LAM_U[3,1] <- 0.7
  LAM_U[4,2] <- 1
  LAM_U[5,2] <- 0.6
  LAM_U[6,2] <- 0.6
  
  dimnames(LAM_U) <- list(Ui_names, Uf_names)
  
  # person-level factor cov matrix
  PHI_U <- matrix(c(0.40, 0.05, 0.05, 0.20), 2, 2)
  
  dimnames(PHI_U) <- list(Uf_names, Uf_names)
  
  # person-level indicator residual cov matrix
  PSI_U <- diag(6)
  PSI_U[1,1] <- 0.2
  PSI_U[2,2] <- 0.2
  PSI_U[3,3] <- 0.2
  PSI_U[4,4] <- 0.1
  PSI_U[1,4] <- PSI_U[4,1] <- 0.05
  PSI_U[5,5] <- 0.1
  PSI_U[6,6] <- 0.1
  PSI_U[3,6] <- PSI_U[6,3] <- 0.03
  
  dimnames(PSI_U) <- list(Ui_names, Ui_names)
  
  ## structure dyad effects
  
  Di_names <- c("V1@AP", "V1@PA", "V2@AP", "V2@PA", "V3@AP", "V3@PA") # dyad-level indicator names
  Df_names <- c("f1@AP", "f1@PA") # dyad-level factor names
  
  # dyad-level lambda
  LAM_D <- matrix(0, 6, 2)
  LAM_D[1,1] <- LAM_D[2,2] <- 1
  LAM_D[3,1] <- LAM_D[4,2] <- 0.8
  LAM_D[5,1] <- LAM_D[6,2] <- 1.4
  
  dimnames(LAM_D) <- list(Di_names, Df_names)
  
  # dyad-level factor cov matrix
  PHI_D <- matrix(c(0.60, 0.15, 0.15, 0.60), 2, 2)
  
  dimnames(PHI_D) <- list(Df_names, Df_names)
  
  # dyad-level indicator residual cov matrix
  PSI_D <- diag(6)
  PSI_D[1,1] <- PSI_D[2,2] <- 0.3
  PSI_D[3,3] <- PSI_D[4,4] <- 0.5
  PSI_D[3,4] <- PSI_D[4,3] <- 0.1
  PSI_D[5,5] <- PSI_D[6,6] <- 0.4
  PSI_D[5,6] <- PSI_D[6,5] <- -0.2
  
  dimnames(PSI_D) <- list(Di_names, Di_names)
  
  # person-level sigma
  SIGMA_U <- LAM_U %*% PHI_U %*% t(LAM_U) + PSI_U
  # dyad-level sigma
  SIGMA_D <- LAM_D %*% PHI_D %*% t(LAM_D) + PSI_D
  
  ## pop values for satmod
  Ufsat_names <- c("f1@A", "f2@A", "f3@A", "f1@P", "f2@P", "f3@P") # person-level satmod factor names
  Dfsat_names <- c("f1@AP", "f1@PA", "f2@AP", "f2@PA", "f3@AP", "f3@PA") # dyad-level satmod factor names
  
  dimnames(SIGMA_U) <- list(Ufsat_names, Ufsat_names)
  dimnames(SIGMA_D) <- list(Dfsat_names, Dfsat_names)
  
  # person-level corr matrix
  R_U <- cov2cor(SIGMA_U)
  # dyad-level corr matrix
  R_D <- cov2cor(SIGMA_D)
  
  mat_list <- list(LAM_U=LAM_U, PHI_U=PHI_U, PSI_U=PSI_U,
                   LAM_D=LAM_D, PHI_D=PHI_D, PSI_D=PSI_D,
                   SIGMA_U = SIGMA_U, SIGMA_D = SIGMA_D, R_U = R_U, R_D = R_D)
  
  if (return_mats)  return(mat_list)
  
  # creating a dataframe of person-level population values
  pop_U.cov <- as.data.frame(as.table(SIGMA_U))
  colnames(pop_U.cov) <- c("row", "col", "pop.cov")
  pop_U.cov$par_names <- paste0(pop_U.cov$row, "~~", pop_U.cov$col)
  pop_U.cov <- pop_U.cov[, c("par_names", "pop.cov")]
  
  # person-level SDs
  mat_U.SD <- diag(sqrt(diag(SIGMA_U)))
  dimnames(mat_U.SD) <- list(Ufsat_names, Ufsat_names)
  pop_U.SD <- as.data.frame(as.table(mat_U.SD))
  pop_U.SD <- pop_U.SD[pop_U.SD$Freq != 0, ]
  colnames(pop_U.SD) <- c("row", "col", "pop.SD")
  pop_U.SD$par_names <- paste0(pop_U.SD$row, "~~", pop_U.SD$col)
  pop_U.SD <- pop_U.SD[, c("par_names", "pop.SD")]
  rownames(pop_U.SD) <- NULL
  
  # person-level correlations
  mat_U.cor <- solve(mat_U.SD) %*% SIGMA_U %*% solve(mat_U.SD)
  pop_U.cor <- as.data.frame(as.table(mat_U.cor))
  colnames(pop_U.cor) <- c("row", "col", "pop.cor")
  pop_U.cor$par_names <- paste0(pop_U.cor$row, "~~", pop_U.cor$col)
  pop_U.cor <- pop_U.cor[, c("par_names", "pop.cor")]
  
  # creating a dataframe of dyad-level population values
  pop_D.cov <- as.data.frame(as.table(SIGMA_D))
  colnames(pop_D.cov) <- c("row", "col", "pop.cov")
  pop_D.cov$par_names <-  paste0(pop_D.cov$row, "~~", pop_D.cov$col)
  pop_D.cov <- pop_D.cov[, c("par_names", "pop.cov")]
  
  # dyad-level SDs
  mat_D.SD <- diag(sqrt(diag(SIGMA_D)))
  dimnames(mat_D.SD) <- list(Dfsat_names, Dfsat_names)
  pop_D.SD <- as.data.frame(as.table(mat_D.SD))
  pop_D.SD <- pop_D.SD[pop_D.SD$Freq != 0, ]
  pop_D.SD <- pop_D.SD[!duplicated(pop_D.SD$Freq), ]
  colnames(pop_D.SD) <- c("row", "col", "pop.SD")
  pop_D.SD$par_names <- paste0(pop_D.SD$row, "~~", pop_D.SD$col)
  pop_D.SD <- pop_D.SD[, c("par_names", "pop.SD")]
  rownames(pop_D.SD) <- NULL
  
  # dyad-level correlations
  mat_D.cor <- solve(mat_D.SD) %*% SIGMA_D %*% solve(mat_D.SD)
  pop_D.cor <- as.data.frame(as.table(mat_D.cor))
  colnames(pop_D.cor) <- c("row", "col", "pop.cor")
  pop_D.cor$par_names <- paste0(pop_D.cor$row, "~~", pop_D.cor$col)
  pop_D.cor <- pop_D.cor[, c("par_names", "pop.cor")]
  
  return(list(pop.cov = rbind(pop_U.cov, pop_D.cov),
              pop.cor = rbind(pop_U.cor, pop_D.cor),
              pop.SD = rbind(pop_U.SD, pop_D.SD)))
}

# getSigma(return_mats = FALSE)
# getSigma(return_mats = TRUE)

#----

# function 1: simulate data for a single group----

genData <- function(n) {
  
  P_Sigma <- getSigma()$SIGMA_U # save person- and dyad-level pop.cov matrices
  D_Sigma <- getSigma()$SIGMA_D
  
  Ui_names <- c("V1@A", "V2@A", "V3@A", "V1@P", "V2@P", "V3@P") # person-level indicator names
  Di_names <- c("V1@AP", "V1@PA", "V2@AP", "V2@PA", "V3@AP", "V3@PA") # dyad-level indicator names
  
  # changing names to generate data per SRM component
  dimnames(P_Sigma) <- list(Ui_names, Ui_names)
  dimnames(D_Sigma) <- list(Di_names, Di_names)
  
  ##MEAN VECTOR-----
  mu <- c(rep(0, 6)) # group-mean centered SRM component variables
  
  ##DATA GENERATION-----
  
  library(rockchalk) # for mvrnorm()
  
  p_dat <- mvrnorm(n = n, mu = mu, Sigma = P_Sigma)
  d_dat <- mvrnorm(n = (n*(n - 1))/2, mu = mu, Sigma = D_Sigma) 
  # Ndyads = n * (n - 1) -- and we already have separate AP and PA columns, 
  # so we need only (n*(n - 1))/2 rows in total
  
  Ydat <- NULL
  Vnames <- c("V1", "V2", "V3")
  
  for(v in Vnames) {
    
    # actor effects
    Aname <- paste0(v, "@A")
    Amat <- matrix(p_dat[, Aname], nrow = n, ncol = n, byrow = FALSE)
    diag(Amat) <- NA
    
    # partner effects
    Pname <- paste0(v, "@P")
    Pmat <- matrix(p_dat[, Pname], nrow = n, ncol = n, byrow = TRUE)
    diag(Pmat) <- NA
    
    # dyad-level effects
    Rmat <- matrix(NA, n, n)
    
    # _ij & _ji effects
    APname <- paste0(v, "@AP")
    PAname <- paste0(v, "@PA")
    
    # build the dyad-level matrix
    foo <- which(lower.tri(Rmat, diag = FALSE), arr.ind = TRUE)
    
    Rmat[foo] <- d_dat[, APname]
    Rmat[foo[, 2:1]] <- d_dat[, PAname]
    
    Y_adj <- Amat + Pmat + Rmat # adjacency matrix for Y
    
    # converting the data to long format
    
    tempY_low <- cbind(foo, Y_adj[foo])
    colnames(tempY_low) <- c("Actor", "Partner", v)
    
    tempY_up <- cbind(foo[, 2:1], Y_adj[foo[, 2:1]])
    colnames(tempY_up) <- c("Actor", "Partner", v)
    
    tempY <- rbind(tempY_low, tempY_up)
    
    # merge with previous variable's long-Y
    if (is.null(Ydat)) {
      Ydat <- tempY
    } else {
      Ydat <- merge(Ydat, tempY) # necessary to generate multigroup data (next function)
    }
  } # end for loop for data generation
  return(Ydat)
  
}

# genData(n = 5)

#----

# function 2: simulate data for multiple groups----

genGroups <- function(MCSampID, n, G) {
  
  # set the seed
  library(parallel)
  library(portableParallelSeeds)
  mySeeds <- seedCreator(nReps = 5000, streamsPerRep = 1, seed = 20505)
  
  setSeeds(projSeeds = mySeeds, run = MCSampID)
  
  gDat <- lapply(1:G, function(g){
    
    gen <- genData(n = n)
    
    gen$Actor <- gen$Actor + g*100
    gen$Partner <- gen$Partner + g*100
    
    cbind(Group = g, gen)
  })
  
  gDat <- do.call("rbind", gDat)
  
  return(gDat)
}

# genGroups(MCSampID = 1, n = 5, G = 3)

#----

# function 3: minimum loss function to optimise over to choose hyperparameters for beta priors----

minLoss <- function(par = c(1.5, 1.5), targetCorr, accept.dev) {
  # targetCorr = population correlation value we want to estimate
  # accept.dev = average deviation from targetCorr we are willing to accept
  
  # default mean (location) of beta distribution (by default, is 0.5)
  defaultM <- par[1] / (par[1] + par[2])
  # default SD of beta distribution (by default, is 0.25)
  defaultSD <- sqrt(par[1]*par[2] / ((par[1]+par[2])^2 * (par[1] + par[2] + 1))) 
  
  
  # target correlation (location of beta) and accepted deviation around it (accept.dev)
  # convert to beta scale
  targetM <- (targetCorr + 1) / 2 # convert targetCorr to beta scale (location of beta)
  targetSD <- accept.dev / 2 # convert accept.dev to beta scale (SD of beta)
  
  mSE <- (defaultM - targetM)^2 # squared error for mean
  sdSE <- (defaultSD - targetSD)^2 # squared error for sd
  
  mSE + sdSE # aim: minimise the sum of squared error (loss function)-- thus, use optim()
}

# optim(par = c(1.5, 1.5), fn = minLoss, targetCorr = 0.3, accept.dev = 0.1,
#       method = "L-BFGS-B", lower = 0)

#----

# function 4: visualise beta (not used in simulation)----

visBeta <- function(a, b, var1, var2, ...) {
  ## variable names available?
  if (!missing(var1) && !missing(var2)) {
    corName <- paste0("Cor(", var1, ", ", var2, ")")
  } else corName <- "Correlation"
  
  ## calculate descriptive summary stats for title
  betaM  <- a/(a+b)
  betaSD <- sqrt((a*b)/((a+b)^2*(a+b+1)))
  corM  <- round(betaM*2 - 1, 2)
  corSD <- round(betaSD*2, 2)
  
  ## now plot
  curve(dbeta((x+1)/2, shape1 = a, shape2 = b), from = -1, to = 1,
        xlab = corName, ylab = "Prior Density",
        main = paste0("M = ", corM, ", SD = ", corSD), ...)
}

#----

# function 5: thoughtful priors----

thoughtful_priors <- function(data, rr.vars, targetCorr, precision) {
  # SDs
  ## m for SDs is already thoughtful/approximate location because lavaan.srm estimates
  ## it based on the data
  priors$rr_in_t$sd <- priors$rr_out_t$sd <- 
    priors$rr_rel_t$sd <- rep(precision, times = 3) # set SD of t-priors to 0.1
  
  # correlations --- same thoughtful targetCorr for all correlations at both case and dyad level
  corr_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss, targetCorr = targetCorr, 
                          accept.dev = precision, method = "L-BFGS-B", lower = 0)
  
  alpha <- corr_hyperpars$par[1]
  beta <- corr_hyperpars$par[2] # save hyperparameters
  
  ## assign case-level hyperpars
  priors$case_beta[lower.tri(priors$case_beta, diag = FALSE)] <- alpha
  priors$case_beta[upper.tri(priors$case_beta, diag = FALSE)] <- beta
  
  ##  assign dyad-level hyperpars
  diag(priors$rr_beta_a) <- priors$rr_beta_a[upper.tri(priors$rr_beta_a)] <- 
    priors$rr_beta_a[lower.tri(priors$rr_beta_a)] <- alpha
  diag(priors$rr_beta_b) <- priors$rr_beta_b[upper.tri(priors$rr_beta_b)] <- 
    priors$rr_beta_b[lower.tri(priors$rr_beta_b)] <- beta
  
  return(priors)
}

# thoughtful_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), 
#                   targetCorr = 0.3, precision = 0.1)

#----


# function xxxx: set customised priors for MCMC stage----

# library(lavaan.srm)

set_priors <- function(data, rr.vars, IDout, IDin, IDgroup, precision, priorType, 
                       multiMLE = FALSE) {
 rr.data <- data
 priors <- srm_priors(rr.data[rr.vars]) # default MCMC priors (diffuse priors)
 
 if (priorType == "default") { # default (diffuse) priors
   priors
 } else if (priorType == "thoughtful") { # thoughtful priors
   
 } else if (priorType == "prophetic") { # prophetic priors
   
 } else if (priorType == "ANOVA") { # method-of-moments priors (ANOVA-based, `TripleR`)
   
 } else if (priorType = "FIML") { # FIML-based priors (`srm`)
   
 }
 return(priors)
}

#----
