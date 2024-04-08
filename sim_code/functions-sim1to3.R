## Aditi M. Bhangale
## Last updated: 8 April 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model
# Simulations 1 to 3

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/sim_code")
# getwd()

# MCSampID = 1; n = 5; G = 3
# rr.data <- genGroups(MCSampID = 1, n = 5, G = 3)
# rr.vars <- c("V1", "V2", "V3")
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"
# precision <- 0.1
# targetCorr <- 0.3
# priorType = "default"
# iter = 10

#####################################
# FUNCTIONS FOR SIMULATIONS 1-3 ----
#####################################

# function 0: generate level-specific (co)variance matrices----

getSigma <- function(return_mats = TRUE){
  
  ## person effects
  
  Inames_c <- c("V1@A", "V2@A", "V3@A", "V1@P", "V2@P", "V3@P") # person-level indicator names
  Fnames_c <- c("f1@A", "f1@P") # person-level factor names
  
  # person-level lambda
  LAM_c <- matrix(0, 6, 2)
  LAM_c[1,1] <- 1
  LAM_c[2,1] <- 1.2
  LAM_c[3,1] <- 0.7
  LAM_c[4,2] <- 1
  LAM_c[5,2] <- 0.6
  LAM_c[6,2] <- 0.6
  
  dimnames(LAM_c) <- list(Inames_c, Fnames_c)
  
  # person-level factor cov matrix
  PHI_c <- matrix(c(0.40, 0.05, 0.05, 0.20), 2, 2)
  
  dimnames(PHI_c) <- list(Fnames_c, Fnames_c)
  
  # person-level indicator residual cov matrix
  PSI_c <- diag(6)
  PSI_c[1,1] <- 0.2
  PSI_c[2,2] <- 0.2
  PSI_c[3,3] <- 0.2
  PSI_c[4,4] <- 0.1
  PSI_c[1,4] <- PSI_c[4,1] <- 0.05
  PSI_c[5,5] <- 0.1
  PSI_c[6,6] <- 0.1
  PSI_c[3,6] <- PSI_c[6,3] <- 0.03
  
  dimnames(PSI_c) <- list(Inames_c, Inames_c)
  
  ## structure dyad effects
  
  Inames_d <- c("V1@AP", "V1@PA", "V2@AP", "V2@PA", "V3@AP", "V3@PA") # dyad-level indicator names
  Fnames_d <- c("f1@AP", "f1@PA") # dyad-level factor names
  
  # dyad-level lambda
  LAM_d <- matrix(0, 6, 2)
  LAM_d[1,1] <- LAM_d[2,2] <- 1
  LAM_d[3,1] <- LAM_d[4,2] <- 0.8
  LAM_d[5,1] <- LAM_d[6,2] <- 1.4
  
  dimnames(LAM_d) <- list(Inames_d, Fnames_d)
  
  # dyad-level factor cov matrix
  PHI_d <- matrix(c(0.60, 0.15, 0.15, 0.60), 2, 2)
  
  dimnames(PHI_d) <- list(Fnames_d, Fnames_d)
  
  # dyad-level indicator residual cov matrix
  PSI_d <- diag(6)
  PSI_d[1,1] <- PSI_d[2,2] <- 0.3
  PSI_d[3,3] <- PSI_d[4,4] <- 0.5
  PSI_d[3,4] <- PSI_d[4,3] <- 0.1
  PSI_d[5,5] <- PSI_d[6,6] <- 0.4
  PSI_d[5,6] <- PSI_d[6,5] <- -0.2
  
  dimnames(PSI_d) <- list(Inames_d, Inames_d)
  
  # person-level sigma
  SIGMA_c <- LAM_c %*% PHI_c %*% t(LAM_c) + PSI_c
  # dyad-level sigma
  SIGMA_d <- LAM_d %*% PHI_d %*% t(LAM_d) + PSI_d
  
  ## pop values for satmod
  FsatNames_c <- c("f1@A", "f2@A", "f3@A", "f1@P", "f2@P", "f3@P") # person-level satmod factor names
  FsatNames_d <- c("f1@AP", "f1@PA", "f2@AP", "f2@PA", "f3@AP", "f3@PA") # dyad-level satmod factor names
  
  dimnames(SIGMA_c) <- list(FsatNames_c, FsatNames_c)
  dimnames(SIGMA_d) <- list(FsatNames_d, FsatNames_d)
  
  # person-level corr matrix
  R_c <- cov2cor(SIGMA_c)
  # dyad-level corr matrix
  R_d <- cov2cor(SIGMA_d)
  
  mat_list <- list(LAM_c=LAM_c, PHI_c=PHI_c, PSI_c=PSI_c,
                   LAM_d=LAM_d, PHI_d=PHI_d, PSI_d=PSI_d,
                   SIGMA_c = SIGMA_c, SIGMA_d = SIGMA_d, R_c = R_c, R_d = R_d)
  
  if (return_mats)  return(mat_list)
  
  # creating a dataframe of person-level population values
  pop_c.cov <- as.data.frame(as.table(SIGMA_c))
  colnames(pop_c.cov) <- c("row", "col", "pop.cov")
  pop_c.cov$par_names <- paste0(pop_c.cov$row, "~~", pop_c.cov$col)
  pop_c.cov <- pop_c.cov[, c("par_names", "pop.cov")]
  
  # person-level SDs
  mat_c.SD <- diag(sqrt(diag(SIGMA_c)))
  dimnames(mat_c.SD) <- list(FsatNames_c, FsatNames_c)
  pop_c.SD <- as.data.frame(as.table(mat_c.SD))
  pop_c.SD <- pop_c.SD[pop_c.SD$Freq != 0, ]
  colnames(pop_c.SD) <- c("row", "col", "pop.SD")
  pop_c.SD$par_names <- paste0(pop_c.SD$row, "~~", pop_c.SD$col)
  pop_c.SD <- pop_c.SD[, c("par_names", "pop.SD")]
  rownames(pop_c.SD) <- NULL
  
  # person-level correlations
  # mat_c.cor <- solve(mat_c.SD) %*% SIGMA_c %*% solve(mat_c.SD) # use R_c instead
  pop_c.cor <- as.data.frame(as.table(R_c))
  colnames(pop_c.cor) <- c("row", "col", "pop.cor")
  pop_c.cor$par_names <- paste0(pop_c.cor$row, "~~", pop_c.cor$col)
  pop_c.cor <- pop_c.cor[, c("par_names", "pop.cor")]
  
  # creating a dataframe of dyad-level population values
  pop_d.cov <- as.data.frame(as.table(SIGMA_d))
  colnames(pop_d.cov) <- c("row", "col", "pop.cov")
  pop_d.cov$par_names <-  paste0(pop_d.cov$row, "~~", pop_d.cov$col)
  pop_d.cov <- pop_d.cov[, c("par_names", "pop.cov")]
  
  # dyad-level SDs
  mat_d.SD <- diag(sqrt(diag(SIGMA_d)))
  dimnames(mat_d.SD) <- list(FsatNames_d, FsatNames_d)
  pop_d.SD <- as.data.frame(as.table(mat_d.SD))
  pop_d.SD <- pop_d.SD[pop_d.SD$Freq != 0, ]
  pop_d.SD <- pop_d.SD[!duplicated(pop_d.SD$Freq), ]
  colnames(pop_d.SD) <- c("row", "col", "pop.SD")
  pop_d.SD$par_names <- paste0(pop_d.SD$row, "~~", pop_d.SD$col)
  pop_d.SD <- pop_d.SD[, c("par_names", "pop.SD")]
  rownames(pop_d.SD) <- NULL
  
  # dyad-level correlations
  # mat_d.cor <- solve(mat_d.SD) %*% SIGMA_d %*% solve(mat_d.SD) # use R_d instead
  pop_d.cor <- as.data.frame(as.table(R_d))
  colnames(pop_d.cor) <- c("row", "col", "pop.cor")
  pop_d.cor$par_names <- paste0(pop_d.cor$row, "~~", pop_d.cor$col)
  pop_d.cor <- pop_d.cor[, c("par_names", "pop.cor")]
  
  return(list(pop.cov = rbind(pop_c.cov, pop_d.cov),
              pop.cor = rbind(pop_c.cor, pop_d.cor),
              pop.SD = rbind(pop_c.SD, pop_d.SD)))
}

# getSigma(return_mats = FALSE)
# getSigma(return_mats = TRUE)

#----

# function 1: simulate data for a single group----

genData <- function(n) {
  
  SIGMA_c <- getSigma()$SIGMA_c # save person- and dyad-level pop.cov matrices
  SIGMA_d <- getSigma()$SIGMA_d
  
  Inames_c <- c("V1@A", "V2@A", "V3@A", "V1@P", "V2@P", "V3@P") # person-level indicator names
  Inames_d <- c("V1@AP", "V1@PA", "V2@AP", "V2@PA", "V3@AP", "V3@PA") # dyad-level indicator names
  
  # changing names to generate data per SRM component
  dimnames(SIGMA_c) <- list(Inames_c, Inames_c)
  dimnames(SIGMA_d) <- list(Inames_d, Inames_d)
  
  ##MEAN VECTOR-----
  mu <- rep(0, 6) # group-mean centered SRM component variables
  
  ##DATA GENERATION-----
  
  library(mnormt) # for rmnorm()
  
  dat_c <- rmnorm(n = n, mean = mu, varcov = SIGMA_c)
  dat_d <- rmnorm(n = (n*(n - 1))/2, mean = mu, varcov = SIGMA_d) 
  # Ndyads = n * (n - 1) -- and we already have separate AP and PA columns, 
  # so we need only (n*(n - 1))/2 rows in total
  
  Ydat <- NULL
  Vnames <- c("V1", "V2", "V3")
  
  for(v in Vnames) {
    
    # actor effects
    Aname <- paste0(v, "@A")
    Amat <- matrix(dat_c[, Aname], nrow = n, ncol = n, byrow = FALSE)
    diag(Amat) <- NA
    
    # partner effects
    Pname <- paste0(v, "@P")
    Pmat <- matrix(dat_c[, Pname], nrow = n, ncol = n, byrow = TRUE)
    diag(Pmat) <- NA
    
    # dyad-level effects
    Rmat <- matrix(NA, n, n)
    
    # _ij & _ji effects
    APname <- paste0(v, "@AP")
    PAname <- paste0(v, "@PA")
    
    # build the dyad-level matrix
    foo <- which(lower.tri(Rmat, diag = FALSE), arr.ind = TRUE)
    
    Rmat[foo] <- dat_d[, APname]
    Rmat[foo[, 2:1]] <- dat_d[, PAname]
    
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

thoughtful_priors <- function(targetCorr, precision, default_prior) {
  
  priors <- default_prior
  
  # for SDs --- t-priors
  ## m for SDs is already thoughtful/approximate location because lavaan.srm estimates
  ## are based on the data
  priors$rr_in_t$sd <- priors$rr_out_t$sd <- 
    priors$rr_rel_t$sd <- rep(precision, times = 3) # set SD of t-priors to 0.1
  
  # for correlations --- beta priors
  ## same thoughtful targetCorr for all correlations at both case and dyad level
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

# thoughtful_priors(targetCorr = 0.3, precision = 0.1,
#                   default_prior = lavaan.srm::srm_priors(data = rr.data[c("V1", "V2", "V3")]))

#----

# function 6: prophetic priors----

prophetic_priors <- function(pop_corMat, pop_SDvec, precision, default_prior) {
  
  priors <- default_prior
  
  # for SDs --- t-priors
  popSD_c <- pop_SDvec$popSD_c
  popSD_d <- pop_SDvec$popSD_d
  
  priors$rr_in_t$m <- popSD_c[4:6] # incoming 
  priors$rr_out_t$m <- popSD_c[1:3] # outgoing
  priors$rr_rel_t$m <- popSD_d[c(1, 3, 5)] # relationship
  
  priors$rr_in_t$sd <- priors$rr_out_t$sd <- 
    priors$rr_rel_t$sd <- rep(precision, times = 3) # set SD of t-priors to 0.1
  
  # for correlations --- beta priors
  ## begin: case-level hyperpars----
  popCorr_c <- pop_corMat$popCorr_c
  Fnames_c <- c("f1@A", "f1@P", "f2@A", "f2@P", "f3@A", "f3@P") # reorder names
  popCorr_c <- popCorr_c[Fnames_c, Fnames_c]
  
  # save  population correlation values
  targetVals_c <- list(gen1 = popCorr_c[2,1], ee21 = popCorr_c[3,1], 
                      ae21 = popCorr_c[4,1], ee31 = popCorr_c[5,1], 
                      ae31 = popCorr_c[6,1], ea21 = popCorr_c[3,2], 
                      aa21 = popCorr_c[4,2], ea31 = popCorr_c[5,2], 
                      aa31 = popCorr_c[6,2], gen2 = popCorr_c[4,3], 
                      ee32 = popCorr_c[5,3], ae32 = popCorr_c[6,3],
                      ea32 = popCorr_c[5,4], aa32 = popCorr_c[6,4],
                      gen3 = popCorr_c[6,5])
  
  newHyperpars_c <- mapply(function(tc, ...) optim(par = c(1.5, 1.5), fn = minLoss, ...)$par, 
                       targetCorr = targetVals_c, accept.dev = precision,
                       method = "L-BFGS-B", lower = 0, SIMPLIFY = FALSE)
  
  # save the new (strongly informative) alpha and beta parameters
  alpha_c <- sapply(newHyperpars_c, "[[", 1)
  beta_c <- sapply(newHyperpars_c, "[[", 2)
  
  # populate alpha hyperpars
  priors$case_beta[lower.tri(priors$case_beta, diag = FALSE)] <- alpha_c
  
  # populate beta parameters
  betamat_c <- matrix(NA, 6, 6)
  betamat_c[lower.tri(betamat_c, diag = FALSE)] <- beta_c
  betamat_c <- t(betamat_c)
  priors$case_beta[upper.tri(priors$case_beta, diag = FALSE)] <- 
    betamat_c[upper.tri(betamat_c, diag = FALSE)]
  ## end: case-level hyperpars----
  
  ## begin: dyad-level hyperpars----
  popCorr_d <- pop_corMat$popCorr_d
  
  # saving the population correlation values
  ## DIAG == dyadic reciprocity -- [2,1]; [4,3]; [6,5]
  ## ABOVE = inter -- [3,2]; [5,2]; [5,4]
  ## BELOW = intra -- [4,2]; [6,2]; [6,4]
  targetVals_d <- list(dyad1 = popCorr_d[2,1], 
                      dyad2 = popCorr_d[4,3], 
                      dyad3 = popCorr_d[6,5],
                      inter21 = popCorr_d[3,2], 
                      inter31 = popCorr_d[5,2], 
                      inter32 = popCorr_d[5,4],
                      intra21 = popCorr_d[4,2], 
                      intra31 = popCorr_d[6,2], 
                      intra32 = popCorr_d[6,4])
  
  newHyperpars_d <- mapply(function(tc, ...) optim(par = c(1.5, 1.5), fn = minLoss, ...)$par, 
                       targetCorr = targetVals_d, accept.dev = precision,
                       method = "L-BFGS-B", lower = 0, SIMPLIFY = FALSE)
  
  # new hyperpars for dyadic reciprocity, interpersonal corr, and intrapersonal corr
  dyad_alpha <- c(newHyperpars_d$dyad1[1], newHyperpars_d$dyad2[1], newHyperpars_d$dyad3[1])
  dyad_beta <- c(newHyperpars_d$dyad1[2], newHyperpars_d$dyad2[2], newHyperpars_d$dyad3[2])
  
  inter_alpha <- c(newHyperpars_d$inter21[1], newHyperpars_d$inter31[1], newHyperpars_d$inter32[1])
  inter_beta <- c(newHyperpars_d$inter21[2], newHyperpars_d$inter31[2], newHyperpars_d$inter32[2])
  
  intra_alpha <- c(newHyperpars_d$intra21[1], newHyperpars_d$intra31[1], newHyperpars_d$intra32[1])
  intra_beta <- c(newHyperpars_d$intra21[2], newHyperpars_d$intra31[2], newHyperpars_d$intra32[2])
  
  # replace in priors tables
  diag(priors$rr_beta_a) <- dyad_alpha
  diag(priors$rr_beta_b) <- dyad_beta # dyadic reciprocities on diagonal
  
  priors$rr_beta_a[upper.tri(priors$rr_beta_a)] <- inter_alpha
  priors$rr_beta_b[upper.tri(priors$rr_beta_b)] <- inter_beta # inter = above
  
  priors$rr_beta_a[lower.tri(priors$rr_beta_a)] <- intra_alpha
  priors$rr_beta_b[lower.tri(priors$rr_beta_b)] <- intra_beta # intra = below
  
  ## end: dyad-level hyperpars----
  
  return(priors)
}

# prophetic_priors(pop_corMat = list(popCorr_c = getSigma()$R_c,
#                                                    popCorr_d = getSigma()$R_d),
#                  pop_SDvec = list(popSD_c = sqrt(diag(getSigma()$SIGMA_c)),
#                  popSD_d = sqrt(diag(getSigma()$SIGMA_d))), precision = 0.1,
#                  default_prior = lavaan.srm::srm_priors(data = rr.data[c("V1", "V2", "V3")]))

#----

# function 7: ANOVA-based methods-of-moments priors----

ANOVA_priors <- function(data, rr.vars, IDout, IDin, IDgroup,
                         # modelG = FALSE, modelM = FALSE,
                         #TODO: TDJ: different precisions for SD/Cor?
                         #      optionally use estimated SE as precision?
                         precision, default_prior) {
  rr.data <- data
  priors <- default_prior # default priors
  
  library(TripleR) #FIXME: TDJ: lavaan.srm should check whether it is available
  
  matNames_c <- paste0(rep(rr.vars, each = 2), c("_out", "_in"))
  covMat_c <- matrix(0, length(rr.vars)*2, length(rr.vars)*2,
                     dimnames = list(matNames_c, matNames_c))
  matNames_d <- paste0(rep(rr.vars, each = 2), c("_ij", "_ji"))
  covMat_d <- matrix(0, length(rr.vars)*2, length(rr.vars)*2,
                     dimnames = list(matNames_d, matNames_d))
  # univariate SRM using TripleR()
  for (i in 1:length(rr.vars)) {
    uv.formula <- paste0(rr.vars[i], "~", IDout, "*", IDin)
    if (!missing(IDgroup)) uv.formula <- paste0(uv.formula, "|", IDgroup, collapse = "")
    uvSRM <- RR(as.formula(uv.formula), data = rr.data)
    uv.gEsts <- cbind(uvSRM$varComp.groups[uvSRM$varComp.groups$type != "error variance",]
                      , rr.var = rr.vars[i])
    uv.varNames <- unique(uv.gEsts$type[grep(" variance", uv.gEsts$type)])
    uv.covNames <- unique(uv.gEsts$type[grep("covariance", uv.gEsts$type)])
    
    # rescale negative variances (inadmissable cases) to 0
    uv.gEsts[uv.gEsts$type %in% uv.varNames,]$estimate <- 
      ifelse(uv.gEsts[uv.gEsts$type %in% uv.varNames,]$estimate < 0, 0, 
             uv.gEsts[uv.gEsts$type %in% uv.varNames,]$estimate)
    
    # save variances per variable (weighted mean)
    uv.var <- sapply(uv.varNames, function(type) {
      n <- uv.gEsts[uv.gEsts$type == type,]$group.size
      #TODO: TDJ: needs to import this function from stats using Roxygen
      #TODO: TDJ: weight dyad level by number of dyads?
      weighted.mean(uv.gEsts[uv.gEsts$type == type,]$estimate, w = n-1)
    })
    
    # save covariances per variable (weighted mean)
    uv.cov <- sapply(uv.covNames, function(type) {
      n <- uv.gEsts[uv.gEsts$type == type,]$group.size
      weighted.mean(uv.gEsts[uv.gEsts$type == type,]$estimate, w = n-1)
    })
    
    # construct covariance matrices (to compute correlation matrices)
    covMat_c[paste0(rr.vars[i], "_out"), paste0(rr.vars[i], "_out")] <- uv.var["actor variance"]
    covMat_c[paste0(rr.vars[i], "_in"), paste0(rr.vars[i], "_in")] <- uv.var["partner variance"]
    covMat_c[paste0(rr.vars[i], "_out"), paste0(rr.vars[i], "_in")] <- covMat_c[paste0(rr.vars[i], "_in"), paste0(rr.vars[i], "_out")] <- 
      uv.cov["actor-partner covariance"]
    
    covMat_d[paste0(rr.vars[i], "_ij"), paste0(rr.vars[i], "_ij")] <- 
      covMat_d[paste0(rr.vars[i], "_ji"), paste0(rr.vars[i], "_ji")] <- uv.var["relationship variance"]
    covMat_d[paste0(rr.vars[i], "_ij"), paste0(rr.vars[i], "_ji")] <- 
      covMat_d[paste0(rr.vars[i], "_ji"), paste0(rr.vars[i], "_ij")] <- uv.cov["relationship covariance"]
    
    ## SD priors
    priors$rr_in_t[rr.vars[i],  "m"] <- sqrt(uv.var["partner variance"])
    priors$rr_out_t[rr.vars[i], "m"] <- sqrt(uv.var["actor variance"])
    priors$rr_rel_t[rr.vars[i], "m"] <- sqrt(uv.var["relationship variance"])
    
    priors$rr_in_t[rr.vars[i],  "sd"] <- priors$rr_out_t[rr.vars[i], "sd"] <- priors$rr_rel_t[rr.vars[i], "sd"] <- precision
  }
  
  # bivariate SRM using TripleR()
  for (i in 2:length(rr.vars)) {
    for (j in 1:(i-1)) {
      bv.formula <- paste0(rr.vars[i], "+", rr.vars[j], "~", IDout, "*", IDin)
      if (!missing(IDgroup)) bv.formula <- paste0(bv.formula, "|", IDgroup, collapse = "")
      bvSRM <- RR(as.formula(paste0(bv.formula)), data = rr.data)
      bv.gEsts <- do.call("rbind", lapply(1:length(bvSRM$groups), 
                                          function (g) cbind(bvSRM$groups[[g]]$bivariate, 
                                                             Group = g, 
                                                             rr.var = paste0(rr.vars[i],",",rr.vars[j]))))
      bv.covNames <- unique(bv.gEsts$type[grep("covariance", bv.gEsts$type)])
      
      ## weighted mean of the covariances
      bv.cov <- sapply(bv.covNames, function(type) {
        #FIXME: TDJ: Observed "n" in data might be too large if listwise deleted by TripleR.
        #       Can TripleR add $group.size to result? Why is it missing for bivariate?
        #TODO: TDJ: weight dyad level by number of dyads?
        n <- sapply(1:length(unique(rr.data$Group)), 
                    function(x) length(unique(rr.data[rr.data$Group == x,]$Actor)), 
                    simplify = TRUE)
        weighted.mean(bv.gEsts[bv.gEsts$type == type,]$estimate, w = n-1)
      })
      
      covMat_c[paste0(rr.vars[i], "_out"), paste0(rr.vars[j], "_out")] <- 
        covMat_c[paste0(rr.vars[j], "_out"), paste0(rr.vars[i], "_out")] <- 
        bv.cov["actor-actor covariance"]
      covMat_c[paste0(rr.vars[i], "_in"), paste0(rr.vars[j], "_in")] <- 
        covMat_c[paste0(rr.vars[j], "_in"), paste0(rr.vars[i], "_in")] <- 
        bv.cov["partner-partner covariance"]
      covMat_c[paste0(rr.vars[i], "_out"), paste0(rr.vars[j], "_in")] <- 
        covMat_c[paste0(rr.vars[j], "_in"), paste0(rr.vars[i], "_out")] <- 
        bv.cov["actor-partner covariance"]
      covMat_c[paste0(rr.vars[i], "_in"), paste0(rr.vars[j], "_out")] <- 
        covMat_c[paste0(rr.vars[j], "_out"), paste0(rr.vars[i], "_in")] <- 
        bv.cov["partner-actor covariance"]
      
      covMat_d[paste0(rr.vars[i], "_ij"), paste0(rr.vars[j], "_ij")] <- 
        covMat_d[paste0(rr.vars[j], "_ij"), paste0(rr.vars[i], "_ij")] <- 
        covMat_d[paste0(rr.vars[i], "_ji"), paste0(rr.vars[j], "_ji")] <- 
        covMat_d[paste0(rr.vars[j], "_ji"), paste0(rr.vars[i], "_ji")] <- 
        bv.cov["intrapersonal relationship covariance"]
      
      covMat_d[paste0(rr.vars[i], "_ij"), paste0(rr.vars[j], "_ji")] <- 
        covMat_d[paste0(rr.vars[j], "_ji"), paste0(rr.vars[i], "_ij")] <- 
        covMat_d[paste0(rr.vars[i], "_ji"), paste0(rr.vars[j], "_ij")] <- 
        covMat_d[paste0(rr.vars[j], "_ij"), paste0(rr.vars[i], "_ji")] <- 
        bv.cov["interpersonal relationship covariance"]
    }
  }
  
  corMat_c <- cov2cor(covMat_c)
  corMat_d <- cov2cor(covMat_d)
  
  # if case-level correlation > abs(0.9) then rescale to +/-0.9 and populate srm_priors matrices
  for (x in 2:nrow(corMat_c)) {
    for (y in 1:(x-1)) {
      corMat_c[x,y] <- ifelse(corMat_c[x,y]> 0.9, 0.9, 
                              ifelse(corMat_c[x,y] < -0.9, -0.9, corMat_c[x,y]))
      hyperpars_c <- optim(par = c(1.5, 1.5), fn = minLoss,
                           targetCorr = corMat_c[x,y], accept.dev = precision,
                           method = "L-BFGS-B", lower = 0)$par
      priors$case_beta[x,y] <- hyperpars_c[1] # lower triangle
      priors$case_beta[y,x] <- hyperpars_c[2] # upper triangle
    }
  }
  
  # if dyad-level correlation > abs(0.9) then rescale to +/-0.9, then populate srm_priors matrices 
  corMat_d[row(corMat_d)!=col(corMat_d)] <- ifelse(corMat_d[row(corMat_d)!=col(corMat_d)] > 0.9, 0.9,
                                                   ifelse(corMat_d[row(corMat_d)!=col(corMat_d)] < -0.9,
                                                          -0.9, corMat_d[row(corMat_d)!=col(corMat_d)]))
  
  for (rr in 1:length(rr.vars)) {
    
    ### dyadic-- diagonal
    dyadic_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                              targetCorr = corMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[rr], "_ji")],
                              accept.dev = precision, 
                              method = "L-BFGS-B", lower = 0)$par
    priors$rr_beta_a[rr,rr] <- dyadic_hyperpars[1]
    priors$rr_beta_b[rr,rr] <- dyadic_hyperpars[2]
    
    if (rr > 1L) for (kk in 1:(rr-1)) { 
      ### intra-- below diagonal
      intra_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                               targetCorr = corMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[kk], "_ij")],
                               accept.dev = precision, 
                               method = "L-BFGS-B", lower = 0)$par
      priors$rr_beta_a[rr,kk] <- intra_hyperpars[1]
      priors$rr_beta_b[rr,kk] <- intra_hyperpars[2]
      ### inter-- above diagonal
      inter_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                               targetCorr = corMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[kk], "_ji")],
                               accept.dev = precision, 
                               method = "L-BFGS-B", lower = 0)$par
      priors$rr_beta_a[kk,rr] <- inter_hyperpars[1]
      priors$rr_beta_b[kk,rr] <- inter_hyperpars[2]
    }
  }
  priors
}

# ANOVA_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), IDout = "Actor",
#              IDin = "Partner", IDgroup = "Group", precision = 0.1,
#              default_prior = lavaan.srm::srm_priors(data = rr.data[c("V1", "V2", "V3")]))

#----

# function 8: FIML-based priors----

FIML_priors <- function(data, rr.vars, IDout, IDin, IDgroup, precision = NULL,
                        multi = FALSE, default_prior) {
  
  rr.data <- data
  priors <- default_prior # default priors
  
  library(srm)
  library(car)
  
  #README SEs as prior precisions only implemented for `multi = FALSE`. not necessary 
  # to implement for `multi = TRUE`, because we are not exploring that anymore. 
  
  if (multi) {
    
    ## case syntax
    rr_c <- paste0(rep(rr.vars, each = 2), c("@A", "@P"))
    rrF_c <- paste0("f_", rr_c)
    rrFLoad_c <- paste0(rrF_c, "=~1*", rr_c, "\n") # ULI constraint - default in `srm`
    rrIcov_c <- c(paste0(rr_c, "~~0*", rr_c, "\n"),
                  paste0(rr_c[grep("@A", rr_c)], "~~0*", rr_c[grep("@P", rr_c)], "\n")) # residual (co)v
    rrFvar_c <- paste0(rrF_c, "~~", rrF_c, "\n") # F vars (A/P var)
    rrFgen_c <- paste0("f_", rr.vars, "@A~~gen", rr.vars, "*f_", rr.vars, "@P\n") # generalised
    rrFcov_c <- combn(rr.vars, m = 2, FUN = function(x) {
      rrFAA_c <- paste0("f_", x[1], "@A~~AA", x[1], x[2], "*f_", x[2], "@A\n") # AA cov
      rrFAP_c <- paste0("f_", x[1], "@A~~AP", x[1], x[2], "*f_", x[2], "@P\n") # AP cov
      rrFPA_c <- paste0("f_", x[1], "@P~~PA", x[1], x[2], "*f_", x[2], "@A\n") # PA cov
      rrFPP_c <- paste0("f_", x[1], "@P~~PP", x[1], x[2], "*f_", x[2], "@P\n") # PP cov
      
      c(rrFAA_c, rrFAP_c, rrFPA_c, rrFPP_c)
    })
    mv.syn_c <- c("%Person\n", rrFLoad_c, rrIcov_c, rrFvar_c,
                  rrFgen_c, rrFcov_c)
    
    ## dyad syntax
    rr_d <- paste0(rep(rr.vars, each = 2), c("@AP", "@PA"))
    rrF_d <- paste0("f_", rr_d)
    rrFLoad_d <- paste0(rrF_d, "=~1*", rr_d, "\n") # ULI constraint
    rrIcov_d <- c(paste0(rr_d, "~~0*", rr_d, "\n"),
                  paste0(rr_d[grep("@AP", rr_d)], "~~0*", rr_d[grep("@PA", rr_d)], "\n")) # residual (co)v
    rrFvarAP_d <- paste0("f_", rr.vars, "@AP~~var", rr.vars, "*f_", rr.vars, "@AP\n") # AP vars
    rrFvarPA_d <- paste0("f_", rr.vars, "@PA~~var", rr.vars, "*f_", rr.vars, "@PA\n") # PA vars
    rrFdyad_d <- paste0("f_", rr.vars, "@AP~~dyad", rr.vars, "*f_", rr.vars, "@PA\n") # dyadic
    rrFintra.inter_d <- combn(rr.vars, m = 2, FUN = function(x) {
      rrFintra1_d <- paste0("f_", x[1], "@AP~~intra", x[1], x[2], "*f_", x[2], "@AP\n") # intra1
      rrFintra2_d <- paste0("f_", x[1], "@PA~~intra", x[1], x[2], "*f_", x[2], "@PA\n") # intra2
      rrFinter1_d <- paste0("f_", x[1], "@AP~~inter", x[1], x[2], "*f_", x[2], "@PA\n") # inter1
      rrFinter2_d <- paste0("f_", x[1], "@PA~~inter", x[1], x[2], "*f_", x[2], "@AP\n") # inter2
      c(rrFintra1_d, rrFintra2_d, rrFinter1_d, rrFinter2_d)
    }) # intra + inter
    mv.syn_d <- c("%Dyad\n", rrFLoad_d, rrIcov_d, rrFvarAP_d, rrFvarPA_d, rrFdyad_d,
                  rrFintra.inter_d)
    
    mv.syn <- paste0(c(mv.syn_c, mv.syn_d))
    
    mv.fit <- srm(model.syntax = mv.syn, data = rr.data, rrgroup_name = IDgroup,
                  person_names = c(IDout, IDin), 
                  fixed.groups = FALSE,
                  verbose = FALSE)
    
    ## populate case-level hyperparameters
    MIcov_c <- mv.fit$sigma$U[[1]] # case-level model-implied covmat
    MIdimNames_c <- c(rrF_c[grep("@A", rrF_c)], rrF_c[grep("@P", rrF_c)]) # `srm` stores in a different order
    dimnames(MIcov_c) <- list(MIdimNames_c, MIdimNames_c)
    MIcov_c <- MIcov_c[rrF_c, rrF_c] # rearrange to same order as srm_priors() function
    for(x in 1:nrow(MIcov_c)) {
      for (y in 1:ncol(MIcov_c)) {
        if (x == y) {
          next
        } else {
          if (MIcov_c[x,x] > 0 && MIcov_c[y,y] > 0) {
            MIcor_c <- MIcov_c[x,y] / sqrt(MIcov_c[x,x]*MIcov_c[y,y]) 
            MIcor_c <- ifelse(MIcor_c > 0.9, 0.9, 
                              ifelse(MIcor_c < -0.9, -0.9, MIcor_c)) # rescale large correlations
            hyperpars_c <- optim(par = c(1.5, 1.5), fn = minLoss,
                                 targetCorr = MIcor_c, accept.dev = precision,
                                 method = "L-BFGS-B", lower = 0)$par # alpha and beta parameters
            
            priors$case_beta[x,y] <- hyperpars_c[1] # lower triangle
            priors$case_beta[y,x] <- hyperpars_c[2] # upper triangle
          } else {
            priors$case_beta[x,y] <- priors$case_beta[y,x] <- 1.5
          } 
        }
      }
    }
    diag(MIcov_c) <- ifelse(diag(MIcov_c) < 0, 0, diag(MIcov_c)) # rescale negative variances before assigning to SD t-priors
    MISDs_c <- sqrt(diag(MIcov_c))
    priors$rr_in_t$m <- MISDs_c[grep("@P", names(MISDs_c))]
    priors$rr_out_t$m <- MISDs_c[grep("@A", names(MISDs_c))]
    priors$rr_in_t$sd <- priors$rr_out_t$sd <- precision # incoming effects priors
    
    ## populate dyad-level hyperparameters
    MIcov_d <- mv.fit$sigma$D[[1]] # dyad-level model-implied covmat
    MIdimNames_d <- gsub("@AP", "_ij", gsub("@PA", "_ji", gsub("f_", "", rrF_d)))
    dimnames(MIcov_d) <- list(MIdimNames_d, MIdimNames_d) # change names---make assignment to prior mats easy
    
    for (rr in 1:length(rr.vars)) {
      if (MIcov_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[rr], "_ij")] <= 0) { #FIXME--- am i doing the ifelse statements incorrectly?
        priors$rr_beta_a[rr,rr] <- priors$rr_beta_b[rr,rr] <- 1.5
      } else {
        ### dyadic-- diagonal
        MIdyadic_d <- MIcov_d[paste0(rr.vars[rr], "_ij"), 
                              paste0(rr.vars[rr], "_ji")] / MIcov_d[paste0(rr.vars[rr], "_ij"), 
                                                                    paste0(rr.vars[rr], "_ij")]
        MIdyadic_d <- ifelse(MIdyadic_d > 0.9, 0.9, 
                             ifelse(MIdyadic_d < -0.9, -0.9, MIdyadic_d)) # rescale large correlations
        dyadic_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                  targetCorr = MIdyadic_d,
                                  accept.dev = precision, 
                                  method = "L-BFGS-B", lower = 0)$par
        priors$rr_beta_a[rr,rr] <- dyadic_hyperpars[1]
        priors$rr_beta_b[rr,rr] <- dyadic_hyperpars[2]
        
        if (rr > 1L) for (kk in 1:(rr-1)) {
          if (MIcov_d[paste0(rr.vars[kk], "_ij"), paste0(rr.vars[kk], "_ij")] <= 0 | 
              MIcov_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[rr], "_ij")] <= 0) {
            priors$rr_beta_a[rr,kk] <- priors$rr_beta_b[rr,kk] <- 
              priors$rr_beta_a[kk,rr] <- priors$rr_beta_b[kk,rr] <- 1.5
          } else {
            ### intra-- below diagonal
            MIintra_d <- MIcov_d[paste0(rr.vars[rr], "_ij"), 
                                 paste0(rr.vars[kk], "_ij")] / sqrt(MIcov_d[paste0(rr.vars[rr], "_ij"), 
                                                                            paste0(rr.vars[rr], "_ij")] * 
                                                                      MIcov_d[paste0(rr.vars[kk], "_ij"), 
                                                                              paste0(rr.vars[kk], "_ij")])
            MIintra_d <- ifelse(MIintra_d > 0.9, 0.9, 
                                ifelse(MIintra_d < -0.9, -0.9, MIintra_d)) # rescale large correlations
            intra_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                     targetCorr = MIintra_d,
                                     accept.dev = precision, 
                                     method = "L-BFGS-B", lower = 0)$par
            priors$rr_beta_a[rr,kk] <- intra_hyperpars[1]
            priors$rr_beta_b[rr,kk] <- intra_hyperpars[2]
            ### inter-- above diagonal
            MIinter_d <- MIcov_d[paste0(rr.vars[rr], "_ij"), 
                                 paste0(rr.vars[kk], "_ji")] / sqrt(MIcov_d[paste0(rr.vars[rr], "_ij"), 
                                                                            paste0(rr.vars[rr], "_ij")] * 
                                                                      MIcov_d[paste0(rr.vars[kk], "_ij"), 
                                                                              paste0(rr.vars[kk], "_ij")])
            MIinter_d <- ifelse(MIinter_d > 0.9, 0.9, 
                                ifelse(MIinter_d < -0.9, -0.9, MIinter_d)) # rescale large correlations
            inter_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                     targetCorr = MIinter_d,
                                     accept.dev = precision, 
                                     method = "L-BFGS-B", lower = 0)$par
            priors$rr_beta_a[kk,rr] <- inter_hyperpars[1]
            priors$rr_beta_b[kk,rr] <- inter_hyperpars[2]
          }
        }
      }
    }
    diag(MIcov_d) <- ifelse(diag(MIcov_d) < 0, 0, diag(MIcov_d)) # rescale negative variances before assigning to SD t-priors
    MISDs_d <- sqrt(diag(MIcov_d))
    priors$rr_rel_t$m <- MISDs_d[grep("_ij", names(MISDs_d))]
    priors$rr_rel_t$sd <- precision
  } else {
    
    for (rr in 1:length(rr.vars)) {
      # uvSRM
      uv.syn <- paste0("%Person\n",
                       "f_", rr.vars[rr], "@A=~1*", rr.vars[rr], "@A\n", 
                       "f_", rr.vars[rr], "@P=~1*", rr.vars[rr], "@P\n", # F loadings
                       
                       rr.vars[rr], "@A~~0*", rr.vars[rr], "@A+0*", rr.vars[rr], "@P\n",
                       rr.vars[rr], "@P~~0*", rr.vars[rr], "@P\n", # residual (co)v
                       
                       "f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@A+f_", rr.vars[rr], "@P\n",
                       "f_", rr.vars[rr], "@P~~f_", rr.vars[rr], "@P\n", # A/P var + gen cov
                       
                       "%Dyad\n",
                       "f_", rr.vars[rr], "@AP=~1*", rr.vars[rr], "@AP\n",
                       "f_", rr.vars[rr], "@PA=~1*", rr.vars[rr], "@PA\n", # F loadings
                       
                       rr.vars[rr], "@AP~~0*", rr.vars[rr], "@AP+0*", rr.vars[rr], "@PA\n",
                       rr.vars[rr], "@PA~~0*", rr.vars[rr], "@PA\n", # residual (co)v
                       
                       "f_", rr.vars[rr], "@AP~~var", rr.vars[rr], "*f_", rr.vars[rr], "@AP+f_", rr.vars[rr], "@PA\n",
                       "f_", rr.vars[rr], "@PA~~var", rr.vars[rr], "*f_", rr.vars[rr], "@PA") # AP/PA var + dyadic cov
      
      uv.fit <- srm(model.syntax = uv.syn, data = rr.data, rrgroup_name = IDgroup,
                    person_names = c(IDout, IDin), 
                    fixed.groups = FALSE,
                    verbose = FALSE)
      
      MIuvcov_c <- uv.fit$sigma$U[[1]]
      MIuvdimNames_c <- paste0(rep(rr.vars[rr], each = 2), c("@A", "@P")) 
      dimnames(MIuvcov_c) <- list(MIuvdimNames_c, MIuvdimNames_c)
      if (MIuvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] >= 0) { # if actor var > 0
        AvarDelta_c <- deltaMethod(uv.fit, g. = paste0("sqrt(`f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`)"))
        priors$rr_out_t[rr.vars[rr], "m"] <- AvarDelta_c$Estimate
        priors$rr_out_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, AvarDelta_c$SE)
      } else {
        priors$rr_out_t[rr.vars[rr], "m"] <- 0
        priors$rr_out_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision,
                                                     uv.fit$parm.table[uv.fit$parm.table$par_names == 
                                                                         paste0("f_", rr.vars[rr],
                                                                                "@A~~f_", rr.vars[rr], "@A"), "se"])
      }
      if (MIuvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] >= 0) { # if partner var > 0
        PvarDelta_c <- deltaMethod(uv.fit, g. = paste0("sqrt(`f_", rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`)"))
        priors$rr_in_t[rr.vars[rr], "m"] <- PvarDelta_c$Estimate
        priors$rr_in_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, PvarDelta_c$SE)
      } else {
        priors$rr_in_t[rr.vars[rr], "m"] <- 0
        priors$rr_in_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, 
                                                    uv.fit$parm.table[uv.fit$parm.table$par_names == 
                                                                        paste0("f_", rr.vars[rr],
                                                                               "@P~~f_", rr.vars[rr], "@P"), "se"])
      }
      if (MIuvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] > 0 && 
          MIuvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] > 0) { # if both actor & partner var >0
        
        ## generalised cor
        uvcorDelta_c <- deltaMethod(uv.fit, g. = paste0("`f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@P`/(sqrt(`f_", 
                                                        rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`*`f_", 
                                                        rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`))"))
        MIuvcor_c <- uvcorDelta_c$Estimate
        MIuvcorSE_c <- uvcorDelta_c$SE
        MIuvcor_c <- ifelse(MIuvcor_c > 0.9, 0.9, ifelse(MIuvcor_c < -0.9, -0.9, MIuvcor_c))
        genR_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = MIuvcor_c, 
                                accept.dev = ifelse(!is.null(precision), precision, MIuvcorSE_c),
                                method = "L-BFGS-B", lower = 0)$par
        priors$case_beta[paste0(rr.vars[rr],"_in"),paste0(rr.vars[rr],"_out")] <- genR_hyperpars[1] # lower triangle
        priors$case_beta[paste0(rr.vars[rr],"_out"),paste0(rr.vars[rr],"_in")] <- genR_hyperpars[2] # upper triangle
        
      } else {
        priors$case_beta[paste0(rr.vars[rr],"_in"),
                         paste0(rr.vars[rr],"_out")] <- priors$case_beta[paste0(rr.vars[rr],"_out"),
                                                                         paste0(rr.vars[rr],"_in")] <- 1.5
      }
      
      MIuvcov_d <- uv.fit$sigma$D[[1]]
      MIuvdimNames_d <- paste0(rep(rr.vars[rr], each = 2), c("@AP", "@PA"))
      dimnames(MIuvcov_d) <- list(MIuvdimNames_d, MIuvdimNames_d)
      if (MIuvcov_d[paste0(rr.vars[rr], "@AP"), paste0(rr.vars[rr], "@AP")] > 0) { # if rel var > 0
        uvcorDelta_d <- deltaMethod(uv.fit, g. = paste0("`f_", rr.vars[rr], "@AP~~f_", rr.vars[rr], "@PA`/`f_", 
                                                        rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`"))
        MIuvcor_d <- uvcorDelta_d$Estimate
        MIuvcorSE_d <- uvcorDelta_d$SE
        MIuvcor_d <- ifelse(MIuvcor_d > 0.9, 0.9, ifelse(MIuvcor_d < -0.9, -0.9, MIuvcor_d))
        dyadicR_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                   targetCorr = MIuvcor_d, 
                                   accept.dev = ifelse(!is.null(precision), precision, MIuvcorSE_d),
                                   method = "L-BFGS-B", lower = 0)$par 
        priors$rr_beta_a[rr,rr] <- dyadicR_hyperpars[1]
        priors$rr_beta_b[rr,rr] <- dyadicR_hyperpars[2]
        
        relvarDelta_c <- deltaMethod(uv.fit, g. = paste0("sqrt(`f_", rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`)"))
        priors$rr_rel_t[rr.vars[rr], "m"] <- relvarDelta_c$Estimate
        priors$rr_rel_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, relvarDelta_c$SE)
      } else {
        priors$rr_beta_a[rr,rr] <- priors$rr_beta_b[rr,rr] <- 1.5
        
        priors$rr_rel_t[rr.vars[rr], "m"] <- 0
        priors$rr_rel_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision,
                                                     uv.fit$parm.table[uv.fit$parm.table$par_names == 
                                                                         paste0("f_", rr.vars[rr],
                                                                                "@AP~~f_", rr.vars[rr], "@AP"), "se"])
      }
      
      if (rr > 1L) for (kk in 1:(rr-1)) {
        # bvSRM
        bv.syn <- paste0("%Person\n",
                         "f_", rr.vars[rr], "@A=~1*", rr.vars[rr], "@A\n",
                         "f_", rr.vars[kk], "@A=~1*", rr.vars[kk], "@A\n",
                         
                         "f_", rr.vars[rr], "@P=~1*", rr.vars[rr], "@P\n",
                         "f_", rr.vars[kk], "@P=~1*", rr.vars[kk], "@P\n",
                         
                         rr.vars[rr], "@A~~0*", rr.vars[rr], "@A+0*", rr.vars[rr], "@P\n",
                         rr.vars[rr], "@P~~0*", rr.vars[rr], "@P\n",
                         rr.vars[kk], "@A~~0*", rr.vars[kk], "@A+0*", rr.vars[kk], "@P\n",
                         rr.vars[kk], "@P~~0*", rr.vars[kk], "@P\n",
                         
                         "f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@A+f_", rr.vars[rr], "@P+f_", rr.vars[kk], "@A+f_", rr.vars[kk], "@P\n",
                         "f_", rr.vars[kk], "@A~~f_", rr.vars[kk], "@A+f_", rr.vars[kk], "@P+f_", rr.vars[rr], "@P\n",
                         "f_", rr.vars[rr], "@P~~f_", rr.vars[rr], "@P+f_", rr.vars[kk], "@P\n",
                         "f_", rr.vars[kk], "@P~~f_", rr.vars[kk], "@P\n",
                         
                         "%Dyad\n",
                         "f_", rr.vars[rr], "@AP=~1*", rr.vars[rr], "@AP\n",
                         "f_", rr.vars[kk], "@AP=~1*", rr.vars[kk], "@AP\n",
                         
                         "f_", rr.vars[rr], "@PA=~1*", rr.vars[rr], "@PA\n",
                         "f_", rr.vars[kk], "@PA=~1*", rr.vars[kk], "@PA\n",
                         
                         rr.vars[rr], "@AP~~0*", rr.vars[rr], "@AP+0*", rr.vars[rr], "@PA\n",
                         rr.vars[rr], "@PA~~0*", rr.vars[rr], "@PA\n",
                         rr.vars[kk], "@AP~~0*", rr.vars[kk], "@AP+0*", rr.vars[kk], "@PA\n",
                         rr.vars[kk], "@PA~~0*", rr.vars[kk], "@PA\n",
                         
                         "f_", rr.vars[rr], "@AP~~var", rr.vars[rr], "*f_", rr.vars[rr], "@AP+dyad", 
                         rr.vars[rr], "*f_", rr.vars[rr], "@PA+intra", rr.vars[rr], 
                         rr.vars[kk], "*f_", rr.vars[kk], "@AP+inter", rr.vars[rr], 
                         rr.vars[kk], "*f_", rr.vars[kk], "@PA\n",
                         "f_", rr.vars[kk], "@AP~~var", rr.vars[kk],"*f_", rr.vars[kk], 
                         "@AP+dyad", rr.vars[kk],"*f_", rr.vars[kk], "@PA+inter", 
                         rr.vars[rr], rr.vars[kk], "*f_", rr.vars[rr], "@PA\n",
                         "f_", rr.vars[rr], "@PA~~var", rr.vars[rr], "*f_", rr.vars[rr], 
                         "@PA+intra", rr.vars[rr], rr.vars[kk], "*f_", rr.vars[kk], "@PA\n",
                         "f_", rr.vars[kk], "@PA~~var", rr.vars[kk],"*f_", rr.vars[kk], "@PA")
        
        bv.fit <- srm(model.syntax = bv.syn, data = rr.data, rrgroup_name = IDgroup,
                      person_names = c(IDout, IDin), 
                      fixed.groups = FALSE,
                      verbose = FALSE)
        
        # case-level hyperparameters
        MIbvcov_c <- bv.fit$sigma$U[[1]]
        MIbvdimNames_c <- c(paste0(c(rr.vars[kk], rr.vars[rr]), "@A"),
                            paste0(c(rr.vars[kk], rr.vars[rr]), "@P"))
        dimnames(MIbvcov_c) <- list(MIbvdimNames_c, MIbvdimNames_c)
        ## AA
        if (MIbvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@A"), paste0(rr.vars[kk], "@A")] > 0) {
          bvAADelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@A~~f_", rr.vars[kk], "@A`/(sqrt(`f_", 
                                                         rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`*`f_", 
                                                         rr.vars[kk], "@A~~f_", rr.vars[kk], "@A`))"))
          AA_c <- bvAADelta_c$Estimate
          AASE_c <- bvAADelta_c$SE
          AA_c <- ifelse(AA_c > 0.9, 0.9, ifelse(AA_c < -0.9, -0.9, AA_c))
          AA_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = AA_c, 
                                accept.dev = ifelse(!is.null(precision), precision, AASE_c),
                                method = "L-BFGS-B", lower = 0)$par ## AA corr
          priors$case_beta[paste0(rr.vars[rr], "_out"), paste0(rr.vars[kk], "_out")] <- AA_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_out"), paste0(rr.vars[rr], "_out")] <- AA_hyperpars[2] # upper triangle
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_out"), 
                           paste0(rr.vars[kk], "_out")] <- priors$case_beta[paste0(rr.vars[kk], "_out"), 
                                                                            paste0(rr.vars[rr], "_out")] <- 1.5
        }
        ## AP
        if (MIbvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@P"), paste0(rr.vars[kk], "@P")] > 0) {
          bvAPDelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@A~~f_", rr.vars[kk], "@P`/(sqrt(`f_", 
                                                         rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`*`f_", 
                                                         rr.vars[kk], "@P~~f_", rr.vars[kk], "@P`))"))
          AP_c <- bvAPDelta_c$Estimate
          APSE_c <- bvAPDelta_c$SE
          AP_c <- ifelse(AP_c > 0.9, 0.9, ifelse(AP_c < -0.9, -0.9, AP_c))
          AP_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = AP_c, 
                                accept.dev = ifelse(!is.null(precision), precision, APSE_c),
                                method = "L-BFGS-B", lower = 0)$par ## AP corr
          priors$case_beta[paste0(rr.vars[rr], "_out"), paste0(rr.vars[kk], "_in")] <- AP_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_in"), paste0(rr.vars[rr], "_out")] <- AP_hyperpars[2] # upper triangle
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_out"), 
                           paste0(rr.vars[kk], "_in")] <- priors$case_beta[paste0(rr.vars[kk], "_in"), 
                                                                           paste0(rr.vars[rr], "_out")] <- 1.5
        }
        ## PA
        if (MIbvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@A"), paste0(rr.vars[kk], "@A")] > 0) {
          bvPADelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[kk], "@A~~f_", rr.vars[rr], "@P`/(sqrt(`f_", 
                                                         rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`*`f_", 
                                                         rr.vars[kk], "@A~~f_", rr.vars[kk], "@A`))"))
          PA_c <- bvPADelta_c$Estimate
          PASE_c <- bvPADelta_c$SE
          PA_c <- ifelse(PA_c > 0.9, 0.9, ifelse(PA_c < -0.9, -0.9, PA_c))
          PA_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = PA_c, 
                                accept.dev = ifelse(!is.null(precision), precision, PASE_c),
                                method = "L-BFGS-B", lower = 0)$par ## PA corr
          priors$case_beta[paste0(rr.vars[rr], "_in"), paste0(rr.vars[kk], "_out")] <- PA_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_out"), paste0(rr.vars[rr], "_in")] <- PA_hyperpars[2] # upper triangle
          
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_in"), 
                           paste0(rr.vars[kk], "_out")] <- priors$case_beta[paste0(rr.vars[kk], "_out"), 
                                                                            paste0(rr.vars[rr], "_in")] <- 1.5
        }
        ## PP
        if (MIbvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@P"), paste0(rr.vars[kk], "@P")] > 0) {
          
          bvPPDelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@P~~f_", rr.vars[kk], "@P`/(sqrt(`f_", 
                                                         rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`*`f_", 
                                                         rr.vars[kk], "@P~~f_", rr.vars[kk], "@P`))"))
          PP_c <- bvPPDelta_c$Estimate
          PPSE_c <- bvPPDelta_c$SE
          PP_c <- ifelse(PP_c > 0.9, 0.9, ifelse(PP_c < -0.9, -0.9, PP_c))
          PP_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = PP_c, 
                                accept.dev = ifelse(!is.null(precision), precision, PPSE_c),
                                method = "L-BFGS-B", lower = 0)$par ## PA corr
          priors$case_beta[paste0(rr.vars[rr], "_in"), paste0(rr.vars[kk], "_in")] <- PP_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_in"), paste0(rr.vars[rr], "_in")] <- PP_hyperpars[2] # upper triangle
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_in"), 
                           paste0(rr.vars[kk], "_in")] <- priors$case_beta[paste0(rr.vars[kk], "_in"), 
                                                                           paste0(rr.vars[rr], "_in")] <- 1.5
        }
        
        # dyad-level hyperparameters
        MIbvcov_d <- bv.fit$sigma$D[[1]]
        MIbvdimNames_d <- c(paste0(rep(rr.vars[kk], times = 2), c("@AP", "@PA")), 
                            paste0(rep(rr.vars[rr], times = 2), c("@AP", "@PA"))) #FIXME maybe just directly specify _ij/ji here?
        dimnames(MIbvcov_d) <- list(MIbvdimNames_d, MIbvdimNames_d)
        if (MIbvcov_d[paste0(rr.vars[rr], "@AP"), paste0(rr.vars[rr], "@AP")] > 0 && 
            MIbvcov_d[paste0(rr.vars[kk], "@AP"), paste0(rr.vars[kk], "@AP")] > 0) {
          ### intra-- below diagonal
          bvintraDelta_d <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@AP~~f_", rr.vars[kk], "@AP`/(sqrt(`f_", 
                                                            rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`*`f_", 
                                                            rr.vars[kk], "@AP~~f_", rr.vars[kk], "@AP`))"))
          intra_d <- bvintraDelta_d$Estimate
          intraSE_d <- bvintraDelta_d$SE
          intra_d <- ifelse(intra_d > 0.9, 0.9, ifelse(intra_d < -0.9, -0.9, intra_d))
          intra_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                   targetCorr = intra_d,
                                   accept.dev = ifelse(!is.null(precision), precision, intraSE_d), 
                                   method = "L-BFGS-B", lower = 0)$par
          priors$rr_beta_a[rr,kk] <- intra_hyperpars[1]
          priors$rr_beta_b[rr,kk] <- intra_hyperpars[2]
          
          ### inter-- above diagonal
          bvinterDelta_d <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[kk], "@AP~~f_", rr.vars[rr], "@PA`/(sqrt(`f_", 
                                                            rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`*`f_", 
                                                            rr.vars[kk], "@AP~~f_", rr.vars[kk], "@AP`))"))
          inter_d <- bvinterDelta_d$Estimate
          interSE_d <- bvinterDelta_d$SE
          inter_d <- ifelse(inter_d > 0.9, 0.9, ifelse(inter_d < -0.9, -0.9, inter_d))
          inter_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                   targetCorr = inter_d,
                                   accept.dev = ifelse(!is.null(precision), precision, interSE_d), 
                                   method = "L-BFGS-B", lower = 0)$par
          priors$rr_beta_a[kk,rr] <- inter_hyperpars[1]
          priors$rr_beta_b[kk,rr] <- inter_hyperpars[2]
        } else {
          priors$rr_beta_a[rr,kk] <- priors$rr_beta_b[rr,kk] <- 
            priors$rr_beta_a[kk,rr] <- priors$rr_beta_b[kk,rr] <- 1.5
        }
      }
    }
  }
  priors
}

# FIML_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), IDout = "Actor",
#              IDin = "Partner", IDgroup = "Group", precision = 0.1, multi = FALSE,
#              default_prior = lavaan.srm::srm_priors(data = rr.data[c("V1", "V2", "V3")]))

#----

# function 9: set customised priors for MCMC stage----

set_priors <- function(data, rr.vars, IDout, IDin, IDgroup, priorType, targetCorr,
                       pop_corMat = list(popCorr_c = getSigma()$R_c,
                                         popCorr_d = getSigma()$R_d),
                       pop_SDvec = list(popSD_c = sqrt(diag(getSigma()$SIGMA_c)),
                                        popSD_d = sqrt(diag(getSigma()$SIGMA_d))),
                       precision, multiMLE = FALSE) {
  prior_env <- new.env()
  prior_env$default_prior <- srm_priors(data = data[rr.vars]) # default MCMC priors (diffuse priors)
 
 if (priorType == "default") { # default (diffuse) priors
   srmPriors <- get("default_prior", envir = prior_env)
 } else if (priorType == "thoughtful") { # thoughtful priors
   srmPriors <- thoughtful_priors(targetCorr = targetCorr, 
                               precision = precision, 
                               default_prior = get("default_prior", envir = prior_env))
   
 } else if (priorType == "prophetic") { # prophetic priors
   srmPriors <- prophetic_priors(pop_corMat = pop_corMat, 
                              pop_SDvec = pop_SDvec, precision = precision, 
                              default_prior = get("default_prior", envir = prior_env))
   
 } else if (priorType == "ANOVA") { # method-of-moments priors (ANOVA-based, `TripleR`)
   srmPriors <- ANOVA_priors(data = data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                             IDgroup = IDgroup, precision = precision,
                             default_prior = get("default_prior", envir = prior_env))
   
 } else if (priorType == "FIML") { # FIML-based priors (`srm`)
   srmPriors <- FIML_priors(data = data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                            IDgroup = IDgroup, precision = precision, multi = multiMLE,
                            default_prior = get("default_prior", envir = prior_env))
 }
 return(srmPriors)
}

# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), priorType = "default")
# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"),priorType = "prophetic",
#            precision = 0.1)
# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), priorType = "thoughtful",
#            targetCorr = 0.3, precision = 0.1)
# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), priorType = "ANOVA",
#            IDout = "Actor", IDin = "Partner", IDgroup = "Group", precision = 0.1)
# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), priorType = "FIML",
#            IDout = "Actor", IDin = "Partner", IDgroup = "Group",
#            multiMLE = FALSE, precision = 0.1)


# README https://stackoverflow.com/questions/18338331/nested-function-environment-selection 
# README https://digitheadslabnotebook.blogspot.com/2011/06/environments-in-r.html
# README http://adv-r.had.co.nz/Environments.html#env-answers 

#----

# function 10: stage 1 in `lavaan.srm`----

s1sat <- function(MCSampID, n, G, rr.vars = c("V1", "V2", "V3"), 
                  IDout = "Actor", IDin = "Partner", IDgroup = "Group", priorType,
                  targetCorr = 0.3, pop_corMat = list(popCorr_c = getSigma()$R_c,
                                    popCorr_d = getSigma()$R_d),
                  pop_SDvec = list(popSD_c = sqrt(diag(getSigma()$SIGMA_c)),
                                   popSD_d = sqrt(diag(getSigma()$SIGMA_d))),
                  precision = 0.1, multiMLE = FALSE, iter = 2000, savefile = FALSE) {
  library(lavaan.srm)
  library(coda) # for gelman.diag()
  library(rstan) # for As.mcmc.list()
  
  s1_env <- new.env()
  s1_env$dat <- genGroups(MCSampID = MCSampID, n = n, G = G)
  s1_env$MCMC_pars <- c("s_rr", "S_p", "r_d2", 
                        paste0("Rp[", combn(1:6, 2, paste0, 
                                            collapse = ","), "]")) # save non-redundant MCMC parameter labels
  t0 <- Sys.time()
  
  if (priorType == "default") { # default
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, priorType = priorType)
    
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    ## compute mPSRF
    mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
    mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
    
    s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                      probs = c(0.025, 0.975))$summary))
    s1long <- s1long[-nrow(s1long), ]
    
    if (mPSRF > 1.05) {
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter*2, priors = s1_priors, seed = 1512, verbose = F)
      
      ## compute mPSRF
      mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
      mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
      
      s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                      probs = c(0.025, 0.975))$summary))
      s1long <- s1long[-nrow(s1long), ]
    }
  } else if (priorType == "thoughtful" && !missing(targetCorr)) { # thoughtful
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, priorType = priorType,
                          targetCorr = targetCorr, precision = precision)
    
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    ## compute mPSRF
    mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
    mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
    
    s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                    probs = c(0.025, 0.975))$summary))
    s1long <- s1long[-nrow(s1long), ]
    
    if (mPSRF > 1.05) {
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter*2, priors = s1_priors, seed = 1512, verbose = F)
      
      ## compute mPSRF
      mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
      mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
      
      s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                      probs = c(0.025, 0.975))$summary))
      s1long <- s1long[-nrow(s1long), ]
    }
  } else if (priorType == "prophetic") { # prophetic
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, priorType = priorType,
                                   precision = precision)
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    ## compute mPSRF
    mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
    mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
    
    s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                    probs = c(0.025, 0.975))$summary))
    s1long <- s1long[-nrow(s1long), ]
    
    if (mPSRF > 1.05) {
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter*2, priors = s1_priors, seed = 1512, verbose = F)
      
      ## compute mPSRF
      mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
      mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
      
      s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                      probs = c(0.025, 0.975))$summary))
      s1long <- s1long[-nrow(s1long), ]
    }
  } else if (priorType == "ANOVA" && !missing(IDout) && !missing(IDin) && !missing(IDgroup)) { # ANOVA
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                               IDgroup = IDgroup, priorType = priorType, precision = precision)
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    ## compute mPSRF
    mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
    mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
    
    s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                    probs = c(0.025, 0.975))$summary))
    s1long <- s1long[-nrow(s1long), ]
    
    if (mPSRF > 1.05) {
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter*2, priors = s1_priors, seed = 1512, verbose = F)
      
      ## compute mPSRF
      mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
      mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
      
      s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                      probs = c(0.025, 0.975))$summary))
      s1long <- s1long[-nrow(s1long), ]
    }
  } else if (priorType == "FIML" && !missing(IDout) && !missing(IDin) && !missing(IDgroup)) { # FIML
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                              IDgroup = IDgroup, priorType = priorType, precision = precision,
                              multi = multiMLE)
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    ## compute mPSRF
    mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
    mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
    
    s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                    probs = c(0.025, 0.975))$summary))
    s1long <- s1long[-nrow(s1long), ]
    
    if (mPSRF > 1.05) {
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter*2, priors = s1_priors, seed = 1512, verbose = F)
      
      ## compute mPSRF
    mcmcList <- As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env))
    mPSRF <- gelman.diag(mcmcList, autoburnin = F)$mpsrf
    
    s1long <- cbind(iter = iter, data.frame(summary(s1ests, as.stanfit = TRUE,
                                                      probs = c(0.025, 0.975))$summary))
    s1long <- s1long[-nrow(s1long), ]
    }
  }
  t1 <- Sys.time()

  #begin: saving results----
  
  # save variable names
  pVarNames <- paste0("f", rep(1:3, each = 2), "@", rep(c("A", "P"), times = 3))
  dVarNames <- paste0("f", rep(1:3, each = 2), "@", rep(c("AP", "PA"), times = 3))
  
  # save IDs for rows that need to be extracted
  p.idx <- outer(1:6, 1:6, FUN = paste, sep = ",")[upper.tri(diag(6), diag = TRUE)]
  d.idx <- c(paste0(rep(1, each = 6), ",", 1:6), # relvar1; dyadcov11; intra/inter12; intra/inter13
             paste0(rep(3, each = 4), ",", 3:6), #relvar2; dyadcov22; intra/inter23
             paste0(rep(5, each = 2), ",", 5:6)) # relvar3; dyadcov33 
  
  # case-level covariances
  pcovElements <- grep(pattern = "pSigma", rownames(s1long))
  pSigma <- s1long[pcovElements, ] # EAPs
  pSigma$row <- rep(pVarNames, each = 6)
  pSigma$col <- rep(pVarNames, times = 6)
  pSigma$par_names <- paste0(pSigma$row, "~~", pSigma$col)
  pSigma <- subset(pSigma, select = -c(row, col))
  pSigmaRows <- paste0("pSigma[", p.idx, "]")
  pSigma <- pSigma[rownames(pSigma) %in% pSigmaRows, ]
  colnames(pSigma) <- c("iter", "Ecov", "Ecov.MCE", "Ecov.SE", 
                        "Ecov.low", "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", 
                        "par_names")
  pSigma$level <- "case"
  
  s1pM.covmat  <- summary(s1ests, component = "case", interval = "hdi",
                          posterior.est = "mode", srm.param = "cov") # MAPs
  s1pM.covEsts <- as.data.frame(as.table(s1pM.covmat$case$cov$mode))
  s1pM.covlow  <- as.data.frame(as.table(s1pM.covmat$group$cov$hdi$lower))
  s1pM.covup   <- as.data.frame(as.table(s1pM.covmat$group$cov$hdi$upper))
  s1pM.covlist <- list(s1pM.covEsts, s1pM.covlow, s1pM.covup)
  s1pM.covtab  <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), s1pM.covlist)
  colnames(s1pM.covtab) <- c("row", "col", "Mcov", "Mcov.low", "Mcov.up") 
  
  s1pM.covtab$row <- as.character(s1pM.covtab$row) # convert into character to replace
  s1pM.covtab$col <- as.character(s1pM.covtab$col)
  s1pM.covtab[s1pM.covtab == "V1_out"] <- "f1@A"; s1pM.covtab[s1pM.covtab == "V1_in"]  <- "f1@P"
  s1pM.covtab[s1pM.covtab == "V2_out"] <- "f2@A"; s1pM.covtab[s1pM.covtab == "V2_in"]  <- "f2@P"
  s1pM.covtab[s1pM.covtab == "V3_out"] <- "f3@A"; s1pM.covtab[s1pM.covtab == "V3_in"]  <- "f3@P"
  s1pM.covtab$par_names <- paste0(s1pM.covtab$row, "~~", s1pM.covtab$col)
  s1pM.covtab <- subset(s1pM.covtab, select = -c(row, col))
  
  pSigma <- merge(pSigma, s1pM.covtab, by = "par_names") # merged dataset with EAPs and MAPs
  
  # dyad-level covariances
  dcovElements <- grep(pattern = "dSigma", rownames(s1long))
  dSigma <- s1long[dcovElements, ] # EAPs
  dSigma$row <- rep(dVarNames, each = 6)
  dSigma$col <- rep(dVarNames, times = 6)
  dSigma$par_names <- paste0(dSigma$row, "~~", dSigma$col)
  dSigma <- subset(dSigma, select = -c(row, col))
  dSigmaRows <- paste0("dSigma[", d.idx, "]")
  dSigma <- dSigma[rownames(dSigma) %in% dSigmaRows, ]
  colnames(dSigma) <- c("iter", "Ecov", "Ecov.MCE", "Ecov.SE", 
                        "Ecov.low", "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", 
                        "par_names")
  dSigma$level <- "dyad"
  
  s1dM.covmat <- summary(s1ests, component = "dyad", interval = "hdi",
                         posterior.est = "mode", srm.param = "cov") # MAPs
  s1dM.covEsts <- as.data.frame(as.table(s1dM.covmat$dyad$cov$mode))
  s1dM.covlow  <- as.data.frame(as.table(s1dM.covmat$group$cov$hdi$lower))
  s1dM.covup   <- as.data.frame(as.table(s1dM.covmat$group$cov$hdi$upper))
  s1dM.covlist <- list(s1dM.covEsts, s1dM.covlow, s1dM.covup)
  s1dM.covtab  <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), s1dM.covlist)
  colnames(s1dM.covtab) <-c("row", "col", "Mcov", "Mcov.low", "Mcov.up") 
  
  s1dM.covtab$row <- as.character(s1dM.covtab$row) # convert into character to replace
  s1dM.covtab$col <- as.character(s1dM.covtab$col)
  s1dM.covtab[s1dM.covtab == "V1_ij"] <- "f1@AP"; s1dM.covtab[s1dM.covtab == "V1_ji"]  <- "f1@PA"
  s1dM.covtab[s1dM.covtab == "V2_ij"] <- "f2@AP"; s1dM.covtab[s1dM.covtab == "V2_ji"]  <- "f2@PA"
  s1dM.covtab[s1dM.covtab == "V3_ij"] <- "f3@AP"; s1dM.covtab[s1dM.covtab == "V3_ji"]  <- "f3@PA"
  s1dM.covtab$par_names <- paste0(s1dM.covtab$row, "~~", s1dM.covtab$col)
  s1dM.covtab <- subset(s1dM.covtab, select = -c(row, col))
  
  dSigma <- merge(dSigma, s1dM.covtab, by = "par_names") # merged dataset with EAPs and MAPs
  
  # case-level correlations
  pcorElements <- grep(pattern = "Rp", rownames(s1long))
  pR <- s1long[pcorElements, ] #EAPs
  pR$row <- rep(pVarNames, each = 6)
  pR$col <- rep(pVarNames, times = 6)
  pR$par_names <- paste0(pR$row, "~~", pR$col)
  pR <- subset(pR, select = -c(row, col))
  pRRows <- paste0("Rp[", p.idx, "]")
  pR <- pR[rownames(pR) %in% pRRows, ]
  colnames(pR) <- c("iter", "Ecor", "Ecor.MCE", "Ecor.SE", "Ecor.low", "Ecor.up", 
                    "Ecor.n_eff", "Ecor.Rhat", "par_names")
  pR$level <- "case"
  
  s1pM.cormat <- summary(s1ests, component = "case", interval = "hdi",
                         posterior.est = "mode", srm.param = "cor") #MAPs
  s1pM.corEsts <- as.data.frame(as.table(s1pM.cormat$case$cor$mode))
  s1pM.corlow  <- as.data.frame(as.table(s1pM.cormat$group$cor$hdi$lower))
  s1pM.corup   <- as.data.frame(as.table(s1pM.cormat$group$cor$hdi$upper))
  s1pM.corlist <- list(s1pM.corEsts, s1pM.corlow, s1pM.corup)
  s1pM.cortab  <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), s1pM.corlist)
  colnames(s1pM.cortab) <- c("row", "col", "Mcor", "Mcor.low", "Mcor.up")
  
  s1pM.cortab$row <- as.character(s1pM.cortab$row) # convert into character to replace
  s1pM.cortab$col <- as.character(s1pM.cortab$col)
  s1pM.cortab[s1pM.cortab == "V1_out"] <- "f1@A"; s1pM.cortab[s1pM.cortab == "V1_in"]  <- "f1@P"
  s1pM.cortab[s1pM.cortab == "V2_out"] <- "f2@A"; s1pM.cortab[s1pM.cortab == "V2_in"]  <- "f2@P"
  s1pM.cortab[s1pM.cortab == "V3_out"] <- "f3@A"; s1pM.cortab[s1pM.cortab == "V3_in"]  <- "f3@P"
  s1pM.cortab$par_names <- paste0(s1pM.cortab$row, "~~", s1pM.cortab$col)
  s1pM.cortab <- subset(s1pM.cortab, select = -c(row, col))
  
  pR <- merge(pR, s1pM.cortab, by = "par_names")
  
  # dyad-level correlations
  dcorElements <- grep(pattern = "Rd2", rownames(s1long))
  dR <- s1long[dcorElements, ] # EAPs
  dR$row <- rep(dVarNames, each = 6)
  dR$col <- rep(dVarNames, times = 6)
  dR$par_names <- paste0(dR$row, "~~", dR$col)
  dR <- subset(dR, select = -c(row, col))
  dRRows <- paste0("Rd2[", d.idx, "]")
  dR <- dR[rownames(dR) %in% dRRows, ]
  colnames(dR) <- c("iter", "Ecor", "Ecor.MCE", "Ecor.SE", "Ecor.low", "Ecor.up", 
                    "Ecor.n_eff", "Ecor.Rhat", "par_names")
  dR$level <- "dyad"
  
  s1dM.cormat <- summary(s1ests, component = "dyad", interval = "hdi",
                         posterior.est = "mode", srm.param = "cor") # MAPs
  ## dyadic reciprocity == DIAGONAL
  ## intra == BELOW
  ## inter == ABOVE
  s1dM.corEsts <- as.data.frame(as.table(s1dM.cormat$dyad$cor$mode)) 
  s1dM.corlow  <- as.data.frame(as.table(s1dM.cormat$group$cor$hdi$lower))
  s1dM.corup   <- as.data.frame(as.table(s1dM.cormat$group$cor$hdi$upper))
  s1dM.corlist <- list(s1dM.corEsts, s1dM.corlow, s1dM.corup)
  s1dM.cortab  <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), s1dM.corlist)
  
  # changing names and adding subscripts (to match ogsrm)
  s1dM.cortab$Var1 <- c("f1@AP", "f1@AP", "f1@AP",
                        "f1@AP", "f2@AP", "f2@AP", 
                        "f1@AP", "f2@AP", "f3@AP")
  s1dM.cortab$Var2 <- c("f1@PA", "f2@PA", "f3@PA", 
                        "f2@AP", "f2@PA", "f3@PA", 
                        "f3@AP", "f3@AP", "f3@PA")
  colnames(s1dM.cortab) <- c("row", "col", "Mcor", "Mcor.low", "Mcor.up")
  
  s1dM.cortab$par_names <- paste0(s1dM.cortab$row, "~~", s1dM.cortab$col)
  s1dM.cortab <- subset(s1dM.cortab, select = -c(row, col))
  
  dR <- merge(dR, s1dM.cortab, by = "par_names", all.x = TRUE)
  
  # case-level SDs
  pSDElements <- grep(pattern = "S_p", rownames(s1long))
  pSD <- s1long[pSDElements, ] # EAPs
  pSD$row <- pVarNames
  pSD$col <- pVarNames
  pSD$par_names <- paste0(pSD$row, "~~", pSD$col)
  pSD <- subset(pSD, select = -c(row, col))
  colnames(pSD) <- c("iter", "Esd", "Esd.MCE", "Esd.SE", "Esd.low", "Esd.up", 
                     "Esd.n_eff", "Esd.Rhat", "par_names")
  pSD$level <- "case"
  
  s1pM.SD <- as.data.frame(as.table(summary(s1ests, component = "case", 
                                            posterior.est = "mode")$case$sd$mode)) #MAPs
  s1pM.SD$Var1 <- as.character(s1pM.SD$Var1)
  s1pM.SD[s1pM.SD == "V1_out"] <- "f1@A"; s1pM.SD[s1pM.SD == "V1_in"]  <- "f1@P"
  s1pM.SD[s1pM.SD == "V2_out"] <- "f2@A"; s1pM.SD[s1pM.SD == "V2_in"]  <- "f2@P"
  s1pM.SD[s1pM.SD == "V3_out"] <- "f3@A"; s1pM.SD[s1pM.SD == "V3_in"]  <- "f3@P"
  s1pM.SD$par_names <- paste0(s1pM.SD$Var1, "~~", s1pM.SD$Var1)
  s1pM.SD <- subset(s1pM.SD, select = c("par_names", "Freq"))
  names(s1pM.SD)[names(s1pM.SD) == "Freq"] <- "Msd"
  s1pM.SD$Msd.low <- summary(s1ests, component = "case", 
                             posterior.est = "mode")$case$sd$central["lower",] # lower hdi limit
  s1pM.SD$Msd.up <- summary(s1ests, component = "case", 
                             posterior.est = "mode")$case$sd$central["upper",] # upper hdi limit
  
  pSD <- merge(pSD, s1pM.SD, by = "par_names")
  
  #dyad-level SDs
  dSDElements <- grep(pattern = "s_rr", rownames(s1long))
  dSD <- s1long[dSDElements, ]
  dSD$row <- dVarNames[c(1, 3, 5)]
  dSD$col <- dVarNames[c(1, 3, 5)]
  dSD$par_names <- paste0(dSD$row, "~~", dSD$col)
  dSD <- subset(dSD, select = -c(row, col))
  colnames(dSD) <- c("iter", "Esd", "Esd.MCE", "Esd.SE", "Esd.low", "Esd.up", 
                     "Esd.n_eff", "Esd.Rhat", "par_names")
  dSD$level <- "dyad"
  
  s1dM.SD <- as.data.frame(as.table(summary(s1ests, component = "dyad", 
                                            posterior.est = "mode")$dyad$sd$mode))
  s1dM.SD$Var1 <- as.character(s1dM.SD$Var1)
  s1dM.SD[s1dM.SD == "V1"] <- "f1@AP"
  s1dM.SD[s1dM.SD == "V2"]  <- "f2@AP"
  s1dM.SD[s1dM.SD == "V3"] <- "f3@AP"
  s1dM.SD$par_names <- paste0(s1dM.SD$Var1, "~~", s1dM.SD$Var1)
  s1dM.SD <- subset(s1dM.SD, select = c("par_names", "Freq"))
  names(s1dM.SD)[names(s1dM.SD) == "Freq"] <- "Msd"
  s1dM.SD$Msd.low <- summary(s1ests, component = "dyad", 
                             posterior.est = "mode")$dyad$sd$central["lower",] # lower hdi limit
  s1dM.SD$Msd.up <- summary(s1ests, component = "dyad", 
                             posterior.est = "mode")$dyad$sd$central["upper",] # upper hdi limit
  
  dSD <- merge(dSD, s1dM.SD, by = "par_names")
  
  #end: saving results----
  
  # begin: compiling final results ----
  
  # rbind() the person and dyad levels
  Sigma <- rbind(pSigma, dSigma)
  R <- rbind(pR, dR)
  SD <- rbind(pSD, dSD)
  
  # remove redundant rows in R (rows with correlations of 1)
  R <- R[!R$Ecor == 1,]
  rownames(R) <- NULL
  
  # merge all values with their population values
  popVals <- getSigma(return_mats = FALSE)
  Sigma <- merge(Sigma, popVals$pop.cov, by = "par_names") # (co)variances
  R <- merge(R, popVals$pop.cor, by = "par_names") # correlations
  SD <- merge(SD, popVals$pop.SD, by = "par_names") # standard deviations
  
  # add MCSampID, n, and G columns to final dataframes
  Sigma$MCSampID <- MCSampID; Sigma$n <- n; Sigma$G <- G; Sigma$condition <- paste0(n, "-", G)
  R$MCSampID <- MCSampID; R$n <- n; R$G <- G; R$condition <- paste0(n, "-", G)
  SD$MCSampID <- MCSampID; SD$n <- n; SD$G <- G; SD$condition <- paste0(n, "-", G)
  
  
  # add analType columns to the final dataframes
  Sigma$analType <- paste0("MCMC-", priorType, "-", ifelse(!missing(precision), precision, "SE"))
  R$analType <- paste0("MCMC-", priorType, "-", ifelse(!missing(precision), precision, "SE"))
  SD$analType <- paste0("MCMC-", priorType, "-", ifelse(!missing(precision), precision, "SE"))
  
  # add prior1 and prior2 columns to R and SD dataframes
  R$prior1 <- NA; R$prior2 <- NA; SD$prior1 <- NA; SD$prior2 <- NA
  
  for(i in 1:length(rr.vars)) {
    ### outgoing SD prior parameters
    SD$prior1[SD$par_names == paste0("f", i, "@A~~f", i, "@A")] <- 
      s1_priors$rr_out_t[paste0("V",i), "m"]
    SD$prior2[SD$par_names == paste0("f", i, "@A~~f", i, "@A")] <- 
      s1_priors$rr_out_t[paste0("V",i), "sd"]
    ### incoming SD prior parameters
    SD$prior1[SD$par_names == paste0("f", i, "@P~~f", i, "@P")] <- 
      s1_priors$rr_in_t[paste0("V",i), "m"]
    SD$prior2[SD$par_names == paste0("f", i, "@P~~f", i, "@P")] <- 
      s1_priors$rr_in_t[paste0("V",i), "sd"]
    ## dyad-level SD prior parameters
    SD$prior1[SD$par_names == paste0("f", i, "@AP~~f", i, "@AP")] <- 
      s1_priors$rr_rel_t[paste0("V",i), "m"]
    SD$prior2[SD$par_names == paste0("f", i, "@AP~~f", i, "@AP")] <- 
      s1_priors$rr_rel_t[paste0("V",i), "sd"]
    ## generalised cov
    R$prior1[R$par_names == paste0("f", i, "@A~~f", i, "@P")] <- 
      s1_priors$case_beta[paste0("V", i, "_in"), paste0("V", i, "_out")]
    R$prior2[R$par_names == paste0("f", i, "@A~~f", i, "@P")] <- 
      s1_priors$case_beta[paste0("V", i, "_out"), paste0("V", i, "_in")]
    
    ## dyad level correlation prior parameters
    ### dyadic--diag
    R$prior1[R$par_names == paste0("f", i, "@AP~~f", i, "@PA")] <- 
      s1_priors$rr_beta_a[paste0("V", i), paste0("V", i)]
    R$prior2[R$par_names == paste0("f", i, "@AP~~f", i, "@PA")] <- 
      s1_priors$rr_beta_b[paste0("V", i), paste0("V", i)]
    
    if (i > 1L) for (j in 1:(i-1)) {
      ### intra--below
      R$prior1[R$par_names == paste0("f", j, "@AP~~f", i, "@AP")] <- 
        s1_priors$rr_beta_a[paste0("V", i), paste0("V", j)]
      R$prior2[R$par_names == paste0("f", j, "@AP~~f", i, "@AP")] <- 
        s1_priors$rr_beta_b[paste0("V", i), paste0("V", j)]
      ### inter--above
      R$prior1[R$par_names == paste0("f", j, "@AP~~f", i, "@PA")] <- 
        s1_priors$rr_beta_a[paste0("V", j), paste0("V", i)]
      R$prior2[R$par_names == paste0("f", j, "@AP~~f", i, "@PA")] <- 
        s1_priors$rr_beta_b[paste0("V", j), paste0("V", i)]
      
      ## case-level correlation prior parameters
      ### AA 
      R$prior1[R$par_names == paste0("f", j, "@A~~f", i, "@A")] <- 
        s1_priors$case_beta[paste0("V", i, "_out"), paste0("V", j, "_out")]
      R$prior2[R$par_names == paste0("f", j, "@A~~f", i, "@A")] <- 
        s1_priors$case_beta[paste0("V", j, "_out"), paste0("V", i, "_out")]
      ### AP
      R$prior1[R$par_names == paste0("f", j, "@A~~f", i, "@P")] <- 
        s1_priors$case_beta[paste0("V", i, "_in"), paste0("V", j, "_out")]
      R$prior2[R$par_names == paste0("f", j, "@A~~f", i, "@P")] <- 
        s1_priors$case_beta[paste0("V", j, "_out"), paste0("V", i, "_in")]
      ### PA
      R$prior1[R$par_names == paste0("f", j, "@P~~f", i, "@A")] <- 
        s1_priors$case_beta[paste0("V", i, "_out"), paste0("V", j, "_in")]
      R$prior2[R$par_names == paste0("f", j, "@P~~f", i, "@A")] <- 
        s1_priors$case_beta[paste0("V", j, "_in"), paste0("V", i, "_out")]
      ### PP
      R$prior1[R$par_names == paste0("f", j, "@P~~f", i, "@P")] <- 
        s1_priors$case_beta[paste0("V", i, "_in"), paste0("V", j, "_in")]
      R$prior2[R$par_names == paste0("f", j, "@P~~f", i, "@P")] <- 
        s1_priors$case_beta[paste0("V", j, "_in"), paste0("V", i, "_in")]
    }
  }
  
  # add run time
  Sigma$RunTime <- difftime(t1, t0, units = "mins")
  R$RunTime <- difftime(t1, t0, units = "mins")
  SD$RunTime <- difftime(t1, t0, units = "mins")
  
  # reorder columns and add redundant columns (those in og) to make postprocessing easier
  Sigma <- Sigma[, c("par_names", "MCSampID", "n", "G", "condition", "level",
                     "analType", "iter", "RunTime", "pop.cov", "Ecov", "Ecov.MCE",
                     "Ecov.SE", "Ecov.low", "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", 
                     "Mcov", "Mcov.low", "Mcov.up")]
  R <- R[, c("par_names", "MCSampID", "n", "G", "condition", "level",
                     "analType", "iter", "RunTime", "pop.cor", "Ecor", "Ecor.MCE",
                     "Ecor.SE", "Ecor.low", "Ecor.up", "Ecor.n_eff", "Ecor.Rhat", 
                     "Mcor", "Mcor.low", "Mcor.up", "prior1", "prior2")]
  SD <- SD[, c("par_names", "MCSampID", "n", "G", "condition", "level",
             "analType", "iter", "RunTime", "pop.SD", "Esd", "Esd.MCE",
             "Esd.SE", "Esd.low", "Esd.up", "Esd.n_eff", "Esd.Rhat", 
             "Msd", "Msd.low", "Msd.up", "prior1", "prior2")]
  
  Sigma$ogcov <- NA; Sigma$ogcov.SE <- NA; Sigma$ogcov.low <- NA; Sigma$ogcov.up <- NA
  R$ogcor <- NA; R$ogcor.SE <- NA; R$ogcor.low <- NA; R$ogcor.up <- NA
  SD$ogsd <- NA; SD$ogsd.SE <- NA; SD$ogsd.low <- NA; SD$ogsd.up <- NA
  
  # combine covariances and correlations + SDs as list
  out <- list(cov = Sigma, cor = R, SD = SD, 
              mPSRF = c(MCSampID = MCSampID, condition = paste0(n, "-", G), 
                        analType = paste0("MCMC-", priorType, "-", 
                                          ifelse(!missing(precision), precision, "SE")),
                        mPSRF = mPSRF))
  
  # end: compiling final results ----
  
  if (savefile) saveRDS(out, 
                        file = paste0("ID", MCSampID, ".nG", G, ".n", n, "-", 
                                      priorType, "-", ifelse(!missing(precision), precision, "SE"), ".rds"))
  
  return(out)

}

### README you had the same environment() issue as above within this function, which 
### you fixed by creating a new environment to store the data and then get()ing that 
### data from the new environment each time

# s1sat(MCSampID = 1, n = 5, G = 3, rr.vars = c("V1", "V2", "V3"), IDout = "Actor",
#       IDin = "Partner", IDgroup = "Group", priorType = "default", precision = 0.1,
#       iter = 100)
# s1sat(MCSampID = 1, n = 5, G = 3, rr.vars = c("V1", "V2", "V3"), IDout = "Actor",
#       IDin = "Partner", IDgroup = "Group", priorType = "thoughtful", targetCorr = 0.3,
#       precision = 0.1, iter = 100)
# s1sat(MCSampID = 1, n = 5, G = 3, rr.vars = c("V1", "V2", "V3"), IDout = "Actor",
#       IDin = "Partner", IDgroup = "Group", priorType = "prophetic",
#       precision = 0.1, iter = 100)
# s1sat(MCSampID = 1, n = 5, G = 3, rr.vars = c("V1", "V2", "V3"), IDout = "Actor",
#       IDin = "Partner", IDgroup = "Group", priorType = "ANOVA",
#       precision = 0.1, iter = 100)
# s1sat(MCSampID = 1, n = 5, G = 3, rr.vars = c("V1", "V2", "V3"), IDout = "Actor",
#       IDin = "Partner", IDgroup = "Group", priorType = "FIML",
#       precision = 0.1, multiMLE = F, iter = 100)

#----

# function 11: fit saturated model in `srm`----

ogsat <- function(MCSampID, n, G, savefile = FALSE) {
  
  library(srm)
  library(car)
  
  dat <- genGroups(MCSampID = MCSampID, n = n, G = G)
  
  t0 <- Sys.time()
  
  # specify satmod 
  satsyn <- ' %Person
  f1@A =~ 1*V1@A
  f2@A =~ 1*V2@A
  f3@A =~ 1*V3@A
  
  f1@P =~ 1*V1@P
  f2@P =~ 1*V2@P
  f3@P =~ 1*V3@P
  
  V1@A ~~ 0*V1@A + 0*V1@P
  V1@P ~~ 0*V1@P
  V2@A ~~ 0*V2@A + 0*V2@P
  V2@P ~~ 0*V2@P
  V3@A ~~ 0*V3@A + 0*V3@P
  V3@P ~~ 0*V3@P
  
  f1@A ~~ f1@A + f2@A + f3@A + f1@P + f2@P + f3@P
  f2@A ~~        f2@A + f3@A + f1@P + f2@P + f3@P
  f3@A ~~               f3@A + f1@P + f2@P + f3@P
  f1@P ~~                      f1@P + f2@P + f3@P
  f2@P ~~                             f2@P + f3@P
  f3@P ~~                                    f3@P
  
  
  %Dyad
  f1@AP =~ 1*V1@AP
  f2@AP =~ 1*V2@AP
  f3@AP =~ 1*V3@AP
  
  # f1@PA =~ 1*V1@PA
  # f2@PA =~ 1*V2@PA
  # f3@PA =~ 1*V3@PA
  
  V1@AP ~~ 0*V1@AP + 0*V1@PA
  V1@PA ~~ 0*V1@PA
  V2@AP ~~ 0*V2@AP + 0*V2@PA
  V2@PA ~~ 0*V2@PA
  V3@AP ~~ 0*V3@AP + 0*V3@PA
  V3@PA ~~ 0*V3@PA
  
  f1@AP ~~ relvar1*f1@AP + intra12*f2@AP + intra13*f3@AP + dyad11*f1@PA  + inter12*f2@PA + inter13*f3@PA
  f2@AP ~~                 relvar2*f2@AP + intra23*f3@AP + inter12*f1@PA +  dyad22*f2@PA + inter23*f3@PA
  f3@AP ~~                                 relvar3*f3@AP + inter13*f1@PA + inter23*f2@PA +  dyad33*f3@PA
  f1@PA ~~                                                 relvar1*f1@PA + intra12*f2@PA + intra13*f3@PA
  f2@PA ~~                                                                 relvar2*f2@PA + intra23*f3@PA
  f3@PA ~~                                                                                 relvar3*f3@PA
  '
  
  satfit <- srm(satsyn, data = dat, rrgroup_name = "Group",
                person_names = c("Actor", "Partner"), 
                fixed.groups = FALSE,
                verbose = FALSE)
  
  t1 <- Sys.time()
  
  if (!satfit$res_opt$converged) return(NULL)
  
  ## cov results
  covVals <- satfit$parm.table
  
  # specifying the levels per the notation we want to use
  covVals$level[covVals$level == "U"] <- "case"
  covVals$level[covVals$level == "D"] <- "dyad"
  
  # indices of parameters we want to extract
  p.idx <- outer(1:6, 1:6, FUN = paste, sep = ",")[upper.tri(diag(6), diag = TRUE)]
  d.idx <- c(paste0(rep(1, each = 6), ",", 1:6), # relvar1; dyadcov11; intra/inter12; intra/inter13
             paste0(rep(3, each = 4), ",", 3:6), #relvar2; dyadcov22; intra/inter23
             paste0(rep(5, each = 2), ",", 5:6)) # relvar3; dyadcov33 
  
  # split the dataframe by level
  p.covVals <- covVals[covVals$level == "case" & covVals$mat == "PHI_U", ]
  d.covVals <- covVals[covVals$level == "dyad" & covVals$mat == "PHI_D", ]
  
  # create an id variable in each dataframe
  p.covVals$id <- paste0(p.covVals$row, ',', p.covVals$col)
  d.covVals$id <- paste0(d.covVals$row, ',', d.covVals$col)
  
  # extract elements based on the id variable
  p.covVals <- p.covVals[p.covVals$id %in% p.idx, ]
  d.covVals <- d.covVals[d.covVals$id %in% d.idx, ]
  
  # replace 3 par_names in p.covVals to make these consistent with lavsat() output
  p.covVals$par_names[p.covVals$par_names == "f2@A~~f1@P"] <- "f1@P~~f2@A"
  p.covVals$par_names[p.covVals$par_names == "f3@A~~f1@P"] <- "f1@P~~f3@A"
  p.covVals$par_names[p.covVals$par_names == "f3@A~~f2@P"] <- "f2@P~~f3@A"
  
  # rbind() the datasets back together
  covVals <- rbind(p.covVals, d.covVals)
  
  # subsetting only the relevant columns
  covVals <- subset(covVals, select = c("par_names", "level", "est", "se"))
  colnames(covVals) <- c("par_names", "level", "ogcov", "ogcov.SE")
  
  # calculating the upper and lower limits of the confidence interval
  covVals$ogcov.low <- covVals$ogcov - 1.96*covVals$ogcov.SE
  covVals$ogcov.up <- covVals$ogcov + 1.96*covVals$ogcov.SE
  
  # creating a final cov dataframe
  popVals_cov <- getSigma(return_mats = FALSE)$pop.cov
  covResult <- merge(covVals, popVals_cov, by = "par_names", sort = FALSE)
  rownames(covResult) <- NULL
  
  ## add columns pertaining to the simulation
  covResult$MCSampID <- MCSampID; covResult$n <- n; covResult$G <- G
  covResult$iter <- satfit$res_opt$iter # save number of MLE iterations
  covResult$RunTime <- difftime(t1, t0, units = "mins")
  covResult$condition <- paste0(n, "-", G)
  covResult$analType <- "1S-FIML"
  
  ## cor results
  p.names <- paste0("f", rep(1:3, times = 2), "@", rep(c("A", "P"), each = 3))
  # d.names <- paste0("f", rep(1:3, times = 2), "@", rep(c("AP", "PA"), each = 3))
  var.names <- c("f1", "f2", "f3")
  
  ## correlation and SD results
  p.SD <- NULL; d.SD <- NULL; p.cor <- NULL; d.cor <- NULL
  
  for (rr in 1:length(p.names)) { # case level
    p.SDval <- deltaMethod(satfit, g. = paste0("sqrt(`", p.names[rr], "~~", p.names[rr], "`)"))
    p.SDval$par_names <- paste0(p.names[rr], "~~", p.names[rr])
    rownames(p.SDval) <- NULL
    if(is.null(p.SD)) {
      p.SD <- p.SDval
    } else {
      p.SD <- rbind(p.SD, p.SDval)
    }
    if (rr > 1L) for (kk in 1:(rr-1)) {
      p.corVal <- deltaMethod(satfit, g. = paste0("`", p.names[kk], "~~", p.names[rr],"`/(sqrt(`",
                                                  p.names[rr], "~~", p.names[rr], "`*`", 
                                                  p.names[kk], "~~", p.names[kk], "`))"))
      p.corVal$par_names <- paste0(p.names[kk], "~~", p.names[rr])
      rownames(p.corVal) <- NULL
      if(is.null(p.cor)) {
        p.cor <- p.corVal
      } else {
        p.cor <- rbind(p.cor, p.corVal)
      }
    } 
  }
  p.SD <- as.data.frame(p.SD)
  p.cor <- as.data.frame(p.cor)
  p.SD$level <- "case"; p.cor$level <- "case"
  
  # replace 3 par_names in p.cor to make these consistent with lavsat() output
  p.cor$par_names[p.cor$par_names == "f2@A~~f1@P"] <- "f1@P~~f2@A"
  p.cor$par_names[p.cor$par_names == "f3@A~~f1@P"] <- "f1@P~~f3@A"
  p.cor$par_names[p.cor$par_names == "f3@A~~f2@P"] <- "f2@P~~f3@A"
  
  for (ll in 1:length(var.names)) {
    ## SDs
    d.SDval <- deltaMethod(satfit, g. = paste0("sqrt(`", var.names[ll], 
                                               "@AP~~", var.names[ll], "@AP`)"))
    d.SDval$par_names <- paste0(var.names[ll], "@AP~~", 
                                var.names[ll], "@AP")
    rownames(d.SDval) <- NULL
    if(is.null(d.SD)) {
      d.SD <- d.SDval
    } else {
      d.SD <- rbind(d.SD, d.SDval)
    }
    ## dyadic cor
    d.dyadicval <- deltaMethod(satfit, g. = paste0("`", var.names[ll], "@AP~~", 
                                                   var.names[ll],"@PA`/`", var.names[ll], 
                                                   "@AP~~", var.names[ll],"@AP`"))
    d.dyadicval$par_names <- paste0(var.names[ll], "@AP~~", var.names[ll], "@PA")
    rownames(d.dyadicval) <- NULL
    if(is.null(d.cor)) {
      d.cor <- d.dyadicval
    } else {
      d.cor <- rbind(d.cor, d.dyadicval)
    }
    if (ll > 1L) for (mm in 1:(ll-1)) {
      ## intra cor
      d.intraval <- deltaMethod(satfit, g. = paste0("`", var.names[mm], "@AP~~", 
                                                    var.names[ll],"@AP`/(sqrt(`", var.names[ll], 
                                                    "@AP~~", var.names[ll], "@AP`*`",
                                                    var.names[mm], "@AP~~", var.names[mm], "@AP`))"))
      d.intraval$par_names <- paste0(var.names[mm], "@AP~~", var.names[ll], "@AP")
      rownames(d.intraval) <- NULL
      if(is.null(d.cor)) {
        d.cor <- d.intraval
      } else {
        d.cor <- rbind(d.cor, d.intraval)
      }
      ## inter cor
      d.interval <- deltaMethod(satfit, g. = paste0("`", var.names[mm], "@AP~~", 
                                                    var.names[ll],"@PA`/(sqrt(`", var.names[ll], 
                                                    "@AP~~", var.names[ll], "@AP`*`",
                                                    var.names[mm], "@AP~~", var.names[mm], "@AP`))"))
      d.interval$par_names <- paste0(var.names[mm], "@AP~~", var.names[ll], "@PA")
      rownames(d.interval) <- NULL
      if(is.null(d.cor)) {
        d.cor <- d.interval
      } else {
        d.cor <- rbind(d.cor, d.interval)
      }
    }
  }
  d.SD <- as.data.frame(d.SD)
  d.cor <- as.data.frame(d.cor)
  d.SD$level <- "dyad"; d.cor$level <- "dyad"
  
  SDResult <- rbind(p.SD, d.SD)
  SDResult <- SDResult[, c("par_names", "level", "Estimate", "SE", "2.5 %", "97.5 %")]
  names(SDResult) <- c("par_names", "level", "ogsd", "ogsd.SE", "ogsd.low", "ogsd.up")
  SDResult$MCSampID <- MCSampID; SDResult$n <- n; SDResult$G <- G; SDResult$condition <- paste0(n, "-", G)
  SDResult$iter <- satfit$res_opt$iter
  SDResult$RunTime <- difftime(t1, t0, units = "mins")
  popVals_SD <- getSigma(return_mats = FALSE)$pop.SD
  SDResult <- merge(SDResult, popVals_SD, by = "par_names", sort = FALSE)
  SDResult$analType <- "1S-FIML"
  
  corResult <- rbind(p.cor, d.cor)
  corResult <- corResult[, c("par_names", "level", "Estimate", "SE", "2.5 %", "97.5 %")]
  names(corResult) <- c("par_names", "level", "ogcor", "ogcor.SE", "ogcor.low", "ogcor.up")
  corResult$MCSampID <- MCSampID; corResult$n <- n; corResult$G <- G; corResult$condition <- paste0(n, "-", G)
  corResult$iter <- satfit$res_opt$iter
  corResult$RunTime <- difftime(t1, t0, units = "mins")
  popVals_cor <- getSigma(return_mats = FALSE)$pop.cor
  corResult <- merge(corResult, popVals_cor, by = "par_names", sort = FALSE)
  corResult$analType <- "1S-FIML"
  
  # add redundant columns (those in s1sat) and reorder columns to make postporcessing easier
  covResult$Ecov <- NA; covResult$Ecov.MCE <- NA; covResult$Ecov.SE <- NA
  covResult$Ecov.low <- NA; covResult$Ecov.up <- NA; covResult$Ecov.n_eff <- NA
  covResult$Ecov.Rhat <- NA; covResult$Mcov <- NA; covResult$Mcov.low <- NA;
  covResult$Mcov.up <- NA
  corResult$Ecor <- NA; corResult$Ecor.MCE <- NA; corResult$Ecor.SE <- NA
  corResult$Ecor.low <- NA; corResult$Ecor.up <- NA; corResult$Ecor.n_eff <- NA
  corResult$Ecor.Rhat <- NA; corResult$Mcor <- NA; corResult$Mcor.low <- NA;
  corResult$Mcor.up <- NA; corResult$prior1 <- NA; corResult$prior2 <- NA
  SDResult$Esd <- NA; SDResult$Esd.MCE <- NA; SDResult$Esd.SE <- NA
  SDResult$Esd.low <- NA; SDResult$Esd.up <- NA; SDResult$Esd.n_eff <- NA
  SDResult$Esd.Rhat <- NA; SDResult$Msd <- NA; SDResult$Msd.low <- NA;
  SDResult$Msd.up <- NA; SDResult$prior1 <- NA; SDResult$prior2 <- NA
  
  covResult <- covResult[, c("par_names", "MCSampID", "n", "G", "condition", "level",
                             "analType", "iter", "RunTime", "pop.cov", "Ecov", "Ecov.MCE",
                             "Ecov.SE", "Ecov.low", "Ecov.up", "Ecov.n_eff", "Ecov.Rhat", 
                             "Mcov", "Mcov.low", "Mcov.up", "ogcov", "ogcov.SE",
                             "ogcov.low", "ogcov.up")]
  corResult <- corResult[, c("par_names", "MCSampID", "n", "G", "condition", "level",
                             "analType", "iter", "RunTime", "pop.cor", "Ecor", "Ecor.MCE",
                             "Ecor.SE", "Ecor.low", "Ecor.up", "Ecor.n_eff", "Ecor.Rhat", 
                             "Mcor", "Mcor.low", "Mcor.up", "prior1", "prior2",
                             "ogcor", "ogcor.SE", "ogcor.low", "ogcor.up")]
  SDResult <- SDResult[, c("par_names", "MCSampID", "n", "G", "condition", "level",
                           "analType", "iter", "RunTime", "pop.SD", "Esd", "Esd.MCE",
                           "Esd.SE", "Esd.low", "Esd.up", "Esd.n_eff", "Esd.Rhat", 
                           "Msd", "Msd.low", "Msd.up", "prior1", "prior2",
                           "ogsd", "ogsd.SE", "ogsd.low", "ogsd.up")]
  
  
  out <- list(cov = covResult, cor = corResult, SD = SDResult)
  
  if (savefile) saveRDS(out, 
                        file = paste0("ID", MCSampID, ".nG", G, ".n", n, 
                                      "-1SFIML.rds"))
  
  return(out) #TODO
}

# ogsat(1, 6, 10)

#----

# function 12: create runsim files----

# README: in this function, 1S-FIML = FIML1S when analType is single-stage FIML

makeRunsim <- function(nSamps, n, G, analType, precision = NULL, sim) {
  runsimfile <- paste0('## Aditi M. Bhangale
## Last updated:', Sys.Date(), 

'\n# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model
# Simulations 1 to 3

# runsim-',analType, ifelse(!is.null(precision), paste0("-", precision), ""), '-', sim,'

source("functions-sim1to3.R")

# specify conditions\n',

analType,'_grid <- expand.grid(MCSampID =1:', nSamps, ', n =', n, ', G = ', G, ',',
                              ifelse(analType != "FIML1S", paste0('priorType = "', analType, '"'), 
                                     paste0('priorType = "NA"')), ',',
                              ifelse(!is.null(precision), paste0('precision = ', precision), 
                                     paste0('priorType = "NA"')),',
                              stringsAsFactors = F)\n',

analType,'_grid$row_num <- 1:nrow(', analType,'_grid)

# prepare parallel processing\n
library(doSNOW)

nClus <- 32
cl <- makeCluster(nClus)
registerDoSNOW(cl)

# run simulation\n',
if (analType == "FIML1S") {
  paste0('ogResult <- foreach(row_num = 1:nrow(',analType,'_grid),
                    .packages = c("mnormt", "parallel", "portableParallelSeeds", 
                                  "srm", "car")) %dopar% {
                                    
                                    out <- try(ogsat(MCSampID = ',analType,'_grid[row_num, ]$MCSampID, 
                                    n = ',analType,'_grid[row_num, ]$n, 
                                                    G = ',analType,'_grid[row_num, ]$G), silent = T)
                                    if(inherits(out, "try-error")) out <- NULL
                                    
                                    return(out)
                                  }
         
   # close cluster\n
   stopCluster(cl)
   
   saveRDS(ogResult, paste0("results_', analType,'-',sim, '-",Sys.Date(),".rds")) #FIXME
         ')
} else {
  paste0('s1Result <- foreach(row_num = 1:nrow(',analType,'_grid),
                    .packages = c("mnormt", "parallel", "portableParallelSeeds",
                                  "TripleR", "srm", "car", "lavaan.srm", "coda")) %dopar% {
                                    
                                    out <- try(s1sat(MCSampID = ',analType,'_grid[row_num, ]$MCSampID, 
                                                     n = ',analType,'_grid[row_num, ]$n, G = ',analType,'_grid[row_num, ]$G,
                                                     priorType = ',analType,'_grid[row_num, ]$priorType, 
                                                     precision = ',analType,'_grid[row_num, ]$precision), silent = T)
                                    if(inherits(out, "try-error")) out <- NULL
                                    
                                    return(out)
                                  }
         
         # close cluster\n
          stopCluster(cl)
         
         saveRDS(s1Result, paste0("results_', analType, '-', precision, '-',
         sim,'-", Sys.Date(),".rds")) #FIXME
         
         ')
}
)
 
  cat(runsimfile, file = paste0("runsim-", analType, ifelse(!is.null(precision), paste0("-", precision), ""), "-", sim, ".R"))
  invisible(NULL)
}

## make runsim files for all conditions (simulations 1 and 2)
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "default",
#            sim = "sim1")
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "prophetic",
#            precision = 0.05, sim = "sim1")
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "prophetic",
#            precision = 0.1, sim = "sim1")
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "prophetic",
#            precision = 0.2, sim = "sim1")
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "FIML1S",
#            sim = "sim1")
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "thoughtful",
#            precision = 0.1, sim = "sim2")
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "ANOVA",
#            precision = 0.1, sim = "sim2")
# makeRunsim(nSamps = 1000, n = "c(6,8,10,20)", G = "c(10, 25)", analType = "FIML",
#            precision = 0.1, sim = "sim2")

#TODO test all the below out with one sample to check if the code works as expected
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "default",
#            sim = "sim1") 
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "prophetic",
#            precision = 0.05, sim = "sim1")
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "prophetic",
#            precision = 0.1, sim = "sim1")
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "prophetic",
#            precision = 0.2, sim = "sim1")
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "FIML1S",
#            sim = "sim1")
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "thoughtful",
#            precision = 0.1, sim = "sim2")
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "ANOVA",
#            precision = 0.1, sim = "sim2")
# makeRunsim(nSamps = 1, n = "c(6,8)", G = 10, analType = "FIML",
#            precision = 0.1, sim = "sim2")

#----

# function 13: create shell files----
makeShell <- function(analType, precision = NULL, sim, wallTime) {
  shell <- paste0('#!/bin/bash

#SBATCH -J ', paste0(analType, ifelse(!is.null(precision), paste0("-", precision), ""), "-", sim),'
#SBATCH -e ', paste0(analType, ifelse(!is.null(precision), paste0("-", precision), ""), "-", sim),'.SERR
#SBATCH -o ', paste0(analType, ifelse(!is.null(precision), paste0("-", precision), ""), "-", sim),'.SOUT
#SBATCH -N 1
#SBATCH --cpus-per-task 32
#SBATCH -t ', wallTime,'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aditibhangale@gmail.com

cd $SLURM_SUBMIT_DIR

module load R/4.2.2
export MKL_NUM_THREADS=1

Rscript --vanilla', paste0("runsim-", analType, ifelse(!is.null(precision), 
                                                       paste0("-", precision), ""), "-", sim, ".R")

)
  cat(shell, file = paste0("shell-", analType, ifelse(!is.null(precision), 
                                                      paste0("-", precision), ""), "-", sim, ".sh"))
  invisible(NULL)
}

makeShell(analType = "FIML1S", sim = "sim1", wallTime = "00:05:00")

#TODO test all the above to check if the code works as expected

#----

#######
