## Aditi M. Bhangale
## Last updated: 25 March 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model
# Simulations 1 and 2

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/sim_code")
# getwd()

# MCSampID = 1; n = 5; G = 3
# rr.data <- genGroups(MCSampID = 1, n = 5, G = 3)
# rr.vars <- c("V1", "V2", "V3")
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"
# precision <- 0.1
# targetCorr <- 0.3

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
  mat_c.cor <- solve(mat_c.SD) %*% SIGMA_c %*% solve(mat_c.SD)
  pop_c.cor <- as.data.frame(as.table(mat_c.cor))
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
  mat_d.cor <- solve(mat_d.SD) %*% SIGMA_d %*% solve(mat_d.SD)
  pop_d.cor <- as.data.frame(as.table(mat_d.cor))
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
  mu <- c(rep(0, 6)) # group-mean centered SRM component variables
  
  ##DATA GENERATION-----
  
  library(rockchalk) # for mvrnorm()
  
  dat_c <- mvrnorm(n = n, mu = mu, Sigma = SIGMA_c)
  dat_d <- mvrnorm(n = (n*(n - 1))/2, mu = mu, Sigma = SIGMA_d) 
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

#FIXME was i right to remove the `data` argument from the `thoughtful_priors()` function?

thoughtful_priors <- function(targetCorr, precision, default_prior) {
  
  priors <- default_prior
  
  # for SDs --- t-priors
  ## m for SDs is already thoughtful/approximate location because lavaan.srm estimates
  ## it based on the data
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

#FIXME was i right to remove the `data` argument from the `prophetic_priors()` function?

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
  Fnames_c <- c("f1@A", "f1@P", "f2@A", "f2@P", "f3@A", "f3@P")
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

#----

# function xxxx: set customised priors for MCMC stage----

# library(lavaan.srm)

set_priors <- function(data, rr.vars, IDout, IDin, IDgroup, priorType, targetCorr,
                       pop_corMat = list(popCorr_c = getSigma()$R_c,
                                         popCorr_d = getSigma()$R_d),
                       pop_SDvec = list(popSD_c = sqrt(diag(getSigma()$SIGMA_c)),
                                        popSD_d = sqrt(diag(getSigma()$SIGMA_d))),
                       precision, multiMLE = FALSE) {
  #TODO load lavaan.srm in the parent function of this one (`s1sat()`, i think?)
  prior_env <- new.env()
  prior_env$default_prior <- srm_priors(rr.data[rr.vars]) # default MCMC priors (diffuse priors)
 
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
   
 }
 return(srmPriors)
}

# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), priorType = "default")
# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"),
#            priorType = "prophetic", precision = 0.1)
# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), priorType = "thoughtful",
#            targetCorr = 0.3, precision = 0.1)
# set_priors(data = rr.data, rr.vars = c("V1", "V2", "V3"), priorType = "ANOVA",
#            IDout = "Actor", IDin = "Partner", IDgroup = "Group", precision = 0.1)

# README https://stackoverflow.com/questions/18338331/nested-function-environment-selection 
# README https://digitheadslabnotebook.blogspot.com/2011/06/environments-in-r.html
# README http://adv-r.had.co.nz/Environments.html#env-answers 

#----
