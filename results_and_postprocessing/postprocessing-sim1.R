## Aditi M. Bhangale
## Last updated: 24 April 2024

# Hyperparameters of empirical Bayes priors for MCMC estimation of the 
# multivariate social relations model

# Postprocessing (outcome variables) and data visualisation for simulation 1

# rm(list = ls())

##README MCSampID 198 in FIML1S 6-10 is not in the output. This is probably because
## that sample did not converge

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage1hyperpars/results_and_postprocessing")
# getwd()

sim1result <- readRDS("sim1result.rds")
for (i in names(sim1result)) {
  sim1result[[i]]$analType <- factor(sim1result[[i]]$analType,
                                     levels = c("MCMC-default-NA", "MCMC-prophetic-0.05",
                                                "MCMC-prophetic-0.1", "MCMC-prophetic-0.2",
                                                "1S-FIML"),
                                     labels = c("Diffuse priors", "Pr-0.05 priors",
                                                "Pr-0.1 priors", "Pr-0.2 priors",
                                                "FIML"))
}

## we only want to compare EAP and ogFIML in this simulation
sim1cov <- sim1result$cov_result[sim1result$cov_result$estType %in% c("EAP", "ogFIML"),]
sim1cor <- sim1result$cor_result[sim1result$cor_result$estType %in% c("EAP", "ogFIML"),]
sim1SD <- sim1result$SD_result[sim1result$cor_result$estType %in% c("EAP", "ogFIML"),]
sim1mPSRF <- sim1result$mPSRF_result

library(ggplot2)

# case level labels, colours, and shapes----
case_pars <- c(# ego var
  "f1@A~~f1@A" = "darkorange1", "f2@A~~f2@A" = "darkorange3", "f3@A~~f3@A" = "darkorange4",
  # alter var
  "f1@P~~f1@P" = "olivedrab1", "f2@P~~f2@P" = "olivedrab3", "f3@P~~f3@P" = "olivedrab4",
  # e-e cov
  "f1@A~~f2@A" = "steelblue1", "f1@A~~f3@A" = "steelblue3", "f2@A~~f3@A" = "steelblue4",
  # a-a cov
  "f1@P~~f2@P" = "hotpink1", "f1@P~~f3@P" = "hotpink3", "f2@P~~f3@P" = "hotpink4",
  #e-a cov
  "f1@A~~f1@P" = "gold1", "f1@A~~f2@P" = "gold3", "f1@A~~f3@P" = "gold4", 
  "f2@A~~f2@P" = "coral1", "f2@A~~f3@P" = "coral3", "f3@A~~f3@P" = "coral4", 
  #a-e cov
  "f1@P~~f2@A" = "slateblue1", "f1@P~~f3@A" = "slateblue3", "f2@P~~f3@A" = "slateblue4")

case_labels <- c(# ego var
  expression(sigma[E[1]]^2), expression(sigma[E[2]]^2), expression(sigma[E[3]]^2),
  # alter var
  expression(sigma[A[1]]^2), expression(sigma[A[2]]^2), expression(sigma[A[3]]^2),
  # e-e cov
  expression(sigma[E[1]][E[2]]), expression(sigma[E[1]][E[3]]), expression(sigma[E[2]][E[3]]),
  # a-a cov
  expression(sigma[A[1]][A[2]]), expression(sigma[A[1]][A[3]]), expression(sigma[A[2]][A[3]]),
  #e-a cov
  expression(sigma[E[1]][A[1]]), expression(sigma[E[1]][A[2]]), expression(sigma[E[1]][A[3]]), 
  expression(sigma[E[2]][A[2]]), expression(sigma[E[2]][A[3]]), expression(sigma[E[3]][A[3]]), 
  expression(sigma[A[1]][E[2]]), expression(sigma[A[1]][E[3]]), expression(sigma[A[2]][E[3]]))

case_shape <- c(# ego var
  "f1@A~~f1@A" = 0, "f2@A~~f2@A" = 0, "f3@A~~f3@A" = 0,
  # alter var
  "f1@P~~f1@P" = 1, "f2@P~~f2@P" = 1, "f3@P~~f3@P" = 1,
  # e-e cov
  "f1@A~~f2@A" = 8, "f1@A~~f3@A" = 8, "f2@A~~f3@A" = 8,
  # a-a cov
  "f1@P~~f2@P" = 12, "f1@P~~f3@P" = 12, "f2@P~~f3@P" = 12,
  #e-a cov
  "f1@A~~f1@P" = 13, "f1@A~~f2@P" = 13, "f1@A~~f3@P" = 13, 
  "f2@A~~f2@P" = 13, "f2@A~~f3@P" = 13, "f3@A~~f3@P" = 13, 
  "f1@P~~f2@A" = 6, "f1@P~~f3@A" = 6, "f2@P~~f3@A" = 6)
#----

# dyad level labels, colours, and shapes----
dyad_pars <- c(# rel var
  "f1@AP~~f1@AP" = "sienna1", "f2@AP~~f2@AP" = "sienna3", "f3@AP~~f3@AP" = "sienna4",
  # dyadic co
  "f1@AP~~f1@PA" = "deeppink1", "f2@AP~~f2@PA" = "deeppink3", "f3@AP~~f3@PA" = "deeppink4",
  # intra cov
  "f1@AP~~f2@AP" = "darkolivegreen1", "f1@AP~~f3@AP" = "darkolivegreen3", "f2@AP~~f3@AP" = "darkolivegreen4",
  # inter cov
  "f1@AP~~f2@PA" = "royalblue1", "f1@AP~~f3@PA" = "royalblue3", "f2@AP~~f3@PA" = "royalblue4")

dyad_labels <- c(# rel var
  expression(sigma[R[1]]^2), expression(sigma[R[2]]^2), expression(sigma[R[3]]^2),
  # dyadic cov
  expression(sigma[R[1]]^dyadic), expression(sigma[R[2]]^dyadic), expression(sigma[R[3]]^dyadic),
  #intra cov
  expression(sigma[R[1]][R[2]]^intra), expression(sigma[R[1]][R[3]]^intra), expression(sigma[R[2]][R[3]]^intra),
  # inter cov
  expression(sigma[R[1]][R[2]]^inter), expression(sigma[R[1]][R[3]]^inter), expression(sigma[R[2]][R[3]]^inter))

dyad_shape <- c(# rel var
  "f1@AP~~f1@AP" = 0, "f2@AP~~f2@AP" = 0, "f3@AP~~f3@AP" = 0,
  # dyadic co
  "f1@AP~~f1@PA" = 1, "f2@AP~~f2@PA" = 1, "f3@AP~~f3@PA" = 1,
  # intra cov
  "f1@AP~~f2@AP" = 2, "f1@AP~~f3@AP" = 2, "f2@AP~~f3@AP" = 2,
  # inter cov
  "f1@AP~~f2@PA" = 5, "f1@AP~~f3@PA" = 5, "f2@AP~~f3@PA" = 5)
#----

#########################################
# POSTPROCSSING FOR CORRELATION ESTIMATES
#########################################

# correlation bias----
sim1cor$cor.bias <- sim1cor$cor.est - sim1cor$pop.cor
sim1cor$cor.RB <- sim1cor$cor.bias / sim1cor$pop.cor

corBias <- aggregate(cbind(pop.cor, cor.bias, cor.RB) ~ par_names + condition + analType + level, 
                     data = sim1cor, FUN = mean)

ggplot(corBias[corBias$level == "Case level",], 
       mapping = aes(x = analType, y = cor.bias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Bias for case-level correlation parameters")

ggplot(corBias[corBias$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.bias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Bias for dyad-level correlation parameters")

#----

# correlation robust bias----
corRobBias <- aggregate(cor.est ~ par_names + pop.cor + condition + analType + level,
                        data = sim1cor, FUN = median)
corRobBias$cor.robbias <- corRobBias$cor.est - corRobBias$pop.cor
corRobBias$cor.robRB <- corRobBias$cor.robbias / corRobBias$pop.cor

ggplot(corRobBias[corRobBias$level == "Case level",], 
       mapping = aes(x = analType, y = cor.robbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Robust Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust bias for case-level correlation parameters")

ggplot(corRobBias[corRobBias$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.robbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Robust Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust bias for dyad-level correlation parameters")

#----

# correlation SE----
corSE_est <- aggregate(cor.SE ~ par_names + pop.cor + condition + analType + level,
                       data = sim1cor, FUN = mean)
names(corSE_est)[names(corSE_est) == "cor.SE"] <- "cor.estSE"

ggplot(corSE_est[corSE_est$level == "Case level",], 
       mapping = aes(x = analType, y = cor.estSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Standard errors of case-level correlation parameters")

ggplot(corSE_est[corSE_est$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.estSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Standard errors of dyad-level correlation parameters")

#----

# correlation SE bias----
corSE_obs <- aggregate(cor.est ~ par_names + pop.cor + condition + analType + level,
                       data = sim1cor, FUN = sd)
names(corSE_obs)[names(corSE_obs) == "cor.est"] <- "cor.obsSE"

corSEbias <- merge(corSE_est, corSE_obs)
corSEbias$cor.SEbias <- corSEbias$cor.estSE - corSEbias$cor.obsSE
corSEbias$cor.SERB <- corSEbias$cor.SEbias / corSEbias$cor.obsSE

ggplot(corSEbias[corSEbias$level == "Case level",], 
       mapping = aes(x = analType, y = cor.SEbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("SE bias for case-level correlation parameters")

ggplot(corSEbias[corSEbias$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.SEbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("SE bias for dyad-level correlation parameters")

#----

# correlation coverage rate----
sim1cor$cor.coverage <- sim1cor$pop.cor >= sim1cor$cor.low & sim1cor$pop.cor <= sim1cor$cor.up

corCR <- aggregate(cor.coverage ~ par_names + pop.cor + condition + analType + level,
                   data = sim1cor, FUN = function(x) mean(x)*100)

ggplot(corCR[corCR$level == "Case level",], 
       mapping = aes(x = analType, y = cor.coverage, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 95) + xlab("Method") + ylab("Coverage Rate") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Coverage rate for case-level correlation parameters")

ggplot(corCR[corCR$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.coverage, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 95) + xlab("Method") + ylab("Coverage Rate") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Coverage rate for dyad-level correlation parameters")


#----

# correlation RMSE----
sim1cor$cor.sqbias <- (sim1cor$cor.bias)^2

corRMSE <- aggregate(cor.sqbias ~ par_names + pop.cor + condition + analType + level,
                     data = sim1cor, FUN = function (x) sqrt(mean(x)))
names(corRMSE)[names(corRMSE) == "cor.sqbias"] <- "cor.RMSE"

#FIXME should I compute RMSE as bias^2 + SE^2

ggplot(corRMSE[corRMSE$level == "Case level",], 
       mapping = aes(x = analType, y = cor.RMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("RMSE for case-level correlation parameters")

ggplot(corRMSE[corRMSE$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.RMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("RMSE for dyad-level correlation parameters")

#----

# correlation median absolute deviation (MAD)----
corMAD <- aggregate(cor.est ~ par_names + pop.cor + condition + analType + level,
                    data = sim1cor, FUN = function (x) 1.4826*median(abs(x - median(x))))
names(corMAD)[names(corMAD) == "cor.est"] <- "cor.MAD"

ggplot(corMAD[corMAD$level == "Case level",], 
       mapping = aes(x = analType, y = cor.MAD, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Deviation") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MAD for case-level correlation parameters")

ggplot(corMAD[corMAD$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.MAD, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Deviation") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MAD for dyad-level correlation parameters")

#----

# correlation robust RMSE (robRMSE)----

### RobRMSE = (RobBias^2 + MAD^2) 
corMAD$cor.MADsq <- (corMAD$cor.MAD)^2
corRobRMSE <- merge(corRobBias, corMAD)

corRobRMSE$cor.RobRMSE <- (corRobRMSE$cor.robbias)^2 + corRobRMSE$cor.MADsq

ggplot(corRobRMSE[corRobRMSE$level == "Case level",], 
       mapping = aes(x = analType, y = cor.RobRMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Robust RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust RMSE for case-level correlation parameters")

ggplot(corRobRMSE[corRobRMSE$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.RobRMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Robust RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust RMSE for dyad-level correlation parameters")

#----

# correlation median absolute error (MedAE)----

### medAE = median(abs(est - pop))
corMedAE <- aggregate(cor.bias ~ par_names + pop.cor + condition + analType + level,
                      data = sim1cor, FUN = function(x) median(abs(x)))
names(corMedAE)[names(corMedAE) == "cor.bias"] <- "cor.MedAE"

ggplot(corMedAE[corMedAE$level == "Case level",], 
       mapping = aes(x = analType, y = cor.MedAE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Error") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MedAE for case-level correlation parameters")

ggplot(corMedAE[corMedAE$level == "Dyad level",], 
       mapping = aes(x = analType, y = cor.MedAE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Error") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MedAE for dyad-level correlation parameters")

#----

################################
# POSTPROCSSING FOR SD ESTIMATES
################################

# SD bias----
sim1SD$SD.bias <- sim1SD$SD.est - sim1SD$pop.SD
sim1SD$SD.RB <- sim1SD$SD.bias / sim1SD$pop.SD

SDBias <- aggregate(cbind(pop.SD, SD.bias, SD.RB) ~ par_names + condition + analType + level, 
                     data = sim1SD, FUN = mean)

ggplot(SDBias[SDBias$level == "Case level",], 
       mapping = aes(x = analType, y = SD.bias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Bias for case-level SD parameters")

ggplot(SDBias[SDBias$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.bias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Bias for dyad-level SD parameters")

#----

# SD robust bias----
SDRobBias <- aggregate(SD.est ~ par_names + pop.SD + condition + analType + level,
                        data = sim1SD, FUN = median)
SDRobBias$SD.robbias <- SDRobBias$SD.est - SDRobBias$pop.SD
SDRobBias$SD.robRB <- SDRobBias$SD.robbias / SDRobBias$pop.SD

ggplot(SDRobBias[SDRobBias$level == "Case level",], 
       mapping = aes(x = analType, y = SD.robbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Robust Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust bias for case-level SD parameters")

ggplot(SDRobBias[SDRobBias$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.robbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Robust Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust bias for dyad-level SD parameters")

#----

# SD SE----
SDSE_est <- aggregate(SD.SE ~ par_names + pop.SD + condition + analType + level,
                       data = sim1SD, FUN = mean)
names(SDSE_est)[names(SDSE_est) == "SD.SE"] <- "SD.estSE"

ggplot(SDSE_est[SDSE_est$level == "Case level",], 
       mapping = aes(x = analType, y = SD.estSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Standard errors of case-level SD parameters")

ggplot(SDSE_est[SDSE_est$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.estSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Standard errors of dyad-level SD parameters")

#----

# SD SE bias----
SDSE_obs <- aggregate(SD.est ~ par_names + pop.SD + condition + analType + level,
                       data = sim1SD, FUN = sd)
names(SDSE_obs)[names(SDSE_obs) == "SD.est"] <- "SD.obsSE"

SDSEbias <- merge(SDSE_est, SDSE_obs)
SDSEbias$SD.SEbias <- SDSEbias$SD.estSE - SDSEbias$SD.obsSE
SDSEbias$SD.SERB <- SDSEbias$SD.SEbias / SDSEbias$SD.obsSE

ggplot(SDSEbias[SDSEbias$level == "Case level",], 
       mapping = aes(x = analType, y = SD.SEbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("SE bias for case-level SD parameters")

ggplot(SDSEbias[SDSEbias$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.SEbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("SE bias for dyad-level SDn parameters")

#----

# SD coverage rate----
sim1SD$SD.coverage <- sim1SD$pop.SD >= sim1SD$SD.low & sim1SD$pop.SD <= sim1SD$SD.up

SDCR <- aggregate(SD.coverage ~ par_names + pop.SD + condition + analType + level,
                   data = sim1SD, FUN = function(x) mean(x)*100)

ggplot(SDCR[SDCR$level == "Case level",], 
       mapping = aes(x = analType, y = SD.coverage, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 95) + xlab("Method") + ylab("Coverage Rate") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Coverage rate for case-level SD parameters")

ggplot(SDCR[SDCR$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.coverage, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 95) + xlab("Method") + ylab("Coverage Rate") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Coverage rate for dyad-level SD parameters")

#----

# SD RMSE----
sim1SD$SD.sqbias <- (sim1SD$SD.bias)^2

SDRMSE <- aggregate(SD.sqbias ~ par_names + pop.SD + condition + analType + level,
                     data = sim1SD, FUN = function (x) sqrt(mean(x)))
names(SDRMSE)[names(SDRMSE) == "SD.sqbias"] <- "SD.RMSE"

#FIXME should I compute RMSE as bias^2 + SE^2

ggplot(SDRMSE[SDRMSE$level == "Case level",], 
       mapping = aes(x = analType, y = SD.RMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("RMSE for case-level SD parameters")

ggplot(SDRMSE[SDRMSE$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.RMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("RMSE for dyad-level SD parameters")


#----

# SD median absolute deviation (MAD)----
SDMAD <- aggregate(SD.est ~ par_names + pop.SD + condition + analType + level,
                   data = sim1SD, FUN = function (x) 1.4826*median(abs(x - median(x))))
names(SDMAD)[names(SDMAD) == "SD.est"] <- "SD.MAD"

ggplot(SDMAD[SDMAD$level == "Case level",], 
       mapping = aes(x = analType, y = SD.MAD, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Deviation") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MAD for case-level SD parameters")

ggplot(SDMAD[SDMAD$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.MAD, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Deviation") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MAD for dyad-level SD parameters")

#----

# SD robust RMSE (robRMSE)----

### RobRMSE = (RobBias^2 + MAD^2) 
SDMAD$SD.MADsq <- (SDMAD$SD.MAD)^2
SDRobRMSE <- merge(SDRobBias, SDMAD)

SDRobRMSE$SD.RobRMSE <- (SDRobRMSE$SD.robbias)^2 + SDRobRMSE$SD.MADsq

ggplot(SDRobRMSE[SDRobRMSE$level == "Case level",], 
       mapping = aes(x = analType, y = SD.RobRMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Robust RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust RMSE for case-level SD parameters")

ggplot(SDRobRMSE[SDRobRMSE$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.RobRMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Robust RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust RMSE for dyad-level SD parameters")

#----

# SD median absolute error (MedAE)----

### medAE = median(abs(est - pop))
SDMedAE <- aggregate(SD.bias ~ par_names + pop.SD + condition + analType + level,
                     data = sim1SD, FUN = function(x) median(abs(x)))
names(SDMedAE)[names(SDMedAE) == "SD.bias"] <- "SD.MedAE"

ggplot(SDMedAE[SDMedAE$level == "Case level",], 
       mapping = aes(x = analType, y = SD.MedAE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Error") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MedAE for case-level SD parameters")

ggplot(SDMedAE[SDMedAE$level == "Dyad level",], 
       mapping = aes(x = analType, y = SD.MedAE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Error") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MedAE for dyad-level SD parameters")

#----

########################################
# POSTPROCSSING FOR COVARIANCE ESTIMATES
########################################

# covariance bias----
sim1cov$cov.bias <- sim1cov$cov.est - sim1cov$pop.cov
sim1cov$cov.RB <- sim1cov$cov.bias / sim1cov$pop.cov

covBias <- aggregate(cbind(pop.cov, cov.bias, cov.RB) ~ par_names + condition + analType + level, 
                     data = sim1cov, FUN = mean)

ggplot(covBias[covBias$level == "Case level",], 
       mapping = aes(x = analType, y = cov.bias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Bias for case-level covariance parameters")

ggplot(covBias[covBias$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.bias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Bias for dyad-level covariance parameters")

#----

# covariance robust bias----
covRobBias <- aggregate(cov.est ~ par_names + pop.cov + condition + analType + level,
                        data = sim1cov, FUN = median)
covRobBias$cov.robbias <- covRobBias$cov.est - covRobBias$pop.cov
covRobBias$cov.robRB <- covRobBias$cov.robbias / covRobBias$pop.cov

ggplot(covRobBias[covRobBias$level == "Case level",], 
       mapping = aes(x = analType, y = cov.robbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust bias for case-level covariance parameters")

ggplot(covRobBias[covRobBias$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.robbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Bias") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust bias for dyad-level covariance parameters")

#----

# covariance SE----
covSE_est <- aggregate(cov.SE ~ par_names + pop.cov + condition + analType + level,
                       data = sim1cov, FUN = mean)
names(covSE_est)[names(covSE_est) == "cov.SE"] <- "cov.estSE"

ggplot(covSE_est[covSE_est$level == "Case level",], 
       mapping = aes(x = analType, y = cov.estSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Standard errors of case-level covariance parameters")

ggplot(covSE_est[covSE_est$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.estSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Standard errors of dyad-level covariance parameters")

#----

# covariance SE bias----
covSE_obs <- aggregate(cov.est ~ par_names + pop.cov + condition + analType + level,
                       data = sim1cov, FUN = sd)
names(covSE_obs)[names(covSE_obs) == "cov.est"] <- "cov.obsSE"

covSEbias <- merge(covSE_est, covSE_obs)
covSEbias$cov.SEbias <- covSEbias$cov.estSE - covSEbias$cov.obsSE
covSEbias$cov.SERB <- covSEbias$cov.SEbias / covSEbias$cov.obsSE

ggplot(covSEbias[covSEbias$level == "Case level",], 
       mapping = aes(x = analType, y = cov.SEbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("SE bias for case-level covariance parameters")

ggplot(covSEbias[covSEbias$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.SEbias, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 0) + xlab("Method") + ylab("Standard Errors") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("SE bias for dyad-level covariance parameters")

#----

# covariance coverage rate----
sim1cov$cov.coverage <- sim1cov$pop.cov >= sim1cov$cov.low & sim1cov$pop.cov <= sim1cov$cov.up

covCR <- aggregate(cov.coverage ~ par_names + pop.cov + condition + analType + level,
                   data = sim1cov, FUN = function(x) mean(x)*100)

ggplot(covCR[covCR$level == "Case level",], 
       mapping = aes(x = analType, y = cov.coverage, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 95) + xlab("Method") + ylab("Coverage Rate") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Coverage rate for case-level covariance parameters")

ggplot(covCR[covCR$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.coverage, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 95) + xlab("Method") + ylab("Coverage Rate") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Coverage rate for dyad-level covariance parameters")

#----

# covariance RMSE----
sim1cov$cov.sqbias <- (sim1cov$cov.bias)^2

covRMSE <- aggregate(cov.sqbias ~ par_names + pop.cov + condition + analType + level,
                     data = sim1cov, FUN = function (x) sqrt(mean(x)))
names(covRMSE)[names(covRMSE) == "cov.sqbias"] <- "cov.RMSE"

#FIXME should I compute RMSE as bias^2 + SE^2

ggplot(covRMSE[covRMSE$level == "Case level",], 
       mapping = aes(x = analType, y = cov.RMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("RMSE for case-level covariance parameters")

ggplot(covRMSE[covRMSE$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.RMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("RMSE for dyad-level covariance parameters")

#----

# covariance median absolute deviation (MAD)----
covMAD <- aggregate(cov.est ~ par_names + pop.cov + condition + analType + level,
                     data = sim1cov, FUN = function (x) 1.4826*median(abs(x - median(x))))
names(covMAD)[names(covMAD) == "cov.est"] <- "cov.MAD"

ggplot(covMAD[covMAD$level == "Case level",], 
       mapping = aes(x = analType, y = cov.MAD, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Deviation") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MAD for case-level covariance parameters")

ggplot(covMAD[covMAD$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.MAD, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Deviation") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MAD for dyad-level covariance parameters")

#----

# covariance robust RMSE (robRMSE)----

### RobRMSE = (RobBias^2 + MAD^2) 
covMAD$cov.MADsq <- (covMAD$cov.MAD)^2
covRobRMSE <- merge(covRobBias, covMAD)

covRobRMSE$cov.RobRMSE <- (covRobRMSE$cov.robbias)^2 + covRobRMSE$cov.MADsq

ggplot(covRobRMSE[covRobRMSE$level == "Case level",], 
       mapping = aes(x = analType, y = cov.RobRMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Robust RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust RMSE for case-level covariance parameters")

ggplot(covRobRMSE[covRobRMSE$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.RobRMSE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Robust RMSE") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("Robust RMSE for dyad-level covariance parameters")

#----

# covariance median absolute error (MedAE)----

### medAE = median(abs(est - pop))
covMedAE <- aggregate(cov.bias ~ par_names + pop.cov + condition + analType + level,
                       data = sim1cov, FUN = function(x) median(abs(x)))
names(covMedAE)[names(covMedAE) == "cov.bias"] <- "cov.MedAE"

ggplot(covMedAE[covMedAE$level == "Case level",], 
       mapping = aes(x = analType, y = cov.MedAE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = case_pars, labels = case_labels) + 
  scale_shape_manual(values = case_shape, labels = case_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Error") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MedAE for case-level covariance parameters")

ggplot(covMedAE[covMedAE$level == "Dyad level",], 
       mapping = aes(x = analType, y = cov.MedAE, group = par_names, color = par_names,
                     shape = par_names)) + geom_point() + 
  facet_wrap(~condition, nrow = 2, ncol = 4) + 
  scale_color_manual(values = dyad_pars, labels = dyad_labels) + 
  scale_shape_manual(values = dyad_shape, labels = dyad_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Method") + ylab("Median Absolute Error") +
  labs(color = "Parameters", shape = "Parameters") + theme_grey(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  ggtitle("MedAE for dyad-level covariance parameters")

#----


# foo[which(foo$cov.MCMC_neff < 100), ]

sim1mPSRF[which(sim1mPSRF$mPSRF > 1.05),] -> foo

# > unique(sim1mPSRF$analType)
# [1] Diffuse priors Pr-0.05 priors Pr-0.1 priors  Pr-0.2 priors 
# Levels: Diffuse priors Pr-0.05 priors Pr-0.1 priors Pr-0.2 priors FIML

dim(foo[foo$condition == "6-10" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "6-10" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "6-10" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "6-10" & foo$analType == "Pr-0.2 priors", ])

dim(foo[foo$condition == "6-25" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "6-25" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "6-25" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "6-25" & foo$analType == "Pr-0.2 priors", ])

dim(foo[foo$condition == "8-10" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "8-10" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "8-10" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "8-10" & foo$analType == "Pr-0.2 priors", ])

dim(foo[foo$condition == "8-25" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "8-25" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "8-25" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "8-25" & foo$analType == "Pr-0.2 priors", ])

dim(foo[foo$condition == "10-10" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "10-10" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "10-10" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "10-10" & foo$analType == "Pr-0.2 priors", ])

dim(foo[foo$condition == "10-25" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "10-25" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "10-25" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "10-25" & foo$analType == "Pr-0.2 priors", ])

dim(foo[foo$condition == "20-10" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "20-10" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "20-10" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "20-10" & foo$analType == "Pr-0.2 priors", ])

dim(foo[foo$condition == "20-25" & foo$analType == "Diffuse priors", ])
dim(foo[foo$condition == "20-25" & foo$analType == "Pr-0.05 priors", ])
dim(foo[foo$condition == "20-25" & foo$analType == "Pr-0.1 priors", ])
dim(foo[foo$condition == "20-25" & foo$analType == "Pr-0.2 priors", ])
