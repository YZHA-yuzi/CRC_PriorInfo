########################################################################
# ------------ R codes used for analyzing 3-stream HIV CRC data ------------- #
# NOTE:
# Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" is located.
########################################################################

library(bbmle)
library(MCMCpack)
library(ggplot2)

## read in R functions ##
source(paste0("FUNs.R"))


## read in the 3-stream CRC data 
df_aggre <- data.frame(n111 = 14, n110 = 66, n101 = 23, n100 = 1072,
                       n011 = 58, n010 = 729, n001 = 319)
nc.curr <- as.numeric(sum(df_aggre))
# prepare data for log-linear model #
df_logl <- gen_dflogl(dat = df_aggre, 
                      nstreams = 3, 
                      predictors = c(paste0("X", 1:3), 
                                     "X1X2", "X1X3", "X2X3", "X1X2X3"))


# fit the log-linear model selected in the paper Poorolajal et al. (2017)
fit.logl <- glm(count ~ X1 + X2 + X3 + X1X2 + X2X3, 
                data = df_logl, family = "poisson")
sum.fit <- summary(fit.logl)$coefficients
se.int <- as.numeric(sum.fit[1,2])
est.int <- as.numeric(sum.fit[1,1])
Nhat.logl <- as.numeric(ceiling(nc.curr + exp(est.int)))
CI.logl <- as.numeric(nc.curr + 
                        c(ceiling(exp(est.int - 1.96*se.int)),
                          ceiling(exp(est.int + 1.96*se.int))))
CI.logl[1] <- max(nc.curr, CI.logl[1])


# ## -------- Fit all possilbe log-linear models to identify the 
# ## model with the lowest AIC among 8 candidate models -------- ##
# # data set: Poorolajal 2017 paper #
# # stream 1 (transfusion center); stream 2 (VCTCs); stream 3 (Prison) 
# dat_S3_p <- data.frame(n111 = 14, n110 = 66, n101 = 23, n100 = 1072,
#                        n011 = 58, n010 = 729, n001 = 319)
# predictors_S3 <- c("X1", "X2", "X3",
#                    "X1X2", "X1X3", "X2X3",
#                    "X1X2X3")
# # (1). transfer aggregated cell counts to data frame used for fitting log-linear
# dflogl_S3_p <- gen_dflogl(dat_S3_p, 3, predictors = predictors_S3)
# # (2). fit possible log-linear models (with and without highest interaction)
# fitall_S3_p <- fit_logl(dat = dflogl_S3_p, nstreams = 3,
#                         predictors = predictors_S3)
# tab_S3_p <- fitall_S3_p$re_mat
# # (3). select the 8 candidate models
# pred.sel = paste0(c("X1, X2, X3"), c("", ", X1X3", ", X1X2", ", X2X3", 
#                                      ", X1X2, X1X3", ", X1X3, X2X3", ", X1X2, X2X3",
#                                      ", X1X2, X1X3, X2X3"))
# tab_S3_p <- subset(tab_S3_p, Model %in% pred.sel)
# # (4). find the model with the lowest AIC
# tab_S3_p.lowAIC <- tab_S3_p[order(tab_S3_p$AIC, decreasing = F), ]
# tab_S3_p.lowAIC[1, ]
# ## -------- END: Fit all possilbe log-linear models -------- ##

# log-linear model with lowest AIC
fit.logl.aic <- glm(count ~ X1 + X2 + X3 + X1X2 + X1X3 + X2X3, 
                    data = df_logl, family = "poisson")
sum.fit.aic <- summary(fit.logl.aic)$coefficients
se.int.aic <- as.numeric(sum.fit.aic[1,2])
est.int.aic <- as.numeric(sum.fit.aic[1,1])
Nhat.logl.aic <- as.numeric(ceiling(nc.curr + exp(est.int.aic)))
CI.logl.aic <- as.numeric(nc.curr + 
                            c(ceiling(exp(est.int.aic - 1.96*se.int.aic)),
                              ceiling(exp(est.int.aic + 1.96*se.int.aic))))
CI.logl.aic[1] <- max(nc.curr, CI.logl.aic[1])




# Fit the alternative framework under the assumption 
## p312bar, or p31bar2 = psi
get.psi.beta <- function(df_aggre, a = 0, b = 0, parm.name){
  if(parm.name == "p31bar2"){
    psi.hat <- (a + df_aggre$n011)/(a + b + df_aggre$n011 + df_aggre$n010)
  }else if(parm.name == "p312bar"){
    psi.hat <- (a + df_aggre$n101)/(a + b + df_aggre$n101 + df_aggre$n100)
  }else if(parm.name == "p312"){
    psi.hat <- (a + df_aggre$n111)/(a + b + df_aggre$n111 + df_aggre$n110)
  }
  return(psi.hat)
}

parm.hat.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = 2))
rownames(parm.hat.mat) <- c("p31bar2", "p312bar")
colnames(parm.hat.mat) <- c("psi_bcB")
for(i in c("p31bar2", "p312bar")){
  parm.hat.mat[i, 1] <- get.psi.beta(df_aggre = df_aggre, 
                                     a = 1, b = 0, parm.name = i)
}
parm.hat.mat$phi <- c(2, 1)

### get point estimate under constraints 
### (1) psi = p312bar/phi, phi = 1
### (2) psi = p31bar2/phi, phi = 2
Nhat.bcB <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
colnames(Nhat.bcB) <- c("Assump", "psi_bcB", "phi", "est")
Nhat.bcB$Assump <- c("p31bar2", "p312bar")
Nhat.bcB[,2:3] <- parm.hat.mat[,1:2]
Nhat.bcB$est <- sapply(1:2, function(i) 
  ceiling(get_Nhatpsi_S3(dat = df_aggre, 
                         psi = parm.hat.mat[,"psi_bcB"][i]/
                           parm.hat.mat[,"phi"][i])))

## get Dirichlet-based CI for estimates under constraints
### (1) psi = p312bar/phi, phi = 1
### (2) psi = p31bar2/phi, phi = 2
### The bias-correction based on Beta(1,0) priors were implemented
set.seed(12345)
CI.phi.bcB.list <- list()
count = 1
for(l in c(3, 2)){
  inter <- Dir_CI_phi_S3_beta(df_aggre = df_aggre, n.post = 10000, 
                              phi.type = paste0("phi", l), 
                              phi = parm.hat.mat$phi[count],
                              a = 1, b = 0) 
  inter$CI[1] <- max(nc.curr, inter$CI[1])
  CI.phi.bcB.list[[count]] <- ceiling(inter$CI)
  count = count + 1
}
names(CI.phi.bcB.list) <- c("p31bar2", "p312bar")
CI.phi.bcB.list$p31bar2
CI.phi.bcB.list$p312bar



###### RR-type assumption ######
## get estimated N under RR-type assumption
## using the alternative model
get.betapost <- function(a, b, df_aggre, r = 1, assump.type = "RR"){
  sum.ab = a + b
  p312.post <- (a + df_aggre$n111)/(sum.ab + 
                                      df_aggre$n111 + df_aggre$n110)
  p31bar2.post <- (a + df_aggre$n011)/(sum.ab + 
                                         df_aggre$n011 + df_aggre$n010)
  p312bar.post <- (a + df_aggre$n101)/(sum.ab + 
                                         df_aggre$n101 + df_aggre$n100)
  if(assump.type == "RR"){
    psi.post <- r*(p31bar2.post*p312bar.post)/p312.post
  }else if(assump.type == "OR"){
    inter.or.post <- r*(p31bar2.post/(1-p31bar2.post))*
      (p312bar.post/(1-p312bar.post))/(p312.post/(1-p312.post))
    psi.post <- inter.or.post/(1 + inter.or.post)
  }
  return(c(p312.post, p312bar.post, p31bar2.post, psi.post))
}

nc.rm = sum(df_aggre[1,1:6])
ratio.val.vec <- c(1, 2)
psi.hat.RR <- Nhat.RR <- rep(NA, length(ratio.val.vec))
CI.RR.mat <- matrix(NA, ncol = 2, nrow = length(ratio.val.vec))
for(i in 1:length(ratio.val.vec)){
  psi.hat.RR[i] <- get.betapost(a = 1, b = 0, df_aggre = df_aggre,
                                r = ratio.val.vec[i],
                                assump.type = "RR")[4]
  Nhat.RR[i] <- as.numeric(ceiling(nc.rm + 
                                     c(df_aggre[1,7])/psi.hat.RR[i]))
  res.CI.alter <- Dir_CI_RR(df_aggre = df_aggre, 
                            n.post = 10000,
                            rr = ratio.val.vec[i], a = 1, b = 0)
  CI.RR <- res.CI.alter$CI
  CI.RR[1] <- max(nc.curr, CI.RR[1])
  CI.RR <- ceiling(CI.RR)
  CI.RR.mat[i, ] <- CI.RR
}



###### OR-type assumption ######
## get estimated N under OR-type assumption
## using the alternative model
nc.rm = sum(df_aggre[1,1:6])
ratio.val.vec <- c(1, 2)
psi.hat.OR <- Nhat.OR <- rep(NA, length(ratio.val.vec))
CI.OR.mat <- matrix(NA, ncol = 2, nrow = length(ratio.val.vec))
for(i in 1:length(ratio.val.vec)){
  psi.hat.OR[i] <- get.betapost(a = 1, b = 0, df_aggre = df_aggre,
                                r = ratio.val.vec[i],
                                assump.type = "OR")[4]
  Nhat.OR[i] <- as.numeric(ceiling(nc.rm + 
                                     c(df_aggre[1,7])/psi.hat.OR[i]))
  res.CI.alter <- Dir_CI_OR(df_aggre = df_aggre, 
                            n.post = 10000,
                            or = ratio.val.vec[i], a = 1, b = 0)
  CI.OR <- res.CI.alter$CI
  CI.OR[1] <- max(nc.curr, CI.OR[1])
  CI.OR <- ceiling(CI.OR)
  CI.OR.mat[i, ] <- CI.OR
}



####### Uncertainty Analysis ########
### under constraints:
### (1) psi = p312bar/phi, phi ~ Unif(0.8, 1.2)
### (2) psi = p31bar2/phi, phi ~ Unif(1.6, 2.4)
a.vec <- c(1.6, 0.8)
b.vec <- c(2.4, 1.2)
CI.bcB.unce.list <- list()
phi.name.vec <- c("p31bar2", "p312bar")
for(l in 1:length(a.vec)){
  res.CI.i <- Dir_CI_phi_S3_unce(df_aggre = df_aggre, 
                                 n.post = 1000, 
                                 a = 1, b = 0,
                                 parm.name = phi.name.vec[l],
                                 assume.dist = "Uniform", 
                                 a.dist = a.vec[l], 
                                 b.dist = b.vec[l])
  res.CI.i$CI[1] <- max(nc.curr, res.CI.i$CI[1])
  CI.bcB.unce.list[[l]] <- ceiling(res.CI.i$CI)
  cat(paste0(phi.name.vec[l], ", a = ", a.vec[l], ", b = ", b.vec[l]))
}
names(CI.bcB.unce.list) <- c("p31bar2", "p312bar")

CI.phi.bcB.list
CI.bcB.unce.list





##### uncertainty analysis RR-type 
a.vec <- c(0.8, 1.6)
b.vec <- c(1.2, 2.4)
CI.RR.unce.list <- list()
for(i in 1:length(a.vec)){
  re.CI.RR <- Dir_CI_RR_unce(df_aggre = df_aggre, 
                             n.post = 1000,
                             assume.dist = "Uniform",
                             a.dist = a.vec[i], 
                             b.dist = b.vec[i])
  CI.RR.unce <- re.CI.RR$CI
  CI.RR.unce[1] <- max(nc.curr, CI.RR.unce[1])
  CI.RR.unce.list[[i]] <- ceiling(CI.RR.unce)
}




##### uncertainty analysis OR-type 
a.vec <- c(0.8, 1.6)
b.vec <- c(1.2, 2.4)
CI.OR.unce.list <- list()
for(i in 1:length(a.vec)){
  re.CI.OR <- Dir_CI_OR_unce(df_aggre = df_aggre, 
                             n.post = 1000,
                             assume.dist = "Uniform",
                             a.dist = a.vec[i], b.dist = b.vec[i])
  CI.OR.unce <- re.CI.OR$CI
  CI.OR.unce[1] <- max(nc.curr, CI.OR.unce[1])
  CI.OR.unce.list[[i]] <- ceiling(CI.OR.unce)
  cat(i, ",")
}

CI.OR.unce.list


### Summarize results ####
est.names <- c("p312bar", "p31bar2", "RR", "OR")
est.names.1 <- c(paste0(est.names[1:2], c("_1", "_2")),
                 c(sapply(est.names[3:4], function(x) 
                   paste0(x, "_", 1:2))))
est.val <- c(Nhat.bcB$est[sapply(1:2, function(x) 
  which(Nhat.bcB$Assump == est.names[x]))],
  Nhat.RR, Nhat.OR)
tab.all <- data.frame(assump = c("logl", "logl (AIC)",
                                 est.names.1,
                                 paste0(est.names.1, "_unce")),
                      N = c(Nhat.logl, Nhat.logl.aic, rep(est.val, 2)),
                      lci = c(CI.logl[1], CI.logl.aic[1],
                              c(t(sapply(CI.phi.bcB.list[c("p312bar", "p31bar2")], 
                                         function(x) x[1]))),
                              CI.RR.mat[,1],
                              CI.OR.mat[,1],
                              c(t(sapply(CI.bcB.unce.list[c("p312bar", "p31bar2")], 
                                         function(x) x[1]))),
                              sapply(CI.RR.unce.list, function(x) x[1]),
                              sapply(CI.OR.unce.list, function(x) x[1])),
                      uci = c(CI.logl[2], CI.logl.aic[2],
                              c(t(sapply(CI.phi.bcB.list[c("p312bar", "p31bar2")], 
                                         function(x) x[2]))),
                              CI.RR.mat[,2],
                              CI.OR.mat[,2],
                              c(t(sapply(CI.bcB.unce.list[c("p312bar", "p31bar2")], 
                                         function(x) x[2]))),
                              sapply(CI.RR.unce.list, function(x) x[2]),
                              sapply(CI.OR.unce.list, function(x) x[2])))

keep.name <- c("logl", "logl (AIC)",
               "p312bar_1", "p312bar_1_unce",
               "p31bar2_2", "p31bar2_2_unce",
               "RR_1", "RR_1_unce",
               "RR_2", "RR_2_unce",
               "OR_1", "OR_1_unce",
               "OR_2", "OR_2_unce")
tab.all.sub <- subset(tab.all, assump %in% keep.name)
tab.all.sub$assump <- factor(tab.all.sub$assump,
                             levels = keep.name)
xlab.vec = c("logl" = "log-linear",
             "logl (AIC)" = "log-linear (AIC)",
             "p312bar_1" = 
               expression(paste("p3|1", bar(2), 
                                "/", psi, "=1",sep = "")),
             "p312bar_1_unce" = 
               expression(paste("p3|1", bar(2), 
                                "/", psi, "~Unif(0.8,1.2)",sep = "")),
             "p31bar2_2" = 
               expression(paste("p3|", bar(1), "2",
                                "/", psi, "=2",sep = "")),
             "p31bar2_2_unce" = 
               expression(paste("p3|", bar(1), "2",
                                "/", psi, "~Unif(1.6,2.4)", sep = "")),
             "RR_1" = "RR=1",
             "RR_1_unce" = "RR~Unif(0.8,1.2)",
             "RR_2" = "RR=2",
             "RR_2_unce" = "RR~Unif(1.6,2.4)",
             "OR_1" = "OR = 1",
             "OR_1_unce" = "OR~Unif(0.8,1.2)",
             "OR_2" = "OR = 2",
             "OR_2_unce" = "OR~Unif(1.6,2.4)")

tab.all.sub <- tab.all.sub[order(tab.all.sub$assump), ]

tab.all.sub.1 <- data.frame(group = c("log", "log_OR_AIC",
                                      rep("p312bar", 2),
                                      rep("p31bar2", 2),
                                      rep("RR_1", 2),
                                      rep("RR_2", 2),
                                      rep("log_OR_AIC", 2),
                                      rep("OR", 2)),
                            tab.all.sub)
pl.3s.color <- ggplot(tab.all.sub.1, 
                      aes(x = N, y = assump, xmin = lci, xmax = uci,
                          color = group)) + 
  geom_pointrange(position = position_dodge(width=0.40)) + 
  geom_errorbar(position = position_dodge(width=0.40), 
                width = 0.1, linewidth = 0.8) + 
  geom_vline(xintercept = tab.all$N[1], linetype = "dashed",
             color = "grey", linewidth = 1) + 
  scale_y_discrete(labels = xlab.vec) + 
  theme_bw() + 
  scale_colour_brewer(palette = "Dark2") + 
  labs(y = "Method", 
       x = "Estimated N and 95% CI") +
  theme(axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.position = "none")
pl.3s.color

### save results
if(!file.exists("./Results")){
  dir.create("Results")
}
ggsave(filename = paste0("./Results/", "Fig_HIV3S.png"),
       height = 7, width = 9, units = "in", dpi = 600,
       plot = pl.3s.color)



