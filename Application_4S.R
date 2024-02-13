########################################################################
# ------------ R codes used for analyzing 4-stream HIV CRC data ------------- #
# NOTE:
# Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" is located.
########################################################################

library(bbmle)
library(MCMCpack)
library(ggplot2)


## read in R functions ##
source(paste0("FUNs.R"))

## read in 4-stream CRC data
dat.abeni = as.data.frame( expand.grid(c(1,0), c(1,0), c(1,0), c(1,0)) )
dat.abeni = dat.abeni[-16, ]
dat.abeni$count = c(0, 3, 1, 33, 0, 20, 6, 403, 3, 35, 10, 545, 11, 621, 205)
dat.abeni <- dat.abeni[order(dat.abeni$Var1, dat.abeni$Var2,
                             dat.abeni$Var3, dat.abeni$Var4, decreasing = T),]
pred_names <- apply(dat.abeni[,1:4], 1, 
                    function(x) paste0("n", paste(x, collapse = "")))
dat_S4_abeni <- as.data.frame(matrix(dat.abeni$count,
                                     nrow = 1))
colnames(dat_S4_abeni) <- apply(dat.abeni[,1:4], 1, 
                                function(x) 
                                  paste0("n", paste(x, collapse = "")))
predictors_S4 <- c(paste0("X", 1:4),
                   combn(1:4, 2, FUN = function(x) 
                     paste0(paste0("X",x), collapse = "")),
                   combn(1:4, 3, FUN = function(x) 
                     paste0(paste0("X",x), collapse = "")),
                   combn(1:4, 4, FUN = function(x) 
                     paste0(paste0("X",x), collapse = "")))
df_logl <- gen_dflogl(dat_S4_abeni, 4, predictors = predictors_S4)
df_aggre <- dat_S4_abeni




## the log-linear model fitted in the Abeni data 
fit.logl <- glm(count ~ X1 + X2 + X3 + X4 + X3X4, 
                data = df_logl, family = "poisson")
sum.fit <- summary(fit.logl)$coefficients
se.int <- as.numeric(sum.fit[1,2])
est.int <- as.numeric(sum.fit[1,1])
nc.curr <- as.numeric(sum(dat_S4_abeni))
Nhat.logl <- as.numeric(ceiling(nc.curr + exp(est.int)))
CI.logl <- as.numeric(nc.curr + 
                        c(ceiling(exp(est.int - 1.96*se.int)),
                          ceiling(exp(est.int + 1.96*se.int))))
CI.logl[1] <- max(nc.curr, CI.logl[1])




###### Analyze this 4-stream CRC data using the proposed model #######

## A function to get credible intervals 
get.psi.beta.S4 <- function(df_aggre, a = 0, b = 0, parm.name){
  if(parm.name == "p4123"){
    nu.val = df_aggre[["n1111"]]
    de.val = df_aggre[["n1110"]]
  }else if(parm.name == "p41bar23"){
    nu.val = df_aggre[["n0111"]]
    de.val = df_aggre[["n0110"]]
  }else if(parm.name == "p412bar3"){
    nu.val = df_aggre[["n1011"]]
    de.val = df_aggre[["n1010"]]
  }else if(parm.name == "p4123bar"){
    nu.val = df_aggre[["n1101"]]
    de.val = df_aggre[["n1100"]]
  }else if(parm.name == "p41bar2bar3"){
    nu.val = df_aggre[["n0011"]]
    de.val = df_aggre[["n0010"]]
  }else if(parm.name == "p41bar23bar"){
    nu.val = df_aggre[["n0101"]]
    de.val = df_aggre[["n0100"]]
  }else if(parm.name == "p412bar3bar"){
    nu.val = df_aggre[["n1001"]]
    de.val = df_aggre[["n1000"]]
  }
  psi.hat <- (nu.val + a)/(a + b + nu.val + de.val)
  return(psi.hat)
}

## get point estimates under constraints 
## psi = p4123bar/phi, phi = 1
## psi = p41bar23bar/phi, phi = 1
## psi = p412bar3bar/phi, phi = 1
## psi = p412bar3/phi, phi = 2
parm.vec <- c("p4123bar", "p41bar23bar", "p412bar3bar", "p412bar3")
parm.hat.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = length(parm.vec)))
rownames(parm.hat.mat) <- parm.vec
colnames(parm.hat.mat) <- c("psi", "psi_bcB")
for(i in parm.vec){
  parm.hat.mat[i, 1] <- get.psi.beta.S4(df_aggre = df_aggre, 
                                        a = 0, b = 0, parm.name = i)
  parm.hat.mat[i, 2] <- get.psi.beta.S4(df_aggre = df_aggre, 
                                        a = 1, b = 0, parm.name = i)
}
parm.hat.mat$phi <- c(1, 1, 1, 2)
phi.val.vec <- c(1, 1, 1, 2)

Nhat.bcB <- as.data.frame(matrix(NA, ncol = 5, nrow = length(parm.vec)))
colnames(Nhat.bcB) <- c("Assump", "psi", "psi_bcB", "phi", "est")
Nhat.bcB[,2:4] <- parm.hat.mat
Nhat.bcB[,1] <- rownames(parm.hat.mat)
nc.rm = as.numeric(sum(dat_S4_abeni[!colnames(dat_S4_abeni) %in% "n0001"]))
Nhat.bcB[, "est"] <- sapply(1:length(parm.vec), function(i) 
  ceiling(nc.rm + as.numeric(dat_S4_abeni$n0001)/
            (parm.hat.mat[,"psi_bcB"][i]/parm.hat.mat[,"phi"][i])))

## get Dirichlet-based intervals under constraints
## psi = p4123bar/phi, phi = 1
## psi = p41bar23bar/phi, phi = 1
## psi = p412bar3bar/phi, phi = 1
## psi = p412bar3/phi, phi = 2
set.seed(12345)
CI.phi.bcB.list <- list()
for(i in 1:length(parm.vec)){
  nc.curr = as.numeric(sum(df_aggre))
  res.CI.l <- Dir_CI_phi_S4_beta(df_aggre = df_aggre, n.post = 10000, 
                                 parm.name = parm.vec[i], 
                                 ratio.val = phi.val.vec[i],
                                 a = 1, b = 0)
  res.CI.l$CI[1] <- max(nc.curr, res.CI.l$CI[1])
  CI.phi.bcB.list[[i]] <- ceiling(res.CI.l$CI)
  cat(i, ",")
}
names(CI.phi.bcB.list) <- parm.vec

#### RR-type assumptions #####
phi.val.vec <- c(2.5)
parm.vec.RR <- c("RR_1bar", "RR_2bar")
Nhat.RR.bcB <- ratio.val.RR <- 
  matrix(NA, ncol = length(parm.vec.RR), nrow = length(phi.val.vec))
colnames(Nhat.RR.bcB) <- colnames(ratio.val.RR) <- 
  paste0(parm.vec.RR, "_bcB")
nc.rm = as.numeric(sum(dat_S4_abeni[!colnames(dat_S4_abeni) %in% "n0001"]))
for(l in 1:length(phi.val.vec)){
  phi.l = phi.val.vec[l]
  for(i in 1:length(parm.vec.RR)){
    RR.i = parm.vec.RR[i]
    if(RR.i == "RR_1bar"){
      p41bar23.i = get.psi.beta.S4(df_aggre = df_aggre, 
                                   a = 1, b = 0, parm.name = "p41bar23")
      p41bar23bar.i = get.psi.beta.S4(df_aggre = df_aggre, 
                                      a = 1, b = 0, parm.name = "p41bar23bar")
      p41bar2bar3.i = get.psi.beta.S4(df_aggre = df_aggre, 
                                      a = 1, b = 0, parm.name = "p41bar2bar3")
      psi.hat = phi.l*p41bar23bar.i*p41bar2bar3.i/p41bar23.i
      ratio.i = p41bar23.i/p41bar23bar.i/phi.l
    }else if(RR.i == "RR_2bar"){
      p412bar3.i = get.psi.beta.S4(df_aggre = df_aggre, 
                                   a = 1, b = 0, parm.name = "p412bar3")
      p412bar3bar.i = get.psi.beta.S4(df_aggre = df_aggre, 
                                      a = 1, b = 0, parm.name = "p412bar3bar")
      p41bar2bar3.i = get.psi.beta.S4(df_aggre = df_aggre, 
                                      a = 1, b = 0, parm.name = "p41bar2bar3")
      psi.hat = phi.l*p412bar3bar.i*p41bar2bar3.i/p412bar3.i
      ratio.i = p412bar3.i/p412bar3bar.i/phi.l
    }
    ratio.val.RR[l, i] <- ratio.i 
    Nhat.bcB.l <- ceiling(nc.rm + as.numeric(dat_S4_abeni$n0001)/psi.hat)
    Nhat.RR.bcB[l, i] <- Nhat.bcB.l
  }
}
rownames(Nhat.RR.bcB) <- paste0("r=", phi.val.vec)

set.seed(12346)
CI.RR.bcB.list <- list()
for(i in 1:length(phi.val.vec)){
  CI.bcB.mat <- matrix(NA, ncol = 2, nrow = length(parm.vec.RR))
  rownames(CI.bcB.mat) <- paste0(parm.vec.RR, "_bcB")
  for(l in 1:length(parm.vec.RR)){
    nc.curr = as.numeric(sum(df_aggre))
    res.CI.l <- Dir_CI_RR_S4_beta(df_aggre = df_aggre, n.post = 10000, 
                                  parm.name = parm.vec.RR[l], 
                                  ratio.val = phi.val.vec[i],
                                  a = 1, b = 0)
    res.CI.l$CI[1] <- max(nc.curr, res.CI.l$CI[1])
    CI.bcB.mat[l, ] <- ceiling(res.CI.l$CI)
  }
  CI.RR.bcB.list[[i]] <- CI.bcB.mat
  cat(i, ",")
}
names(CI.RR.bcB.list) <- paste0("r=", phi.val.vec)


#### Uncertainty analysis ######
set.seed(12347)
CI.phi.bcB.unce.list <- list()
a.vec <- c(rep(0.8, 3), 1.6)
b.vec <- c(rep(1.2, 3), 2.4)
for(i in 1:length(parm.vec)){
  nc.curr = as.numeric(sum(df_aggre))
  res.CI.l <- Dir_CI_phi_S4_unce(df_aggre = df_aggre, n.post = 1000, 
                                 parm.name = parm.vec[i], 
                                 a = 1, b = 0,
                                 assume.dist = "Uniform",
                                 a.dist = a.vec[i], b.dist = b.vec[i])
  res.CI.l$CI[1] <- max(nc.curr, res.CI.l$CI[1])
  CI.phi.bcB.unce.list[[i]] <- ceiling(res.CI.l$CI)
  cat(i, ",")
}
names(CI.phi.bcB.unce.list) <- parm.vec



### RR-type assumption uncertainty 
CI.RR.bcB.unce.list <- list()
a.vec <- c(2.25)
b.vec <- c(2.75)
for(i in 1:length(a.vec)){
  CI.bcB.mat <- matrix(NA, ncol = 2, nrow = length(parm.vec.RR))
  rownames(CI.bcB.mat) <- paste0(parm.vec.RR, "_bcB")
  for(l in 1:length(parm.vec.RR)){
    nc.curr = as.numeric(sum(df_aggre))
    res.CI.l <- Dir_CI_RR_S4_unce(df_aggre = df_aggre, n.post = 1000, 
                                  parm.name = parm.vec.RR[l], 
                                  a = 1, b = 0,
                                  assume.dist = "Uniform",
                                  a.dist = a.vec[i], b.dist = b.vec[i])
    res.CI.l$CI[1] <- max(nc.curr, res.CI.l$CI[1])
    CI.bcB.mat[l, ] <- ceiling(res.CI.l$CI)
    cat(l, ",")
  }
  CI.RR.bcB.unce.list[[i]] <- CI.bcB.mat
  cat(i, ",")
}
names(CI.RR.bcB.unce.list) <- paste0("r=", phi.val.vec)


### generates plots ####
est.names.phi <- c("p4123bar", "p41bar23bar", "p412bar3bar", "p412bar3")
est.names.phi.1 <- paste0(est.names.phi, "_", c(rep(1, 3), 2))
est.names.RR <- c("RR_1bar", "RR_2bar")
num.RR.assum <- as.character(c(2.5))
est.names.RR.1 <- paste0(est.names.RR, "_", num.RR.assum)

tab.all <- data.frame(assump = c("logl", 
                                 c(est.names.phi.1, est.names.RR.1),
                                 paste0(c(est.names.phi.1, est.names.RR.1), 
                                        "_unce")),
                      N = c(Nhat.logl, 
                            rep(c(c(Nhat.bcB$est),
                                  c(Nhat.RR.bcB)), 2)),
                      lci = c(CI.logl[1], 
                              c(t(sapply(CI.phi.bcB.list, function(x) x[1]))),
                              c(t(sapply(CI.RR.bcB.list, function(x) x[,1]))),
                              c(t(sapply(CI.phi.bcB.unce.list, function(x) x[1]))),
                              c(t(sapply(CI.RR.bcB.unce.list, function(x) x[,1])))),
                      uci = c(CI.logl[2], 
                              c(t(sapply(CI.phi.bcB.list, function(x) x[2]))),
                              c(t(sapply(CI.RR.bcB.list, function(x) x[,2]))),
                              c(t(sapply(CI.phi.bcB.unce.list, function(x) x[2]))),
                              c(t(sapply(CI.RR.bcB.unce.list, function(x) x[,2])))))

### Summarize results ####
keep.name.1 <- c("logl",
                 "p4123bar_1", "p4123bar_1_unce",
                 "p41bar23bar_1", "p41bar23bar_1_unce",
                 "p412bar3bar_1", "p412bar3bar_1_unce",
                 "p412bar3_2", "p412bar3_2_unce")
keep.name.RR <- c("RR_1bar_2.5", "RR_1bar_2.5_unce",
                  "RR_2bar_2.5", "RR_2bar_2.5_unce")
keep.name <- c(keep.name.1, keep.name.RR)

tab.all.sub <- subset(tab.all, assump %in% keep.name)
tab.all.sub$assump <- factor(tab.all.sub$assump,
                             levels = keep.name)
tab.all.sub <- tab.all.sub[order(tab.all.sub$assump), ]
tab.all.sub.1 <- data.frame(group = c("log",
                                      rep("logl_like", 6),
                                      rep("p412bar3", 2),
                                      rep("RR_1bar_2.5", 2),
                                      rep("RR_2bar_2.5", 2)),
                            tab.all.sub)

xlab.vec = c("logl" = "log-linear",
             "p4123bar_1" = 
               expression(paste("p4|12", bar(3), 
                                "/", psi, "=1",sep = "")),
             "p4123bar_1_unce" = 
               expression(paste("p4|12", bar(3), 
                                "/", psi, "~Unif(0.8,1.2)",sep = "")),
             "p41bar23bar_1" = 
               expression(paste("p4|", bar(1), "2", bar(3),
                                "/", psi, "=1",sep = "")),
             "p41bar23bar_1_unce" = 
               expression(paste("p4|", bar(1), "2", bar(3),
                                "/", psi, "~Unif(0.8,1.2)", sep = "")),
             "p412bar3bar_1" = 
               expression(paste("p4|1", bar(2), bar(3),
                                "/", psi, "=1", sep = "")),
             "p412bar3bar_1_unce" = 
               expression(paste("p4|1", bar(2), bar(3),
                                "/", psi, "~Unif(0.8,1.2)", sep = "")),
             "p412bar3_2" = 
               expression(paste("p4|1", bar(2), "3",
                                "/", psi, "=2", sep = "")),
             "p412bar3_2_unce" = 
               expression(paste("p4|1", bar(2), "3",
                                "/", psi, "~Unif(1.6,2.4)", sep = "")),
             "RR_1bar_2.5" = "RR_1bar=2.5",
             "RR_1bar_2.5_unce" = "RR_1bar~Unif(2.25,2.75)",
             "RR_2bar_2.5" = "RR_2bar=2.5",
             "RR_2bar_2.5_unce" = "RR_2bar~Unif(2.25,2.75)")

pl.4s.color <- ggplot(tab.all.sub.1, 
                      aes(x = N, y = assump, xmin = lci, xmax = uci,
                          color = group)) + 
  geom_pointrange(position = position_dodge(width=0.40)) + 
  geom_errorbar(position = position_dodge(width=0.40), 
                width = 0.1, size = 0.8) + 
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

### save results
if(!file.exists("./Results")){
  dir.create("Results")
}
ggsave(filename = paste0("./Results/", "Fig_HIV4S.png"),
       height = 7, width = 9, units = "in", dpi = 600,
       plot = pl.4s.color)








