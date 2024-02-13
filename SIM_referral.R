########################################################################
# --------------- R codes used in the simulation study --------------- #
# * Data were generated under the referral scenarios
# (1) p312 = p31bar2; psi = p312bar and 
#     referral of a proportion q individuals from S1 to S3
# (2) NO corresponding log-linear models
## simulation results are given in Table S2 
# (1) Fit the true model
# (2) NO corresponding log-linear models; 
#     fit 8 possible candidate log-linear models, 
#     and select the one with the lowest AIC
# NOTE:
# Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" is located.
########################################################################

args <-  commandArgs(trailingOnly = TRUE)
sce_index = eval( parse(text=args[1]) )
B_index = eval( parse(text=args[2]) )

## Load in packages ##
library(bbmle); library(MCMCpack)
library(numDeriv); library(doParallel); library(Rcapture)

## sce_index = 1, ..., 8 
## (index of simulation scenarios under the constraint p31bar2 = psi)
## B_index = 1, ..., 1000/A, 
## (where A is the number of simulations ran in one R job)
## e.g., When set A = 100, the total of 1,000 simulation can be
## ran by implementing parallel computation. 
## As a result, 10 jobs are submitted to run simultaneously.  

sce.setup <- read.csv(paste0("./Data/Setup_referral.csv"))

## read in R functions ##
source(paste0("FUNs.R"))

## read in simulated data ##
load(paste0("./Data/Data_referral_S", sce.setup$Sce[sce_index], ".rda")) 

########################## BEGIN SIMULATION ###############################
set.seed(34567) 
# the number of simulations ran in one R job
nsims.curr = 10

# ---------- I. Initialize Vectors/Matrices for Storing Results --------- #
Nhat.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = nsims.curr))
est.names <- c("logl", "true")
colnames(Nhat.mat) <- paste0("N_", c(est.names))


parm.hat.mat <- SE.hat.mat <- 
  as.data.frame(matrix(NA, ncol = 8, nrow = nsims.curr))
colnames(parm.hat.mat) <- colnames(SE.hat.mat) <- 
  c("p1", "p21", "p21bar", "p312", "p312bar", "p31bar2", "psi", "q")

CI.list <- replicate(3, matrix(NA, ncol = 2, nrow = nsims.curr), simplify = F)
names(CI.list) <- c("logl_exp", "logl_prof", "true")

nc.vec <- rep(NA, nsims.curr)
profCI.vec <- rep(NA, nsims.curr)

sel.logl.models <- rep(NA, nsims.curr)
psi.logl.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = nsims.curr))
colnames(psi.logl.mat) <- c("psi.hat.obs", "psi.hat.fitted")
# ------ I. END Initialize Vectors/Matrices for Storing Results --------- #


# -------------------- II. BEGIN ESTIMATION ---------------- #
B.vec = (1+nsims.curr*(B_index-1)):(nsims.curr*B_index)
start = proc.time()[3]
for(i in 1:nsims.curr){
  
  k = B.vec[i]
  df_aggre <- dat.sim[k, ][,-8]
  
  # prepare data for log-linear model #
  df_logl <- gen_dflogl(dat = df_aggre, 
                        nstreams = 3, 
                        predictors = c(paste0("X", 1:3), 
                                       "X1X2", "X1X3", "X2X3", "X1X2X3"))
  
  ## (1) Fit the true model
  ## get estimated N under referral scenario
  ## using the alternative model
  fit.alter <- fit_referral(dat = df_aggre)
  
  parm.hat.i <- coef(fit.alter)
  se.hat.i <- sqrt(diag(fit.alter@vcov))
  
  psi.hat <- coef(fit.alter)["p312bar"]
  parm.hat.mat[i, ] <- as.numeric(c(parm.hat.i[1:5],
                                    parm.hat.i[4:6]))
  SE.hat.mat[i, ] <- as.numeric(c(se.hat.i[1:5],
                                  se.hat.i[4:6]))
  Nhat.mat[i, "N_true"] <- as.numeric(ceiling(sum(df_aggre[1,1:6]) + 
                                                 c(df_aggre[1,7])/psi.hat))
  
  ## get Dirichlet-based CI ##
  res.CI.alter <- Dir_CI_referral(df_aggre = df_aggre, n.post = 1000)
  nc.curr <- sum(df_aggre)
  nc.vec[i] <- nc.curr
  res.CI.alter$CI[1] <- max(nc.curr, res.CI.alter$CI[1])
  CI.list[["true"]][i, ] <- ceiling(res.CI.alter$CI)
  
  ### fit log-linear 8 possible log-linear models
  gene_input <- function(x){
    paste0("count ~ X1+X2+X3+", paste0(x, collapse = "+"))
  }
  predictors_S3 <- c("X1X2", "X1X3", "X2X3")
  models_form_all <- unlist(lapply(1:3, function(n)
    combn(predictors_S3, n, FUN = gene_input)))
  models_form_all <- c("count ~ X1+X2+X3", models_form_all)
  # fit possible log-linear models (with and without highest interaction)
  fitall_S3_p <- lapply(models_form_all,
                        function(y) glm(y, family = "poisson", data = df_logl))
  sel.logl.index <- which.min(sapply(fitall_S3_p, AIC))
  sel.logl.models[i] <- models_form_all[sel.logl.index]
  
  fit.logl <- fitall_S3_p[[sel.logl.index]]
  
  sum.fit <- summary(fit.logl)$coefficients
  se.int <- as.numeric(sum.fit[1,2])
  est.int <- as.numeric(sum.fit[1,1])
  Nhat.mat[i, "N_logl"] <- as.numeric(ceiling(nc.curr + exp(est.int)))
  CI.list[["logl_exp"]][i, ] <- as.numeric(nc.curr + 
                                             c(ceiling(exp(est.int - 1.96*se.int)),
                                               ceiling(exp(est.int + 1.96*se.int))))
  CI.list[["logl_exp"]][i, 1] <- max(nc.curr, CI.list[["logl_exp"]][i, 1])
  
  ## psihat from the AIC-selected log-linear model
  X.design0 <- model.matrix(fit.logl)
  X.design <- rbind(X.design0,
                    c(1, rep(0, ncol(X.design0)-1)))
  n.fitted <- exp(X.design %*% matrix(coef(fit.logl), ncol = 1))
  psi.hat.logl.fitted <- n.fitted[7]/(sum(n.fitted[7:8]))
  psi.hat.logl <- df_aggre[1, 7]/((df_aggre[1, 7] + n.fitted[8]))
  
  psi.logl.mat[i, ] <- as.numeric(c(psi.hat.logl, psi.hat.logl.fitted))
  ## compute profile likelihood CI
  dat.logl.profCI <- do.call(rbind.data.frame,
                             apply(df_logl[,c(paste0("X", 1:3), "count")],
                                   1, function(x)
                                     matrix(rep(x[1:3], x[4]), 
                                            ncol = 3, byrow = T)))
  colnames(dat.logl.profCI) <- paste0("c", 1:3)
  skip_to_next <- FALSE
  if(is.finite(AIC(fit.logl))){
    tryCatch(re.logl.CI <- profileCI(dat.logl.profCI,
                                     mX = model.matrix(fit.logl)[,-1]),
             error = function(e) { skip_to_next <<- TRUE})
  }else{
    CI.list[["logl_prof"]][i, ] <- c(NA, NA)
    skip_to_next <- TRUE
  }
  
  if(skip_to_next == FALSE){
    CI.list[["logl_prof"]][i, ] <- as.numeric(round(re.logl.CI$results[2:3]))
    CI.list[["logl_prof"]][i, 1] <- max(CI.list[["logl_prof"]][i, 1], nc.curr)
  }
  
  profCI.vec[i] <- !skip_to_next
  
  if(i < 10){cat(k, ",")}
  if(k %% 100 == 0){cat(k, ",")}
  # cat(i, "\r")
  # print((proc.time()[3] - start))
  # print(paste0( (proc.time()[3] - start)/60, "min"))
  
}
print(paste0( (proc.time()[3] - start)/60, "min"))
# ---------- save results ---------- #
results = list(Nhat.mat = Nhat.mat,
               parm.hat.mat = parm.hat.mat,
               SE.hat.mat = SE.hat.mat,
               CI.list = CI.list,
               sel.loglmodels = sel.logl.models,
               psi.logl.mat = psi.logl.mat,
               sce = sce.setup[sce_index, ],
               nc.vec = nc.vec,
               profCI.vec = profCI.vec)

### save results
if(!file.exists("./Results")){
  dir.create("Results")
}
save(results, 
     file = paste0("./Results/Res_referral_S", sce.setup$Sce[sce_index],
                   "_B", B_index, ".rda"))
#############################################################################





