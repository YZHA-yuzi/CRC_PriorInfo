########################################################################
# --------------- R codes used in the simulation study --------------- #
# * Data were generated under the OR-type constraint 
# (p312/(1-p312))/(p31bar2/(1-p31bar2)) = (p312bar/(1-p312bar))/(psi/(1-psi))
# (1) Fit the true model
# (2) Fit the log-linear model which imposes the OR-type constraint 
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

sce.setup <- read.csv(paste0("./Data/Setup_OR.csv"))

## read in R functions ##
source(paste0("FUNs.R"))

## read in simulated data ##
load(paste0("./Data/Data_OR_S", sce.setup$Sce[sce_index], ".rda")) 

########################## BEGIN SIMULATION ###############################
set.seed(34567) 
# the number of simulations ran in one R job
nsims.curr = 10

# ---------- I. Initialize Vectors/Matrices for Storing Results --------- #
Nhat.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = nsims.curr))
est.names <- c("logl", paste0("true", c("", "_bcB")))
colnames(Nhat.mat) <- paste0("N_", c(est.names))


parm.hat.mat <- as.data.frame(matrix(NA, ncol = 4*2, nrow = nsims.curr))
colnames(parm.hat.mat) <- c(c("p312", "p312bar", "p31bar2", "psi"),
                            paste0(c("p312", "p312bar", "p31bar2", "psi"), 
                                   "_B"))

CI.list <- replicate(length(est.names) + 1, 
                     matrix(NA, ncol = 2, nrow = nsims.curr), simplify = F)
names(CI.list) <- c(paste0(est.names[1], c("_exp", "_prof")),
                    est.names[-1])
nc.vec <- rep(NA, nsims.curr)
profCI.vec <- rep(NA, nsims.curr)
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
  
  ## (1). Fit the true model 
  ## (1.1) Get point estimate
  ## NO bias-corrections
  p312.hat <- df_aggre$n111/(df_aggre$n111 + df_aggre$n110)
  p31bar2.hat <- df_aggre$n011/(df_aggre$n011 + df_aggre$n010)
  p312bar.hat <- df_aggre$n101/(df_aggre$n101 + df_aggre$n100)
  inter.or <- (p31bar2.hat/(1-p31bar2.hat))*
    (p312bar.hat/(1-p312bar.hat))/(p312.hat/(1-p312.hat))
  psi.hat <- inter.or/(1 + inter.or)
  
  ## Bias-corrections based on Beta(1,0) priors
  get.betapost <- function(a, b, df_aggre){
    sum.ab = a + b
    p312.post <- (a + df_aggre$n111)/(sum.ab + 
                                        df_aggre$n111 + df_aggre$n110)
    p31bar2.post <- (a + df_aggre$n011)/(sum.ab + 
                                           df_aggre$n011 + df_aggre$n010)
    p312bar.post <- (a + df_aggre$n101)/(sum.ab + 
                                           df_aggre$n101 + df_aggre$n100)
    inter.or.post <- (p31bar2.post/(1-p31bar2.post))*
      (p312bar.post/(1-p312bar.post))/(p312.post/(1-p312.post))
    psi.post <- inter.or.post/(1 + inter.or.post)
    return(c(p312.post, p312bar.post, p31bar2.post, psi.post))
  }
  parm.names <- c("p312", "p312bar", "p31bar2", "psi")
  parm.hat.mat[i, parm.names] <- c(p312.hat, p312bar.hat, p31bar2.hat, psi.hat)
  parm.hat.mat[i, paste0(parm.names, "_B")] <- 
    get.betapost(a = 1, b = 0, df_aggre = df_aggre)

  nc.rm = sum(df_aggre[1,1:6])
  Nhat.mat[i, "N_true"] <- 
    as.numeric(ceiling(nc.rm + c(df_aggre[1,7])/psi.hat))
  Nhat.mat[i, "N_true_bcB"] <- 
    as.numeric(ceiling(nc.rm +c(df_aggre[1,7])/parm.hat.mat$psi_B[i]))

  
  ## (1.2). get Dirichlet-based intervals #
  # no bias-correction 
  res.CI.alter <- Dir_CI_OR(df_aggre = df_aggre, 
                            n.post = 1000, or = 1)
  nc.curr <- sum(df_aggre)
  nc.vec[i] <- nc.curr
  res.CI.alter$CI[1] <- max(nc.curr, res.CI.alter$CI[1])
  CI.list[["true"]][i, ] <- round(res.CI.alter$CI)
  
  # bias-corrections, beta(1, 0) priors
  a.vec = c(1); b.vec = c(0)
  names.vec = c("B")
  for(ll in 1){
    res.CI.alter <- Dir_CI_OR(df_aggre = df_aggre, 
                              n.post = 1000,
                              or = 1, a = a.vec[ll], b = b.vec[ll])
    res.CI.alter$CI[1] <- max(nc.curr, res.CI.alter$CI[1])
    CI.list[[paste0("true_bc", 
                    names.vec[ll])]][i, ] <- round(res.CI.alter$CI)
  }
  
  ## (2). fit log-linear model under the true constraint  
  fit.logl <- glm(count ~ X1 + X2 + X3 + X1X2 + X1X3 + X2X3, 
                  data = df_logl, family = "poisson")
  sum.fit <- summary(fit.logl)$coefficients
  se.int <- as.numeric(sum.fit[1,2])
  est.int <- as.numeric(sum.fit[1,1])
  Nhat.mat[i, "N_logl"] <- as.numeric(ceiling(nc.curr + exp(est.int)))
  CI.list[["logl_exp"]][i, ] <- as.numeric(nc.curr + 
                                             c(ceiling(exp(est.int - 1.96*se.int)),
                                               ceiling(exp(est.int + 1.96*se.int))))
  CI.list[["logl_exp"]][i, 1] <- max(nc.curr, CI.list[["logl_exp"]][i, 1])
  
  dat.logl.profCI <- do.call(rbind.data.frame,
                             apply(df_logl[,c(paste0("X", 1:3), "count")],
                                   1, function(x)
                                     matrix(rep(x[1:3], x[4]), 
                                            ncol = 3, byrow = T)))
  colnames(dat.logl.profCI) <- paste0("c", 1:3)
  
  skip_to_next <- FALSE
  if(is.finite(AIC(fit.logl))){
    tryCatch(re.logl.CI <- profileCI(dat.logl.profCI,
                                     mX = df_logl[,c(paste0("X", 1:3),
                                                     "X1X2", "X1X3", "X2X3")]),
             error = function(e) { skip_to_next <<- TRUE})
  }else{
    CI.list[["logl_prof"]][i, ] <- c(NA, NA)
    skip_to_next <- TRUE
  }
  
  if(skip_to_next == FALSE){
    CI.list[["logl_prof"]][i, ] <- 
      as.numeric(ceiling(re.logl.CI$results[2:3]))
    CI.list[["logl_prof"]][i, 1] <- 
      max(CI.list[["logl_prof"]][i, 1], nc.curr)
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
               CI.list = CI.list,
               sce = sce.setup[sce_index, ],
               nc.vec = nc.vec,
               profCI.vec = profCI.vec)

### save results
if(!file.exists("./Results")){
  dir.create("Results")
}
save(results, 
     file = paste0("./Results/Res_OR_S", sce.setup$Sce[sce_index],
                   "_B", B_index, ".rda"))
###################################################################################################################





