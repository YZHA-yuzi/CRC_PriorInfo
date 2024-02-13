########################################################################
# --------------- R codes used in the simulation study --------------- #
# * Data were generated under constraint (A) p31bar2 = psi
# (1) Fit the true model
# (2) Fit the log-linear model which imposes the constraint (A) 
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

sce.setup <- read.csv(paste0("./Data/Setup_p31bar2.csv"))

## read in R functions ##
source(paste0("FUNs.R"))

## read in simulated data ##
load(paste0("./Data/Data_p31bar2_S", sce.setup$Sce[sce_index], ".rda")) 


########################## BEGIN SIMULATION ###############################
set.seed(34567) 
# the number of simulations ran in one R job
nsims.curr = 10

# ---------- I. Initialize Vectors/Matrices for Storing Results --------- #
Nhat.mat <- as.data.frame(matrix(NA, ncol = 4, nrow = nsims.curr))
est.names <- c("logl", paste0("true", c("", "_bc2", "_bcB")))
colnames(Nhat.mat) <- paste0("N_", c(est.names))

parm.hat.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = nsims.curr))
colnames(parm.hat.mat) <- c("psi", "psi_bcB")

CI.list <- replicate(length(est.names) + 1, 
                     matrix(NA, ncol = 2, nrow = nsims.curr), simplify = F)
names(CI.list) <- c(paste0(est.names[1], c("_exp", "_prof")),
                    est.names[-1])
nc.vec <- rep(NA, nsims.curr)
profCI.vec <- rep(NA, nsims.curr)
# ------ I. END Initialize Vectors/Matrices for Storing Results ----------- #


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
  get.psi.beta <- function(df_aggre, a = 0, b = 0, r = 1){
    p31bar2.hat <- (a + df_aggre$n011)/(a + b + df_aggre$n011 + df_aggre$n010)
    return(p31bar2.hat/r)
  }
  # get point estimates 
  parm.hat.mat[i, "psi"] <- 
    get.psi.beta(df_aggre = df_aggre)
  parm.hat.mat[i, "psi_bcB"] <- 
    get.psi.beta(df_aggre = df_aggre, a = 1, b = 0)

  Nhat.mat[i, "N_true"] <- 
    ceiling(get_Nhatpsi_S3(dat = df_aggre, psi = parm.hat.mat$psi[i]))
  Nhat.mat[i, "N_true_bcB"] <- 
    ceiling(get_Nhatpsi_S3(dat = df_aggre, psi = parm.hat.mat$psi_bcB[i]))
  re.bc.phi3 <- get_Nhatphi3_BC_S3(dat = df_aggre, phi = 1)
  Nhat.mat[i, "N_true_bc2"] <- ceiling(re.bc.phi3$BC2)
  
  ## (1). get Dirichlet-based intervals #
  nc.curr <- sum(df_aggre)
  nc.vec[i] <- nc.curr
  bc.type.vec = c("NO", "BC2")
  ci.names.list <- list(paste0("true", c("", "_bc2")))
  val.forloop <- c(1)
  for(lll in 1){
    phi.lll = val.forloop[lll]
    count.ci = 1
    ci.names = ci.names.list[[lll]]
    for(lb in 1:length(bc.type.vec)){
      res.CI.l <- Dir_CI_phi_S3(df_aggre = df_aggre, n.post = 1000, 
                                bias.correction = bc.type.vec[lb], 
                                phi.type = "phi3", phi = phi.lll)
      res.CI.l$CI[1] <- max(nc.curr, res.CI.l$CI[1])
      CI.list[[ci.names[count.ci]]][i, ] <- ceiling(res.CI.l$CI)
      count.ci = count.ci + 1
    }
  }
  # ad-hoc bias-correction, beta(1, 0) prior
  a.vec = c(1); b.vec = c(0)
  names.vec = c("B")
  val.forloop <- c(1)
  ci.names <- c("true")
  for(lll in 1:length(ci.names)){
    phi.lll = val.forloop[lll]
    for(ll in 1:1){
      res.CI.alter <- Dir_CI_phi_S3_beta(df_aggre = df_aggre, 
                                         n.post = 1000,
                                         phi.type = "phi3",
                                         phi = phi.lll,
                                         a = a.vec[ll], b = b.vec[ll])
      res.CI.alter$CI[1] <- max(nc.curr, res.CI.alter$CI[1])
      CI.list[[paste0(ci.names[lll], "_bc", 
                      names.vec[ll])]][i, ] <- round(res.CI.alter$CI)
    }
  }
  
  
  ## (2). fit log-linear model under the true constraint psi = p31bar2 
  fit.logl <- glm(count ~ X1 + X2 + X3 + X1X2 + X1X3 + X1X2X3, 
                  data = df_logl, family = "poisson")
  sum.fit <- summary(fit.logl)$coefficients
  se.int <- as.numeric(sum.fit[1,2])
  est.int <- as.numeric(sum.fit[1,1])
  Nhat.mat[i, "N_logl"] <- as.numeric(ceiling(nc.curr + exp(est.int)))
  ### compute Wald-type interval 
  CI.list[["logl_exp"]][i, ] <- as.numeric(nc.curr + 
                                             c(ceiling(exp(est.int - 1.96*se.int)),
                                               ceiling(exp(est.int + 1.96*se.int))))
  CI.list[["logl_exp"]][i, 1] <- max(nc.curr, CI.list[["logl_exp"]][i, 1])
  
  ## compute profile likelihood-based interval using the R function
  ## "profileCI" in the R package "Rcapture"
  dat.logl.profCI <- do.call(rbind.data.frame, 
                             apply(df_logl[,c(paste0("X", 1:3), "count")], 
                                   1, function(x) 
                                     matrix(rep(x[1:3], x[4]), ncol = 3, byrow = T)))
  colnames(dat.logl.profCI) <- paste0("c", 1:3)
  
  cond.i = ifelse(Nhat.mat[i, "N_logl"] - nc.curr <= 2, FALSE, TRUE)
  skip_to_next <- FALSE
  if(is.finite(AIC(fit.logl)) & cond.i){
    
    tryCatch(re.logl.CI <- Rcapture::profileCI(dat.logl.profCI, 
                                               mX = model.matrix(fit.logl)[,-1]),
             error = function(e) { skip_to_next <<- TRUE})
    
  }else{
    CI.list[["logl_prof"]][i, ] <- c(NA, NA)
    skip_to_next <- TRUE
  }
  
  if(skip_to_next == FALSE){
    CI.list[["logl_prof"]][i, ] <- 
      as.numeric(round(re.logl.CI$results[2:3]))
    CI.list[["logl_prof"]][i, 1] <- 
      max(CI.list[["logl_prof"]][i, 1], nc.curr)
  }
  profCI.vec[i] <- !skip_to_next
  
  if(i < 10){cat(k, ",")}
  if(k %% 100 == 0){cat(k, ",")}
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
     file = paste0("./Results/Res_p31bar2_S", sce.setup$Sce[sce_index],
                   "_B", B_index, ".rda"))
############################################################################
