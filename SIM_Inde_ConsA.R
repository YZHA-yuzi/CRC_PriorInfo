########################################################################
# --------------- R codes used in the simulation study --------------- #
# * Focuse on data generated under the independence assumption and 
# the constraint psi = p31bar2 
# (1) Fit the alternative model under the independence assumption
# (2) Fit the alternative model under the constraint psi = p31bar2 
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

## sce_index = 1, ..., 4 
## (index of simulation scenarios presented in Table 3)
## B_index = 1, ..., 1000/A, 
## (where A is the number of simulations ran in one R job)
## e.g., When set A = 100, the total of 1,000 simulation can be
## ran by implementing parallel computation. 
## As a result, 10 jobs are submitted to run simultaneously.  

sce.setup.ind <- read.csv(paste0("./Data/Setup_inde.csv"))

sce.setup.p31bar2.0 <- read.csv(paste0("./Data/Setup_p31bar2.csv"))
sce.setup.p31bar2 <- subset(sce.setup.p31bar2.0, psi == 0.1)

## read in R functions ##
source(paste0("FUNs.R"))

## read in data generated under 
## (1). the independence assumption and 
## (2). the constraint p31bar2 = psi ##
load(paste0("./Data/Data_inde_S", sce.setup.ind$Sce[sce_index], ".rda")) 
dat.sim.inde <- dat.sim
load(paste0("./Data/Data_p31bar2_S", sce.setup.p31bar2$Sce[sce_index], ".rda")) 
dat.sim.p31bar2 <- dat.sim

########################## BEGIN SIMULATION ###############################
set.seed(34567) 
# the number of simulations ran in one R job
nsims.curr = 10

# ---------- I. Initialize Vectors/Matrices for Storing Results --------- #
Nhat.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = nsims.curr))
est.names <- c("true", "mis")
colnames(Nhat.mat) <- paste0("N_", c(est.names))
Nhat.list <- list("inde" = Nhat.mat, "p31bar2" = Nhat.mat)

CI.list <- replicate(length(est.names), 
                     matrix(NA, ncol = 2, nrow = nsims.curr), simplify = F)
names(CI.list) <- c("true", "mis")
CI.list.all <- list("inde" = CI.list, "p31bar2" = CI.list)
nc.vec <- rep(NA, nsims.curr)
# ------ I. END Initialize Vectors/Matrices for Storing Results ----------- #


# -------------------- II. BEGIN ESTIMATION ---------------- #
B.vec = (1+nsims.curr*(B_index-1)):(nsims.curr*B_index)
start = proc.time()[3]
for(i in 1:nsims.curr){
  
  k = B.vec[i]
  
  ## Loop over different assumptions that used for generating data
  ## l = 1: the independence assumption
  ## l = 2: psi = p31bar2 
  assump.vec <- c("inde", "p31bar2")
  for(l in 1:2){
    if(l == 1){
      df_aggre <- dat.sim.inde[k, ][,-8]
    }else if(l == 2){
      df_aggre <- dat.sim.p31bar2[k, ][,-8]
    }
    ## (1). Fit the alternative framework under the independence assumption
    initial_parm_inde = c(p1 = 0.5, p2 = 0.5, p3 = 0.5)
    fit.alter.inde <- mle2(ll_inde, start = initial_parm_inde, 
                           data = list(df_aggre = df_aggre),
                           method = "L-BFGS-B", 
                           lower = c(p1 = 1e-6, p2 = 1e-6, p3 = 1e-6),
                           upper = c(p1 = 1-1e-6, p2 = 1-1e-6, p3 = 1-1e-6),
                           skip.hessian = T)
    parm.hat.vec <- coef(fit.alter.inde)
    psi.hat <- parm.hat.vec["p3"]
    if(l == 1){
      Nhat.list[[assump.vec[l]]][i, "N_true"] <- 
        as.numeric(ceiling(sum(df_aggre[1,1:6]) + 
                             c(df_aggre[1,7])/psi.hat))
    }else if(l == 2){
      Nhat.list[[assump.vec[l]]][i, "N_mis"] <- 
        as.numeric(ceiling(sum(df_aggre[1,1:6]) + 
                             c(df_aggre[1,7])/psi.hat))
    }
    ### Dirichlet-based interval
    res.CI.inde <- Dir_CI_inde(df_aggre = df_aggre, 
                               n.post = 1000, 
                               ll = ll_inde, initial = initial_parm_inde)
    nc.curr <- sum(df_aggre)
    nc.vec[i] <- nc.curr
    res.CI.inde$CI[1] <- max(nc.curr, res.CI.inde$CI[1])
    if(l == 1){
      CI.list.all[[assump.vec[l]]][["true"]][i, ] <- ceiling(res.CI.inde$CI)
    }else if(l == 2){
      CI.list.all[[assump.vec[l]]][["mis"]][i, ] <- ceiling(res.CI.inde$CI)
    }
    
    ## (2). Fit the alternative framework under the constraint p31bar2 = psi
    ### NOTE: skip fitting the proposed model under the constraint p31bar2 = psi
    ### for data generated under the constraint p31bar2 = psi, 
    ### since we would have the results by running the R script SIM_ConsA.R
    if(l == 1){
      re.bc.phi3 <- get_Nhatphi3_BC_S3(dat = df_aggre, phi = 1)
      Nhat.list[[assump.vec[l]]][i, "N_mis"] <- ceiling(re.bc.phi3$BC2)
      ### Dirichlet-based interval
      res.CI.bc2 <- Dir_CI_phi_S3(df_aggre = df_aggre, n.post = 1000, 
                                  bias.correction = "BC2", 
                                  phi.type = "phi3", phi = 1)
      res.CI.bc2$CI[1] <- max(nc.curr, res.CI.bc2$CI[1])
      CI.list.all[[assump.vec[l]]][["mis"]][i, ] <- ceiling(res.CI.bc2$CI)
    }
  } # end foor loop 
  if(i < 10){cat(k, ",")}
  if(k %% 100 == 0){cat(k, ",")}
}
print(paste0( (proc.time()[3] - start)/60, "min"))
# ---------- save results ---------- #
results = list(Nhat.mat = Nhat.list,
               CI.list = CI.list.all,
               sce.inde = sce.setup.ind[sce_index, ],
               sce.p31bar2 = sce.setup.p31bar2[sce_index, ],
               nc.vec = nc.vec)

### save results
if(!file.exists("./Results")){
  dir.create("Results")
}
save(results, 
     file = paste0("./Results/Res_inde_p31bar2_S", 
                   sce.setup.ind$Sce[sce_index], "_B", B_index, ".rda"))
############################################################################
