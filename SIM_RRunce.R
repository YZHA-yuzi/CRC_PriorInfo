########################################################################
# --------------- R codes used in the simulation study --------------- #
# * Focuse on data generated under following constraints 
# p312/p31bar2 = r*p312bar/psi
# (a) r ~ Unif(0.75, 1.25)
# (b) r = 0.9
# (c) r = 1.1
# 
# Fit the alternative model under following constraints 
# p312/p31bar2 = r*p312bar/psi
# (1) r = 1
# (2) r ~ N(1, 0.06^2)
# (3) r ~ Unif(0.75, 1.25)
# NOTE:
# Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" is located.
########################################################################


args <-  commandArgs(trailingOnly = TRUE)
sce_index = eval( parse(text=args[1]) )
B_index = eval( parse(text=args[2]) )

## load in packages ##
library(bbmle)
library(MCMCpack)
library(numDeriv)
library(doParallel)
library(Rcapture)

## sce_index = 1, ..., 5 
## (index of simulation scenarios presented in Table 4)
## N = 1,000, 2,000, 5,000, 10,000, 20,000
## B_index = 1, ..., 1000/A, 
## (where A is the number of simulations ran in one R job)
## e.g., When set A = 100, the total of 1,000 simulation can be
## ran by implementing parallel computation. 
## As a result, 10 jobs are submitted to run simultaneously.  

sce.setup.unce <- read.csv(paste0("./Data/Setup_RRunce_fixN.csv"))
sce.setup.fix <- read.csv(paste0("./Data/Setup_RRratio_forunce.csv"))
sce.setup.fix1 <- subset(sce.setup.fix, RR == 0.9)
sce.setup.fix2 <- subset(sce.setup.fix, RR == 1.1)

## read in R functions ##
source(paste0("FUNs.R"))


## read in data generated under following constraints 
# p312/p31bar2 = r*p312bar/psi
# (a) r ~ Unif(0.75, 1.25)
# (b) r = 0.9
# (c) r = 1.1
load(paste0("./Data/Data_RRunce_fixN_S", 
            sce.setup.unce$Sce[sce_index], ".rda"))
dat.sim.unce <- dat.sim.all$dat
load(paste0("./Data/Data_RRratio_forunce_S", 
            sce.setup.fix1$Sce[sce_index], ".rda"))
dat.sim.fix1 <- dat.sim
load(paste0("./Data/Data_RRratio_forunce_S", 
            sce.setup.fix2$Sce[sce_index], ".rda"))
dat.sim.fix2 <- dat.sim

########################## BEGIN SIMULATION ###############################
set.seed(1234 + sce_index)
# the number of simulations ran in one R job
nsims.curr = 10

# ---------- I. Initialize Vectors/Matrices for Storing Results --------- #
Nhat.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = nsims.curr))
est.names <- c(paste0("N_", c("r1", "Normal", "Unif")))
colnames(Nhat.mat) <- est.names
Nhat.list <- list("unif" = Nhat.mat, "r09" = Nhat.mat, "r11" = Nhat.mat)

CI.list <- replicate(3, matrix(NA, ncol = 2, nrow = nsims.curr), simplify = F)
names(CI.list) <- est.names
CI.list.all <- list("unif" = CI.list, "r09" = CI.list, "r11" = CI.list)

nc.vec <- rep(NA, nsims.curr)
# ------ I. END Initialize Vectors/Matrices for Storing Results --------- #


# -------------------- II. BEGIN ESTIMATION ---------------- #
B.vec = (1+nsims.curr*(B_index-1)):(nsims.curr*B_index)
start = proc.time()[3]
for(i in 1:nsims.curr){
  
  k = B.vec[i]
  ## loop over different assumptions that used for generating data 
  ## "unif": Uniform(0.75, 1.25)
  ## "r09": r = 0.9
  ## "r11": r = 1.1
  assump.vec <- c("unif", "r09", "r11")
  for(l in 1:3){
    
    if(l == 1){
      df_aggre <- dat.sim.unce[k, ][,-8]
    }else if(l == 2){
      df_aggre <- dat.sim.fix1[k, ][,-8]
    }else if(l == 3){
      df_aggre <- dat.sim.fix2[k, ][,-8]
    }
    nc.curr = as.numeric(sum(df_aggre))
    nc.vec[i] <- nc.curr
    
    # (1) fit the proposed model under the RR-type constraint while assuming r = 1
    ratio.true <- 1
    get.betapost <- function(a, b, df_aggre, ratio.true){
      sum.ab = a + b
      p312.post <- (a + df_aggre$n111)/(sum.ab + 
                                          df_aggre$n111 + df_aggre$n110)
      p31bar2.post <- (a + df_aggre$n011)/(sum.ab + 
                                             df_aggre$n011 + df_aggre$n010)
      p312bar.post <- (a + df_aggre$n101)/(sum.ab + 
                                             df_aggre$n101 + df_aggre$n100)
      psi.post <- p312bar.post*p31bar2.post/p312.post
      return(c(p312.post, p312bar.post, p31bar2.post, psi.post))
    }
    nc.rm = sum(df_aggre[1,1:6])
    # implement bias-corrections relying on Beta(1,0) priors
    psi.hat.bcB <- get.betapost(a = 1, b = 0, df_aggre = df_aggre,
                                ratio.true = ratio.true)[4]
    Nhat.bcB.i <- as.numeric(ceiling(nc.rm + c(df_aggre[1,7])/psi.hat.bcB))
    Nhat.list[[assump.vec[l]]][i, "N_r1"] <- Nhat.bcB.i
      
   
    res.CI.alter <- Dir_CI_RR(df_aggre = df_aggre, 
                              n.post = 1000, rr = ratio.true, a = 1, b = 0)
    res.CI.alter$CI[1] <- max(nc.curr, res.CI.alter$CI[1])
    CI.list.all[[assump.vec[l]]][["N_r1"]][i, ] <- round(res.CI.alter$CI)
    
    # (2) uncertainty analysis r ~ N(1, 0.06^2)
    ## the point estimate is the same as assuming r = 1
    Nhat.list[[assump.vec[l]]][i, "N_Normal"] <- Nhat.bcB.i
    re.CI.dist.n <- Dir_CI_RR_unce(df_aggre = df_aggre, 
                                   n.post = 100,
                                   assume.dist = "Normal",
                                   a.dist = 1, b.dist = 0.06,
                                   a = 1, b = 0)
    re.CI.dist.n$CI[1] <- max(nc.curr, re.CI.dist.n$CI[1])
    CI.list.all[[assump.vec[l]]][["N_Normal"]][i, ] <- round(re.CI.dist.n$CI)
    
    # (3) uncertainty analysis r ~ Unif(0.75, 1.25)
    ## the point estimate is the same as assuming r = 1, since the 
    ## assumed Uniform dist is centered at 1 
    Nhat.list[[assump.vec[l]]][i, "N_Unif"] <- Nhat.bcB.i
    re.CI.dist.u <- Dir_CI_RR_unce(df_aggre = df_aggre, 
                                   n.post = 100,
                                   assume.dist = "Uniform",
                                   a.dist = 0.75, b.dist = 1.25,
                                   a = 1, b = 0)
    re.CI.dist.u$CI[1] <- max(nc.curr, re.CI.dist.u$CI[1])
    CI.list.all[[assump.vec[l]]][["N_Unif"]][i, ] <- round(re.CI.dist.u$CI)
  } # end loop over different assumptions

  if(i < 10){cat(k, ",")}
  if(k %% 100 == 0){cat(k, ",")}
  
}
print(paste0( (proc.time()[3] - start)/60, "min"))
# ---------- save results ---------- #
results = list(Nhat.mat = Nhat.list,
               CI.list = CI.list.all,
               sce.unce = sce.setup.unce[sce_index, ],
               sce.r09 = sce.setup.fix1[sce_index, ],
               sce.r11 = sce.setup.fix2[sce_index, ],
               nc.vec = nc.vec)

### save results
if(!file.exists("./Results")){
  dir.create("Results")
}
save(results, 
     file = paste0("./Results/Res_RRunce_S", 
                   sce.setup.unce$Sce[sce_index], "_B", B_index, ".rda"))
###############################################################################





