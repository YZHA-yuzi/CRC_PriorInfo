########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# Generating data used in simulation studies under various constraints
# Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" is located.
########################################################################

# Load packages 
library(bbmle); library(MCMCpack)
library(numDeriv); library(doParallel); library(Rcapture)

# source self-defined R functions 
source("FUNs.R")

# create a folder named "Data" under the working directory to store 
# simulated datasets 
dir.create("Data")

# 1. DATA GENERATION 
## 1.1 generate data under constraint (A) p3|1bar2 = \psi 
### set random seed 
set.seed(12345)

### explore different combinations of N and 
### psi = p31bar2 (different values of psi lead to different pc, 
### where pc is the probability of being identified at least once) 
N.vec = c(500, 1000, 2000, 5000)
p31bar2.vec = c(0.1, 0.3)
### specify other conditional probabilities 
### (i.e., p1, p21, p21bar, p312, p312bar)
sce.setup <- expand.grid(N.vec, p31bar2.vec)
sce.setup <- sce.setup[order(sce.setup$Var1), ]
sce.setup <- data.frame(S = 1:nrow(sce.setup), sce.setup)
colnames(sce.setup) <- c("Sce", "N", "p31bar2")
sce.setup$p1 = 0.35
sce.setup$p21 = 0.4
sce.setup$p21bar = 0.3
sce.setup$p312 = 0.3
sce.setup$p312bar = 0.2
sce.setup$pc = 1 - (1-sce.setup$p1)*(1-sce.setup$p21bar)*(1-sce.setup$p31bar2)
sce.setup$psi = p31bar2.vec
### BEGIN generate data based on the multinomial distribution model in eqn. (1)
for(i in 1:nrow(sce.setup)){
  parm.names <- c("p1", "p21", "p21bar", 
                  "p312", "p312bar", "p31bar2", "p31bar2bar")
  parm.vec <- c(sce.setup$p1[i], sce.setup$p21[i], sce.setup$p21bar, 
                sce.setup$p312[i], sce.setup$p312bar[i], 
                sce.setup$p31bar2[i], sce.setup$p31bar2[i])
  names(parm.vec) <- parm.names
  prob.vec <- comp.prob(x = parm.vec)
  N = sce.setup$N[i]
  print(paste0("N = ", N, ", p31bar2 = ", sce.setup$p31bar2[i]))
  nsims = 1000
  dat.sim <- gen.data(prob.vec = prob.vec, N = N, nsims = nsims)
  save(dat.sim, file = paste0("./Data/Data_p31bar2_S", i, ".rda"))
}
write.csv(sce.setup, file = paste0("./Data/Setup_p31bar2.csv"), row.names = F)


## 1.2 generate data under the independence assumption
## p1, p21 = p21bar; p312 = p312bar = p31bar2 = p31bar2bar
### set random seed 
set.seed(23456)
N.vec = c(500, 1000, 2000, 5000)
p3.vec = c(0.1)

sce.setup <- expand.grid(N.vec, p3.vec)
sce.setup <- sce.setup[order(sce.setup$Var1), ]
sce.setup <- data.frame(S = 1:nrow(sce.setup), sce.setup)
colnames(sce.setup) <- c("Sce", "N", "p3")
sce.setup$p1 = 0.35
sce.setup$p2 = 0.25
sce.setup$pc = 1 - (1-sce.setup$p1)*(1-sce.setup$p2)*(1-sce.setup$p3)
for(i in 1:nrow(sce.setup)){
  parm.names <- c("p1", "p21", "p21bar", 
                  "p312", "p312bar", "p31bar2", "p31bar2bar")
  parm.vec <- c(sce.setup$p1[i], rep(sce.setup$p2[i], 2), 
                rep(sce.setup$p3[i], 4))
  names(parm.vec) <- parm.names
  prob.vec <- comp.prob(x = parm.vec)
  N = sce.setup$N[i]
  print(paste0("N = ", N, ", p3 = ", sce.setup$p3[i]))
  nsims = 1000
  dat.sim <- gen.data(prob.vec = prob.vec, 
                      N = N, nsims = nsims)
  save(dat.sim, file = paste0("./Data/Data_inde_S", i, ".rda"))
}
write.csv(sce.setup, file = paste0("./Data/Setup_inde.csv"), row.names = F)


## 1.3 generate data under simulation scenarios presented in Table 4
## p312/p31bar2 = r*p312bar/psi
## 1.3 (a) generated data under while assuming r ~ Unif(0.75, 1.25)
N.vec = c(1000, 2000, 5000, 10000, 20000)
dist.assump = c("Uniform")
sce.setup <- expand.grid(N.vec, dist.assump)
sce.setup <- sce.setup[order(sce.setup$Var2), ]
sce.setup$Var2 <- as.character(sce.setup$Var2)
sce.setup <- data.frame(S = 1:nrow(sce.setup), sce.setup)
colnames(sce.setup) <- c("Sce", "N", "Dist_assumption")
sce.setup <- data.frame(sce.setup,
                        p1 = 0.2, p21 = 0.4, p21bar = 0.2,
                        p312 = 0.25, p31bar2 = 0.2, 
                        p312bar = 0.2,
                        psi.avg = NA, pc.avg = NA)
for(i in 1:nrow(sce.setup)){
  nsims = 1000
  N = sce.setup$N[i]
  assump.i = sce.setup$Dist_assumption[i]
  print(paste0("N = ", N, ", Dist Assumption = ", 
               sce.setup$Dist_assumption[i]))
  if(assump.i == "Uniform"){
    r.vals.vec <- runif(nsims, 0.75, 1.25)
  }
  dat.sim <- NULL
  pc.vec <- c(); psi.vec <- c()
  for(j in 1:nsims){
    parm.names <- c("p1", "p21", "p21bar", 
                    "p312", "p312bar", "p31bar2", "p31bar2bar")
    parm.vec <- c(sce.setup$p1[i], sce.setup$p21[i], sce.setup$p21bar[i],
                  get.psi.RR(p312 = sce.setup$p312[i], 
                             p31bar2 = sce.setup$p31bar2[i], 
                             p312bar = sce.setup$p312bar[i],
                             ratio = r.vals.vec[j]))
    names(parm.vec) <- parm.names
    psi.vec <- c(psi.vec, parm.vec["p31bar2bar"])
    prob.vec <- comp.prob(x = parm.vec)
    pc.vec <- c(pc.vec, sum(prob.vec[-8]))
    dat.sim.i <- gen.data(prob.vec = prob.vec, 
                          N = N, nsims = 1)
    dat.sim <- rbind.data.frame(dat.sim, dat.sim.i)
  }
  
  dat.sim.all <- list(set = data.frame(psi = psi.vec,
                                       r.val = r.vals.vec,
                                       pc = pc.vec),
                      dat = dat.sim)
  
  sce.setup$pc.avg[i] <- mean(pc.vec)
  sce.setup$psi.avg[i] <- mean(psi.vec)
  
  save(dat.sim.all, file = paste0("./Data/Data_RRunce_fixN_S", i, ".rda"))
}
write.csv(sce.setup, file = paste0("./Data/Setup_RRunce_fixN.csv"),
          row.names = F)

## 1.3 (b) generated data under while assuming r = 0.9 and r = 1.1
set.seed(24567)
N.vec = c(1000, 2000, 5000, 10000, 20000)
ratio.vec = c(0.9, 1.1)
sce.setup <- expand.grid(N.vec, ratio.vec)
sce.setup <- sce.setup[order(sce.setup$Var2), ]
sce.setup <- data.frame(S = 1:nrow(sce.setup), sce.setup)
colnames(sce.setup) <- c("Sce", "N", "RR")
sce.setup.parm <- NULL
pc.vec <- c()
prob.mat <- NULL
for(i in 1:nrow(sce.setup)){
  parm.names <- c("p1", "p21", "p21bar", 
                  "p312", "p312bar", "p31bar2", "p31bar2bar")
  
  rr = sce.setup$RR[i]
  p1 = 0.35; p21 = 0.4; p21bar = 0.2
  p312 = 0.25; p312bar = 0.25; p31bar2 = 0.2
  psi = rr*(p312bar*p31bar2/p312)
  
  parm.vec <- c(p1, p21, p21bar, p312, p312bar, p31bar2, psi)
  names(parm.vec) <- parm.names
  sce.setup.parm <- rbind(sce.setup.parm, parm.vec)
  
  prob.vec <- comp.prob(x = parm.vec)
  pc.vec <- c(pc.vec, sum(prob.vec[-8]))
  prob.mat <- rbind(prob.mat, prob.vec)
  N = sce.setup$N[i]
  print(paste0("N = ", N, ", RR = ", sce.setup$RR[i]))
  nsims = 1000
  dat.sim <- gen.data(prob.vec = prob.vec, 
                      N = N, nsims = nsims)
  save(dat.sim, file = paste0("./Data/Data_RRratio_forunce_S", i, ".rda"))
}
sce.setup.parm <- as.data.frame(sce.setup.parm)
colnames(sce.setup.parm) <- c("p1", "p21", "p21bar", 
                              "p312", "p312bar", "p31bar2", "p31bar2bar")
sce.setup <- cbind.data.frame(sce.setup,
                              sce.setup.parm)
sce.setup$pc <- pc.vec
write.csv(sce.setup, file = paste0("./Data/Setup_RRratio_forunce.csv"),
          row.names = F)



## 1.4 generate data under the OR-type constraint
## (p312/(1-p312))/(p31bar2/(1-p31bar2)) = (p312bar/(1-p312bar))/(psi/(1-psi))
## simulation results are given in Table S1 
# ------- I. DATA GENERATION under OR-type assumption -------- #
set.seed(23456)
get.psi.paircorr <- function(p312, p31bar2, or){
  inter.p312bar = (p312/(1-p312))/or
  p312bar = inter.p312bar/(1 + inter.p312bar)
  inter.psi = (p31bar2/(1-p31bar2))/or
  psi = inter.psi/(1 + inter.psi)
  return(c(p312 = p312, p312bar = p312bar, 
           p31bar2 = p31bar2, p31bar2bar = psi))
}

N.vec = c(500, 1000, 2000, 5000)
or.vec = c(0.8, 2)
sce.setup <- expand.grid(N.vec, or.vec)
sce.setup <- sce.setup[order(sce.setup$Var1), ]
sce.setup <- data.frame(S = 1:nrow(sce.setup), sce.setup)
colnames(sce.setup) <- c("Sce", "N", "OR")

sce.setup.parm <- NULL
pc.vec <- c()
prob.mat <- NULL
for(i in 1:nrow(sce.setup)){
  parm.names <- c("p1", "p21", "p21bar", 
                  "p312", "p312bar", "p31bar2", "p31bar2bar")
  parm.vec <- c(0.35, 0.4, 0.2,
                get.psi.paircorr(p312 = 0.25, p31bar2 = 0.25, 
                                 or = sce.setup$OR[i]))
  names(parm.vec) <- parm.names
  sce.setup.parm <- rbind(sce.setup.parm, parm.vec)
  
  prob.vec <- comp.prob(x = parm.vec)
  pc.vec <- c(pc.vec, sum(prob.vec[-8]))
  prob.mat <- rbind(prob.mat, prob.vec)
  N = sce.setup$N[i]
  print(paste0("N = ", N, ", OR = ", sce.setup$OR[i]))
  nsims = 1000
  dat.sim <- gen.data(prob.vec = prob.vec, N = N, nsims = nsims)
  save(dat.sim, file = paste0("./Data/Data_OR_S", i, ".rda"))
}
sce.setup.parm <- as.data.frame(sce.setup.parm)
colnames(sce.setup.parm) <- c("p1", "p21", "p21bar", 
                              "p312", "p312bar", "p31bar2", "p31bar2bar")
sce.setup <- cbind.data.frame(sce.setup,
                              sce.setup.parm)
sce.setup$pc <- pc.vec
write.csv(sce.setup, file = paste0("./Data/Setup_OR.csv"),
          row.names = F)


## 1.5 generate data under the referral scenarios
# (1) p312 = p31bar2; psi = p312bar and 
#     referral of a proportion q individuals from S1 to S3
# (2) NO corresponding log-linear models
## simulation results are given in Table S2 

set.seed(12345)
## A function to compute capture probabilities under referral scenario #####
comp.prob.referral <- function(parm.vec){
  p1 = parm.vec["p1"]
  p21 = parm.vec["p21"]
  p21bar = parm.vec["p21bar"]
  p312 = parm.vec["p312"]
  p312bar = parm.vec["p312bar"]
  p31bar2 = parm.vec["p31bar2"]
  psi = parm.vec["psi"]
  q = parm.vec["q"]
  
  p111 = p1*p21*(p312 + q*(1-p312))
  p110 = p1*p21*(1-p312)*(1-q)
  p101 = p1*(1-p21)*(p312bar + q*(1-p312bar))
  p100 = p1*(1-p21)*(1-p312bar)*(1-q)
  p011 = (1-p1)*p21bar*p31bar2
  p010 = (1-p1)*p21bar*(1-p31bar2)
  p001 = (1-p1)*(1-p21bar)*psi
  p000 = (1-p1)*(1-p21bar)*(1-psi)
  return(c(p111, p110, p101, p100, p011, p010, p001, p000))
}

N.vec = c(500, 1000, 2000, 5000)
q.vec = c(0.1, 0.3)
sce.setup <- expand.grid(N.vec, q.vec)
sce.setup <- sce.setup[order(sce.setup$Var1), ]
sce.setup <- data.frame(S = 1:nrow(sce.setup), sce.setup)
colnames(sce.setup) <- c("Sce", "N", "q")

sce.setup.parm <- NULL
pc.vec <- c()
prob.mat <- NULL
for(i in 1:nrow(sce.setup)){
  
  p1 = 0.35 
  p21 = 0.4
  p21bar = 0.3
  p312 = 0.1 
  p312bar = 0.3
  p31bar2 = p312
  psi = p312bar
  q = sce.setup[i, "q"]
  
  parm.vec <- c(p1, p21, p21bar, p312, p312bar, p31bar2, psi, q)
  parm.names <- c("p1", "p21", "p21bar", 
                  "p312", "p312bar", "p31bar2", "psi", "q")
  names(parm.vec) <- parm.names
  sce.setup.parm <- rbind(sce.setup.parm, parm.vec[-8])
  
  prob.vec <- comp.prob.referral(parm.vec = parm.vec)
  pc.vec <- c(pc.vec, sum(prob.vec[-8]))
  prob.mat <- rbind(prob.mat, prob.vec)
  N = sce.setup$N[i]
  print(paste0("N = ", N, ", q = ", sce.setup$q[i]))
  
  nsims = 1000
  dat.sim <- gen.data(prob.vec = prob.vec, N = N, nsims = nsims)
  save(dat.sim, file = paste0("./Data/Data_referral_S", i, ".rda"))
}
sce.setup.parm <- as.data.frame(sce.setup.parm)
colnames(sce.setup.parm) <- c("p1", "p21", "p21bar", 
                              "p312", "p312bar", "p31bar2", "p31bar2bar")
sce.setup <- cbind.data.frame(sce.setup,
                              sce.setup.parm)
sce.setup$pc <- pc.vec
write.csv(sce.setup, file = paste0("./Data/Setup_referral.csv"), row.names = F)




