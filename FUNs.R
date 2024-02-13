############################################################################
#                      SELF-DEFINED FUNCTIONS 
############################################################################
# A function to compute prob of each capture profile based on parameters vector
# p1, p21, p21bar, p312, p312bar, p31bar2, psi 
comp.prob <- function(x){
  p1 = x[1]; p21 = x[2]; p21bar = x[3];
  p312 = x[4]; p312bar = x[5]; p31bar2 = x[6]; psi = x[7]
  
  p111 = p1*p21*p312; p110 = p1*p21*(1-p312)
  p101 = p1*(1-p21)*p312bar; p100 = p1*(1-p21)*(1-p312bar)
  p011 = (1-p1)*p21bar*p31bar2; p010 = (1-p1)*p21bar*(1-p31bar2)
  p001 = (1-p1)*(1-p21bar)*psi; p000 = (1-p1)*(1-p21bar)*(1-psi)
  prob.vec = c(p111, p110, p101, p100, p011, p010, p001, p000)
  return(prob.vec)
}

# A function to generate data from multinomial distribution
# the parameters are probabilities of each capture histories 
gen.data <- function(prob.vec, N, nsims){
  # prob.vec: p111, p110, ..., p000
  dat.gen <- rmultinom(nsims, size = N, prob = prob.vec)
  dat.gen <- t(dat.gen)
  dat.gen <- as.data.frame(dat.gen)
  colnames(dat.gen) <- c("n111", "n110", "n101", "n100", 
                         "n011", "n010", "n001", "n000")
  dat.gen[dat.gen == 0] <- 0.5
  return(dat.gen)
}

# A function to generate dataframe used in fitting log-linear models 
gen_dflogl <- function(dat, nstreams, predictors){
  # INPUTE: 
  # dat: aggregated data at the population level
  # nstreams: the number of data streams
  # predictors: name of all possible predictors 
  df = expand.grid(rep(list(c(1,0)), nstreams))[-(2^nstreams), ]
  for(i in nstreams:1){
    df = df[order(df[,i], decreasing = T), ]
  }
  inter.list <- NULL
  for(i in 2:nstreams){
    inter.list[[i-1]] <- apply(apply(combn(1:nstreams, i, FUN = paste0), 2, 
                                     as.numeric), 
                               2, function(x) apply(df[,x], 1, prod))
  }
  df <- cbind.data.frame(df, do.call(cbind.data.frame, inter.list))
  colnames(df) <- predictors
  df$count = as.numeric(dat[1, ])
  return(df)
  # OUTPUT: data frame used for fitting log-linear model with capture indicators 
}

# MLE of N with known psi = p31bar2bar 
get_Nhatpsi_S3 <- function(dat, psi){
  return(sum(dat[1,1:6]) + dat$n001/psi)
}

# phi1 = p312/psi
get_Nhatphi1_S3 <- function(dat, phi){
  return(sum(dat[1,1:6]) + dat$n001*(dat$n111+dat$n110)/dat$n111*phi)
}

# phi2 = p312bar/psi
get_Nhatphi2_S3 <- function(dat, phi){
  return(sum(dat[1,1:6]) + dat$n001*(dat$n101+dat$n100)/dat$n101*phi)
}

# phi3 = p31bar2/psi
get_Nhatphi3_S3 <- function(dat, phi){
  return(sum(dat[1,1:6]) + dat$n001*(dat$n011+dat$n010)/dat$n011*phi)
}


#### Functions to compute bias-corrected estimator for 3-S CRC data
## phi1 = p312/psi
get_Nhatphi1_BC_S3 <- function(dat, phi){
  cor.term.1 = (dat$n110*dat$n001)/(dat$n111)^2
  cor.term.2 = (dat$n110*dat$n001)/(dat$n111 + 0.5)^2
  Nhat = get_Nhatphi1_S3(dat = dat, phi = phi)
  re = list(BC = Nhat - cor.term.1, BC2 = Nhat - cor.term.2)
  return(re)
}

## phi2 = p312bar/psi
get_Nhatphi2_BC_S3 <- function(dat, phi){
  cor.term.1 = (dat$n100*dat$n001)/(dat$n101)^2
  cor.term.2 = (dat$n100*dat$n001)/(dat$n101 + 0.5)^2
  Nhat = get_Nhatphi2_S3(dat = dat, phi = phi)
  re = list(BC = Nhat - cor.term.1, BC2 = Nhat - cor.term.2)
  return(re)
}

## phi3 = p31bar2/psi
get_Nhatphi3_BC_S3 <- function(dat, phi){
  cor.term.1 = (dat$n010*dat$n001)/(dat$n011)^2
  cor.term.2 = (dat$n010*dat$n001)/(dat$n011 + 0.5)^2
  Nhat = get_Nhatphi3_S3(dat = dat, phi = phi)
  re = list(BC = Nhat - cor.term.1, BC2 = Nhat - cor.term.2)
  return(re)
}

## A function to compute dirichlet credible interval 
## under the assumption p312/psi = phi1 ##
Dir_CI_phi_S3 <- function(df_aggre, n.post, 
                          bias.correction = "NO", 
                          phi.type = "phi1",
                          phi){
  # phi.type: phi1, phi2, phi3
  # phi1 = p312/psi
  # phi2 = p312bar/psi
  # phi3 = p31bar2/psi
  # bias.correction = "NO", "BC1", "BC2"
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  if(phi.type == "phi1"){
    f <- function(nc, parm, bias.correction, phi){
      dat_simed <- as.data.frame(t(nc*parm))
      colnames(dat_simed) <- c("n111", "n110", "n101", 
                               "n100", "n011", "n010", "n001")
      if(bias.correction == "NO"){
        Nhat <- get_Nhatphi1_S3(dat = dat_simed, phi = phi)
      }else if(bias.correction == "BC1"){
        re.inter <- get_Nhatphi1_BC_S3(dat = dat_simed, phi = phi)
        Nhat <- re.inter$BC
      }else if(bias.correction == "BC2"){
        re.inter <- get_Nhatphi1_BC_S3(dat = dat_simed, phi = phi)
        Nhat <- re.inter$BC2
      }
      pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
      pcj[pcj > 1] <- 1
      return(pcj)
    }
  }else if(phi.type == "phi2"){
    f <- function(nc, parm, bias.correction, phi){
      dat_simed <- as.data.frame(t(nc*parm))
      colnames(dat_simed) <- c("n111", "n110", "n101", 
                               "n100", "n011", "n010", "n001")
      if(bias.correction == "NO"){
        Nhat <- get_Nhatphi2_S3(dat = dat_simed, phi = phi)
      }else if(bias.correction == "BC1"){
        re.inter <- get_Nhatphi2_BC_S3(dat = dat_simed, phi = phi)
        Nhat <- re.inter$BC
      }else if(bias.correction == "BC2"){
        re.inter <- get_Nhatphi2_BC_S3(dat = dat_simed, phi = phi)
        Nhat <- re.inter$BC2
      }
      pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
      pcj[pcj > 1] <- 1
      return(pcj)
    }
  }else if(phi.type == "phi3"){
    f <- function(nc, parm, bias.correction, phi){
      dat_simed <- as.data.frame(t(nc*parm))
      colnames(dat_simed) <- c("n111", "n110", "n101", 
                               "n100", "n011", "n010", "n001")
      if(bias.correction == "NO"){
        Nhat <- get_Nhatphi3_S3(dat = dat_simed, phi = phi)
      }else if(bias.correction == "BC1"){
        re.inter <- get_Nhatphi3_BC_S3(dat = dat_simed, phi = phi)
        Nhat <- re.inter$BC
      }else if(bias.correction == "BC2"){
        re.inter <- get_Nhatphi3_BC_S3(dat = dat_simed, phi = phi)
        Nhat <- re.inter$BC2
      }
      pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
      pcj[pcj > 1] <- 1
      return(pcj)
    }
  }
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               bias.correction = bias.correction,
               phi = phi)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}

## A function to compute dirichlet credible interval 
## under the assumption p312/psi = phi1 ##
## ad-hoc bias-correction 
Dir_CI_phi_S3_beta <- function(df_aggre, n.post, 
                               a = 0, b = 0,
                               phi.type = "phi1",
                               phi){
  # phi.type: phi1, phi2, phi3
  # phi1 = p312/psi
  # phi2 = p312bar/psi
  # phi3 = p31bar2/psi
  # a, b are two parameters for beta prior, 
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  if(phi.type == "phi1"){
    f <- function(nc, parm, a, b, phi){
      dat_simed <- as.data.frame(t(nc*parm))
      colnames(dat_simed) <- c("n111", "n110", "n101", 
                               "n100", "n011", "n010", "n001")
      p312.hat <- (dat_simed$n111 + a)/(a + b + dat_simed$n111 + dat_simed$n110)
      estpsi <- p312.hat/phi
      Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
      pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
      pcj[pcj > 1] <- 1
      return(pcj)
    }
  }else if(phi.type == "phi2"){
    f <- function(nc, parm, a, b, phi){
      dat_simed <- as.data.frame(t(nc*parm))
      colnames(dat_simed) <- c("n111", "n110", "n101", 
                               "n100", "n011", "n010", "n001")
      p312bar.hat <- (dat_simed$n101 + a)/
        (a + b + dat_simed$n101 + dat_simed$n100)
      estpsi <- p312bar.hat/phi
      Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
      pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
      pcj[pcj > 1] <- 1
      return(pcj)
    }
  }else if(phi.type == "phi3"){
    f <- function(nc, parm, a, b, phi){
      dat_simed <- as.data.frame(t(nc*parm))
      colnames(dat_simed) <- c("n111", "n110", "n101", 
                               "n100", "n011", "n010", "n001")
      p31bar2.hat <- (dat_simed$n011 + a)/
        (a + b + dat_simed$n011 + dat_simed$n010)
      estpsi <- p31bar2.hat/phi
      Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
      pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
      pcj[pcj > 1] <- 1
      return(pcj)
    }
  }
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               a = a, b = b,
               phi = phi)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}



# ------ Under the independence assumption -------- # 
## (1). define log-likelihood 
ll_inde <- function(parms){
  p1 = parms[['p1']]; 
  p21 = parms[['p2']]; p21bar = parms[['p2']]; 
  p312 = parms[['p3']]; p312bar = parms[['p3']]; 
  p31bar2 = parms[['p3']]; psi = parms[['p3']]
  pc = max(1E-120, 1-(1-p1)*(1-p21bar)*(1-psi))
  
  lamb111 = max(p1*p21*p312, 1e-10)
  lamb110 = max(p1*p21*(1-p312), 1e-10)
  lamb101 = max(p1*(1-p21)*p312bar, 1e-10)
  lamb100 = max(p1*(1-p21)*(1-p312bar), 1e-10)
  lamb011 = max((1-p1)*p21bar*p31bar2, 1e-10)
  lamb010 = max((1-p1)*p21bar*(1-p31bar2), 1e-10)
  lamb001 = max((1-p1)*(1-p21bar)*psi, 1e-10)
  
  n111 = df_aggre[1, 1]; n110 = df_aggre[1, 2]
  n101 = df_aggre[1, 3]; n100 = df_aggre[1, 4]
  n011 = df_aggre[1, 5]; n010 = df_aggre[1, 6]
  n001 = df_aggre[1, 7]
  
  ll = n111*log(lamb111) + n110*log(lamb110) + n101*log(lamb101)+
    n100*log(lamb100) + n011*log(lamb011)+
    n010*log(lamb010) + n001*log(lamb001) - 
    (n111 + n110 + n101 + n100 + n011 + n010 + n001)*log(pc)
  return(-1*ll)
}
parnames(ll_inde) <- c("p1", "p2", "p3")

# (2). functions to get Dirichlet credible interval under 
# the independence assumption for three-stream cases #
# STEPS: 
# Firstly, based on those sampled conditional probs and nc 
# (total observed cell counts in dataset), generate 'simulated' observed data nobs 
# Secondly, compute MLE of N for 'simulated' observed data
# Lastly, use nobs/Nhat to estimate unconditional probs 
# Calculate pc from those unconditional probs, 
# Nhatnew = nc/pc, then sample a new nc from Binomial(Nhatnew, pc)
# Finally, posterior distribution of Nhat is formed by (new nc)/(pc) 
## (2.1). define a function used for converting generated conditional 
## probs to unconditional probs ##
## (2.2). function returns posterior dist and quantiles ##
Dir_CI_inde <- function(df_aggre, n.post, ll, initial){
  ## INPUT: 
  ## df_aggre: observed 3-stream CRC data
  ## n.post: the number of samples for estimated N
  ## ll: likelihood function defined under the independence assumption
  ## initial: initial value of parameters (p1, p2, p3)
  ### A function to obtain pc under the independence assumption
  find_uncond_est_inde <- function(nc, parm, ll, initial){
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- c("n111", "n110", "n101", 
                             "n100", "n011", "n010", "n001")
    fit <- mle2(ll, start = initial, 
                data = list(df_aggre = dat_simed), 
                method = "L-BFGS-B", 
                lower = c(p1 = 1e-6, p2 = 1e-6, p3 = 1e-6),
                upper = c(p1 = 1-1e-6, p2 = 1-1e-6, p3 = 1-1e-6),
                skip.hessian = T)
    estparm <- coef(fit)
    estpsi <- estparm["p3"]
    Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
    pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
    pcj[pcj > 1] <- 1 
    return(pcj)
  }
  
  ### do parallel computation
  no_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', unset="1"))
  registerDoParallel(cores = no_cores)
  cl <- makeCluster(no_cores, type="FORK")
  pstarcond <- rdirichlet(n.post,
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001
  pcj <- parApply(cl = cl, pstarcond, 1, find_uncond_est_inde,
                  nc = as.numeric(sum(df_aggre)),
                  ll = ll, initial = initial)
  stopCluster(cl)
  
  # pcj <- apply(pstarcond, 1, find_uncond_est_inde, 
  #              nc = as.numeric(sum(df_aggre)),
  #              ll = ll, initial = initial)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  # not condition on nc 
  ncnew <- rbinom(length(Nnew), size = Nnew, prob = pcj)
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}


## A function to compute psi given p312, p31bar2, p312bar, and r
get.psi.RR <- function(p312, p31bar2, p312bar, ratio = 1){
  psi = ratio*(p31bar2*p312bar)/p312
  return(c(p312 = p312, p312bar = p312bar, 
           p31bar2 = p31bar2, p31bar2bar = psi))
}


## A function to compute dirichlet credible interval 
## under the RR-type constraint psi = r*p31bar2*p312bar/p312  ##
Dir_CI_RR <- function(df_aggre, n.post, 
                      rr = 1,
                      a = 0, b = 0){
  
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  
  f <- function(nc, parm, rr, a, b){
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- c("n111", "n110", "n101", 
                             "n100", "n011", "n010", "n001")
    p312.hat = (a + dat_simed$n111)/(a + b + dat_simed$n111 + dat_simed$n110)
    p312bar.hat = (a + dat_simed$n101)/(a + b + dat_simed$n101 + dat_simed$n100)
    p31bar2.hat = (a + dat_simed$n011)/(a + b + dat_simed$n011 + dat_simed$n010)
    
    estpsi <- rr*(p312bar.hat*p31bar2.hat/p312.hat)
    Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
    pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)), rr = rr,
               a = a, b = b)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}

## A function to conduct the proposed uncertainty analysis under 
## the RR-type constraint while allowing r to follow parametric assumptions
## this function only allows for r to follow normal or uniform distributions 
## with pre-specified parameters 
Dir_CI_RR_unce <- function(df_aggre, n.post, 
                           a = 0, b = 0,
                           assume.dist, a.dist, b.dist){
  # INPUT:
  # df_aggre: data framework contains observed 3-stream CRC data
  # n.post: the number of posterior samples drawn for each given value of r
  # a, b: parameters in Beta priors; 
  # e.g., a = b = 0 (unadjusted estimator)
  # e.g., a = 1, b = 0 (bias-corrected estimator based on Beta(1,0) priors)
  # assum.dist = "Uniform" or "Normal"
  # a.dist: mean in the normal dist or lower limit in the uniform dist
  # b.dist: sd in the normal dist or upper limit in the uniform dist
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  
  f <- function(nc, parm, assume.dist, a.dist, b.dist, a, b){
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- c("n111", "n110", "n101", 
                             "n100", "n011", "n010", "n001")
    p312.hat = (dat_simed$n111 + a)/(a + b + dat_simed$n111 + dat_simed$n110)
    p312bar.hat = (dat_simed$n101 + a)/(a + b + dat_simed$n101 + dat_simed$n100)
    p31bar2.hat = (dat_simed$n011 + a)/(a + b + dat_simed$n011 + dat_simed$n010)
    if(assume.dist == "Uniform"){
      rr <- runif(100, a.dist, b.dist)
    }else if(assume.dist == "Normal"){
      rr <- rnorm(100, a.dist, b.dist)
    }
    estpsi <- rr*(p312bar.hat*p31bar2.hat)/p312.hat
    Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
    pcj <- sapply(Nhat, function(y) sum(as.numeric(dat_simed[1, ]/y)))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               assume.dist = assume.dist, a = a, b = b,
               a.dist = a.dist, b.dist = b.dist)
  pcj <- as.numeric(pcj)
  
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}


# A function to fit all possible log-linear models 
fit_logl <- function(dat, nstreams, predictors){
  # INPUT: 
  # dat: counts and capture history indicators in data frame 
  # nstreams: the number of data streams
  # predictors: names of all possible predictors 
  
  nc = sum(dat$count)
  
  ## generate formula used in glm function #
  gene_input <- function(x){
    paste0("count ~ ", paste0(x, collapse = "+"))
  }
  models_form_all <- c("count ~ 1",
                       unlist(lapply(1:(2^nstreams-2), 
                                     function(n) 
                                       combn(predictors, n, FUN = gene_input))))
  ## generate names of used predictors  
  f <- function(x) paste0(x, collapse = ", ")
  predictors_all <- unlist(lapply(1:(2^nstreams-2), 
                                  function(n) combn(predictors, n, FUN = f)))
  predictors_all <- c("Intercept only", predictors_all)
  
  ## fit all possible models #
  fits_all <- lapply(models_form_all, 
                     function(y) glm(y, family = "poisson", data = dat))
  
  ## A function to obtain AIC, log-lik, estimated unobserved cell counts and N 
  get_res <- function(y){
    re = summary(y)
    re_vec = c(exp(re$coefficients[1, 1]), 
               exp(re$coefficients[1, 1] + 
                     c(-1.96, 1.96)*re$coefficients[1, 2]))
    re_vec = c(re_vec, nc + re_vec, AIC(y), as.numeric(logLik(y)))
    re_vec = re_vec[c(4:8, 1:3)]
    return(re_vec)
  }
  re_mat <- as.data.frame(do.call(rbind, lapply(fits_all, get_res)))
  colnames(re_mat) <- c("N", "N_lci", "N_uci", "AIC", "loglik", 
                        "n_unobs", "n_unobs_lci", "n_unobs_uci")
  re_mat <- data.frame(Model = predictors_all, re_mat)
  
  return(list(fits_all = fits_all, re_mat = re_mat))
  
}



## A function to compute Dirichlet-based credible interval 
## under the assumption OR31|2 = OR31|2bar ##
## i.e., [p312/(1-p312)]/[p312bar/(1-p312bar)] = 
## [p31bar2/(1-p31bar2)]/[psi/(1-psi)]
Dir_CI_OR <- function(df_aggre, n.post, 
                      or = 1, a = 0, b = 0){
  
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  
  f <- function(nc, parm, or, a, b){
    # dat_simed <- as.data.frame(t(ceiling(nc*parm)))
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- c("n111", "n110", "n101", 
                             "n100", "n011", "n010", "n001")
    p312.hat = (dat_simed$n111 + a)/(a + b + dat_simed$n111 + dat_simed$n110)
    p312bar.hat = (dat_simed$n101 + a)/(a + b + dat_simed$n101 + dat_simed$n100)
    p31bar2.hat = (dat_simed$n011 + a)/(a + b + dat_simed$n011 + dat_simed$n010)
    or.inter <- or*(p312bar.hat/(1-p312bar.hat))*
      (p31bar2.hat/(1-p31bar2.hat))/(p312.hat/(1-p312.hat))
    estpsi <- or.inter/(1 + or.inter)
    Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
    pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               or = or, a = a, b = b)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}




## A function to compute Dirichlet credible interval 
## under the assumption psi = p312, p31bar2, p312bar 
## while allowing the ratio ~ Unif/Normal ##
Dir_CI_phi_S3_unce <- function(df_aggre, n.post, 
                               a = 0, b = 0,
                               parm.name,
                               assume.dist, a.dist, b.dist){
  
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  
  f <- function(nc, parm, assume.dist, 
                a.dist, b.dist, a, b,
                parm.name){
    # dat_simed <- as.data.frame(t(ceiling(nc*parm)))
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- c("n111", "n110", "n101", 
                             "n100", "n011", "n010", "n001")
    if(parm.name == "p312"){
      psi.hat = (dat_simed$n111 + a)/(a + b + dat_simed$n111 + dat_simed$n110)
    }else if(parm.name == "p31bar2"){
      psi.hat = (dat_simed$n011 + a)/(a + b + dat_simed$n011 + dat_simed$n010)
    }else if(parm.name == "p312bar"){
      psi.hat = (dat_simed$n101 + a)/(a + b + dat_simed$n101 + dat_simed$n100)
    }
    if(assume.dist == "Uniform"){
      rr <- runif(100, a.dist, b.dist)
    }else if(assume.dist == "Normal"){
      rr <- rnorm(100, a.dist, b.dist)
    }
    estpsi <- psi.hat/rr
    Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
    pcj <- sapply(Nhat, function(y) sum(as.numeric(dat_simed[1, ]/y)))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               assume.dist = assume.dist, a = a, b = b,
               a.dist = a.dist, b.dist = b.dist,
               parm.name = parm.name)
  pcj <- as.numeric(pcj)
  
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}



Dir_CI_OR_unce <- function(df_aggre, n.post, 
                           a = 0, b = 0,
                           assume.dist, a.dist, b.dist){
  
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  
  f <- function(nc, parm, assume.dist, a.dist, b.dist, a, b){
    # dat_simed <- as.data.frame(t(ceiling(nc*parm)))
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- c("n111", "n110", "n101", 
                             "n100", "n011", "n010", "n001")
    p312.hat = (dat_simed$n111 + a)/(a + b + dat_simed$n111 + dat_simed$n110)
    p312bar.hat = (dat_simed$n101 + a)/(a + b + dat_simed$n101 + dat_simed$n100)
    p31bar2.hat = (dat_simed$n011 + a)/(a + b + dat_simed$n011 + dat_simed$n010)
    if(assume.dist == "Uniform"){
      or <- runif(100, a.dist, b.dist)
    }else if(assume.dist == "Normal"){
      or <- rnorm(100, a.dist, b.dist)
    }
    or.inter <- or*(p312bar.hat/(1-p312bar.hat))*
      (p31bar2.hat/(1-p31bar2.hat))/(p312.hat/(1-p312.hat))
    estpsi <- or.inter/(1 + or.inter)
    Nhat <- ceiling(sum(dat_simed[1,1:6]) + c(dat_simed[1,7])/estpsi)
    pcj <- sapply(Nhat, function(y) sum(as.numeric(dat_simed[1, ]/y)))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               assume.dist = assume.dist, a = a, b = b,
               a.dist = a.dist, b.dist = b.dist)
  pcj <- as.numeric(pcj)
  
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}



##### Functions for 4-stream CRC data
## A function to compute dirichlet credible interval 
## under the assumptions e.g., psi = p4123
## Beta(1,0) bias-correction 
Dir_CI_phi_S4_beta <- function(df_aggre, n.post, 
                               a = 0, b = 0,
                               parm.name = "p4123",
                               ratio.val){
  
  hist.name = as.data.frame( expand.grid(c(1,0), c(1,0), c(1,0), c(1,0)) )
  hist.name = hist.name[-16, ]
  hist.name <- hist.name[order(hist.name$Var1, hist.name$Var2,
                               hist.name$Var3, hist.name$Var4, decreasing = T),]
  hist.name <- apply(hist.name, 1, 
                     function(x) paste0("n", paste(x, collapse = "")))
  
  # a, b are two parameters for beta prior, 
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  f <- function(nc, parm, a, b, parm.name, ratio.val){
    # dat_simed <- as.data.frame(t(round(nc*parm)))
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- hist.name
    if(parm.name == "p4123"){
      nu.val = dat_simed[["n1111"]]
      de.val = dat_simed[["n1110"]]
    }else if(parm.name == "p41bar23"){
      nu.val = dat_simed[["n0111"]]
      de.val = dat_simed[["n0110"]]
    }else if(parm.name == "p412bar3"){
      nu.val = dat_simed[["n1011"]]
      de.val = dat_simed[["n1010"]]
    }else if(parm.name == "p4123bar"){
      nu.val = dat_simed[["n1101"]]
      de.val = dat_simed[["n1100"]]
    }else if(parm.name == "p41bar2bar3"){
      nu.val = dat_simed[["n0011"]]
      de.val = dat_simed[["n0010"]]
    }else if(parm.name == "p41bar23bar"){
      nu.val = dat_simed[["n0101"]]
      de.val = dat_simed[["n0100"]]
    }else if(parm.name == "p412bar3bar"){
      nu.val = dat_simed[["n1001"]]
      de.val = dat_simed[["n1000"]]
    }
    psi.hat <- (nu.val + a)/(a + b + nu.val + de.val)
    estpsi <- psi.hat/ratio.val
    Nhat <- ceiling(sum(dat_simed[1,1:14]) + c(dat_simed[1,15])/estpsi)
    pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               a = a, b = b,
               ratio.val = ratio.val,
               parm.name = parm.name)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}


Dir_CI_phi_S4_unce <- function(df_aggre, n.post, 
                               a = 0, b = 0,
                               parm.name,
                               assume.dist, a.dist, b.dist){
  
  hist.name = as.data.frame( expand.grid(c(1,0), c(1,0), c(1,0), c(1,0)) )
  hist.name = hist.name[-16, ]
  hist.name <- hist.name[order(hist.name$Var1, hist.name$Var2,
                               hist.name$Var3, hist.name$Var4, decreasing = T),]
  hist.name <- apply(hist.name, 1, 
                     function(x) paste0("n", paste(x, collapse = "")))
  
  # a, b are two parameters for beta prior, 
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  f <- function(nc, parm, a, b, 
                assume.dist, a.dist, b.dist,
                parm.name){
    # dat_simed <- as.data.frame(t(round(nc*parm)))
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- hist.name
    if(parm.name == "p4123"){
      nu.val = dat_simed[["n1111"]]
      de.val = dat_simed[["n1110"]]
    }else if(parm.name == "p41bar23"){
      nu.val = dat_simed[["n0111"]]
      de.val = dat_simed[["n0110"]]
    }else if(parm.name == "p412bar3"){
      nu.val = dat_simed[["n1011"]]
      de.val = dat_simed[["n1010"]]
    }else if(parm.name == "p4123bar"){
      nu.val = dat_simed[["n1101"]]
      de.val = dat_simed[["n1100"]]
    }else if(parm.name == "p41bar2bar3"){
      nu.val = dat_simed[["n0011"]]
      de.val = dat_simed[["n0010"]]
    }else if(parm.name == "p41bar23bar"){
      nu.val = dat_simed[["n0101"]]
      de.val = dat_simed[["n0100"]]
    }else if(parm.name == "p412bar3bar"){
      nu.val = dat_simed[["n1001"]]
      de.val = dat_simed[["n1000"]]
    }
    psi.hat <- (nu.val + a)/(a + b + nu.val + de.val)
    if(assume.dist == "Uniform"){
      rr <- runif(100, a.dist, b.dist)
    }else if(assume.dist == "Normal"){
      rr <- rnorm(100, a.dist, b.dist)
    }
    estpsi <- psi.hat/rr
    Nhat <- ceiling(sum(dat_simed[1,1:14]) + c(dat_simed[1,15])/estpsi)
    pcj <- sapply(Nhat, function(y) sum(as.numeric(dat_simed[1, ]/y)))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               a = a, b = b,
               assume.dist = assume.dist,
               a.dist = a.dist, b.dist = b.dist,
               parm.name = parm.name)
  pcj <- as.numeric(pcj)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}




## A function to compute dirichlet credible interval 
## under the RR-type assumptions for 4-catch stream
## Beta(1,0) bias-correction 
Dir_CI_RR_S4_beta <- function(df_aggre, n.post, 
                              a = 0, b = 0,
                              parm.name = "RR_1bar",
                              ratio.val){
  # parm.name:
  # RR_1bar: p41bar23/p41bar23bar = rr*p41bar2bar3/p41bar2bar3bar
  # RR_2bar: p412bar3/p412bar3bar = rr*p41bar2bar3/p41bar2bar3bar
  
  hist.name = as.data.frame( expand.grid(c(1,0), c(1,0), c(1,0), c(1,0)) )
  hist.name = hist.name[-16, ]
  hist.name <- hist.name[order(hist.name$Var1, hist.name$Var2,
                               hist.name$Var3, hist.name$Var4, decreasing = T),]
  hist.name <- apply(hist.name, 1, 
                     function(x) paste0("n", paste(x, collapse = "")))
  
  # a, b are two parameters for beta prior, 
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  f <- function(nc, parm, a, b, parm.name, ratio.val){
    # dat_simed <- as.data.frame(t(round(nc*parm)))
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- hist.name
    if(parm.name == "RR_1bar"){
      p41bar23.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                   a = 1, b = 0, parm.name = "p41bar23")
      p41bar23bar.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p41bar23bar")
      p41bar2bar3.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p41bar2bar3")
      psi.hat = ratio.val*p41bar23bar.i*p41bar2bar3.i/p41bar23.i
    }else if(parm.name == "RR_2bar"){
      p412bar3.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                   a = 1, b = 0, parm.name = "p412bar3")
      p412bar3bar.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p412bar3bar")
      p41bar2bar3.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p41bar2bar3")
      psi.hat = ratio.val*p412bar3bar.i*p41bar2bar3.i/p412bar3.i
    }
    estpsi <- psi.hat
    Nhat <- ceiling(sum(dat_simed[1,1:14]) + c(dat_simed[1,15])/estpsi)
    pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               a = a, b = b,
               ratio.val = ratio.val,
               parm.name = parm.name)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}


Dir_CI_RR_S4_unce <- function(df_aggre, n.post, 
                              a = 0, b = 0,
                              parm.name,
                              assume.dist, a.dist, b.dist){
  
  hist.name = as.data.frame( expand.grid(c(1,0), c(1,0), c(1,0), c(1,0)) )
  hist.name = hist.name[-16, ]
  hist.name <- hist.name[order(hist.name$Var1, hist.name$Var2,
                               hist.name$Var3, hist.name$Var4, decreasing = T),]
  hist.name <- apply(hist.name, 1, 
                     function(x) paste0("n", paste(x, collapse = "")))
  
  # a, b are two parameters for beta prior, 
  pstarcond <- rdirichlet(n.post, 
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001 
  f <- function(nc, parm, a, b, 
                assume.dist, a.dist, b.dist,
                parm.name){
    # dat_simed <- as.data.frame(t(round(nc*parm)))
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- hist.name
    if(parm.name == "RR_1bar"){
      p41bar23.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                   a = 1, b = 0, parm.name = "p41bar23")
      p41bar23bar.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p41bar23bar")
      p41bar2bar3.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p41bar2bar3")
      psi.hat = p41bar23bar.i*p41bar2bar3.i/p41bar23.i
    }else if(parm.name == "RR_2bar"){
      p412bar3.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                   a = 1, b = 0, parm.name = "p412bar3")
      p412bar3bar.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p412bar3bar")
      p41bar2bar3.i = get.psi.beta.S4(df_aggre = dat_simed, 
                                      a = 1, b = 0, parm.name = "p41bar2bar3")
      psi.hat = p412bar3bar.i*p41bar2bar3.i/p412bar3.i
    }
    if(assume.dist == "Uniform"){
      rr <- runif(100, a.dist, b.dist)
    }else if(assume.dist == "Normal"){
      rr <- rnorm(100, a.dist, b.dist)
    }
    estpsi <- psi.hat*rr
    Nhat <- ceiling(sum(dat_simed[1,1:14]) + c(dat_simed[1,15])/estpsi)
    pcj <- sapply(Nhat, function(y) sum(as.numeric(dat_simed[1, ]/y)))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  pcj <- apply(pstarcond, 1, f, 
               nc = as.numeric(sum(df_aggre)),
               a = a, b = b,
               assume.dist = assume.dist,
               a.dist = a.dist, b.dist = b.dist,
               parm.name = parm.name)
  pcj <- as.numeric(pcj)
  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = as.numeric(quantile(Npost, c(0.025, 0.975), na.rm = T))))
}


## A function to fit the alternative framework incoporating referral ##
fit_referral <- function(dat, se = TRUE){
  ll_s3 <- function(parms){
    # q: proportion of individuals referred from S1 to S3
    # referral of a proportion q individuals from S1 to S3
    p1 = parms[['p1']]
    p21 = parms[['p21']]
    p21bar = parms[['p21bar']]
    p312 = parms[['p312']]
    p312bar = parms[['p312bar']]
    q = parms[['q']]
    
    psi = p312bar
    p31bar2 = p312
    
    # n111 contains two part, one is the cases captured by all three streams,
    # another is the cases captured by S1, S2 and not captured by S3 
    # but is referred from S1 to S3
    # p111 = p1*p21*p312 + p1*p21*(1-p312)*q
    lamb111 = max(p1*p21*(p312 + q*(1-p312)), 1e-10)
    lamb110 = max(p1*p21*(1-p312)*(1-q), 1e-10)
    lamb101 = max(p1*(1-p21)*(p312bar + q*(1-p312bar)), 1e-10)
    lamb100 = max(p1*(1-p21)*(1-p312bar)*(1-q), 1e-10)
    lamb011 = max((1-p1)*p21bar*p31bar2, 1e-10)
    lamb010 = max((1-p1)*p21bar*(1-p31bar2), 1e-10)
    lamb001 = max((1-p1)*(1-p21bar)*psi, 1e-10)
    
    n111 = df_aggre[1, 1]; n110 = df_aggre[1, 2]
    n101 = df_aggre[1, 3]; n100 = df_aggre[1, 4]
    n011 = df_aggre[1, 5]; n010 = df_aggre[1, 6]
    n001 = df_aggre[1, 7]
    
    ll = n111*log(lamb111) + n110*log(lamb110) + n101*log(lamb101)+
      n100*log(lamb100) + n011*log(lamb011)+
      n010*log(lamb010) + n001*log(lamb001) - 
      (n111 + n110 + n101 + n100 + 
         n011 + n010 + n001)*log(1-(1-p1)*(1-p21bar)*(1-psi))
    
    return(-1*ll)
    
  }
  parnames(ll_s3) <- c("p1", "p21", "p21bar", "p312", "p312bar", "q")
  
  initial_parm = c(p1 = 0.5, p21 = 0.5, p21bar = 0.5, 
                   p312 = 0.5, p312bar = 0.5, q = 0.5) 
  if(se == TRUE){
    fit <- mle2(ll_s3, start = initial_parm, data = list(df_aggre = dat), 
                method="L-BFGS-B", 
                lower=c(p1 = 1e-6, p21 = 1e-6, p21bar = 1e-6,
                        p312 = 1e-6, p312bar = 1e-6, q = 1e-6),
                upper = c(p1 = 1-1e-6, p21 = 1-1e-6, p21bar = 1-1e-6, 
                          p312 = 1-1e-6, p312bar = 1-1e-6, q = 1-1e-6))
  }else{
    fit <- mle2(ll_s3, start = initial_parm, data = list(df_aggre = dat), 
                method="L-BFGS-B", 
                lower=c(p1 = 1e-6, p21 = 1e-6, p21bar = 1e-6,
                        p312 = 1e-6, p312bar = 1e-6, q = 1e-6),
                upper = c(p1 = 1-1e-6, p21 = 1-1e-6, p21bar = 1-1e-6, 
                          p312 = 1-1e-6, p312bar = 1-1e-6, q = 1-1e-6),
                skip.hessian = T)
  }
  
  return(fit)
}


## A function to compute dirichlet credible interval 
## under the assumption  with referral ##
Dir_CI_referral <- function(df_aggre, n.post){
  
  f <- function(nc, parm){
    dat_simed <- as.data.frame(t(nc*parm))
    colnames(dat_simed) <- c("n111", "n110", "n101", 
                             "n100", "n011", "n010", "n001")
    
    fit.inter = fit_referral(dat = dat_simed, se = FALSE)
    Nhat = get_Nhatpsi_S3(dat = dat_simed, psi = coef(fit.inter)["p312bar"])
    pcj <- sum(as.numeric(dat_simed[1, ]/Nhat))
    pcj[pcj > 1] <- 1
    return(pcj)
  }
  
  no_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', unset="1"))
  registerDoParallel(cores = no_cores)
  cl <- makeCluster(no_cores, type="FORK")
  pstarcond <- rdirichlet(n.post,
                          alpha = as.numeric(df_aggre)+0.5) # p111, ..., p001
  pcj <- parApply(cl = cl, pstarcond, 1, f,
                  nc = as.numeric(sum(df_aggre)))
  stopCluster(cl)

  Nnew <- ceiling(rowSums(df_aggre)/pcj)
  ncnew <- apply(cbind(Nnew, pcj), 1, function(x) 
    rbinom(1, x[1], x[2])) # not condition on nc 
  Npost <- ncnew/pcj 
  return(list(post = Npost, 
              median = median(Npost, na.rm = T),
              mean = mean(Npost, na.rm = T),
              CI = quantile(Npost, c(0.025, 0.975), na.rm = T)))
}
