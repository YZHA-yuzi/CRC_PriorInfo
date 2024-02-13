README
================

In this document, we illustrate the modeling framework proposed in the
paper titled “A capture-recapture modeling framework emphasizing expert
opinion in disease surveillance” with a real three-stream HIV CRC data.
We also describe how to produce simulation results included in the paper
in this document.

# Analyze the three-stream HIV CRC data using the proposed method

read in self-defined functions

``` r
source("FUNs.R")
```

Input the three-stream HIV CRC data

``` r
df_aggre <- data.frame(n111 = 14, n110 = 66, n101 = 23, n100 = 1072,
                       n011 = 58, n010 = 729, n001 = 319)
df_aggre
```

    ##   n111 n110 n101 n100 n011 n010 n001
    ## 1   14   66   23 1072   58  729  319

``` r
# compute the number of cases identified at least once
nc.curr <- as.numeric(sum(df_aggre))
```

Fit the proposed model under the constraint
![\psi = p\_{3\|1\bar{2}}/\phi](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C1%5Cbar%7B2%7D%7D%2F%5Cphi "\psi = p_{3|1\bar{2}}/\phi"),
where
![\phi=1](https://latex.codecogs.com/png.latex?%5Cphi%3D1 "\phi=1"),
while implementing the bias-correction based on Beta(1,0) priors.

``` r
# A function to get estimated value of estimable parameters while
# implementing the bias-correction based on Beta(1,0) priors.
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
```

Compute the point estimate

- `get_Nhatpsi_S3` is a function to compute the plug-in estimator given
  in eqn. (5) in the manuscript.

``` r
res.mat <- as.data.frame(matrix(NA, ncol = 5, nrow = 1))
colnames(res.mat) <- c("Est", "lci", "uci", "lci_unce", "uci_unce")
phi.assume <- 1
pis.bc.hat <- get.psi.beta(df_aggre = df_aggre, a = 1, b = 0, 
                           parm.name = "p312bar")/phi.assume
res.mat$Est <- ceiling(get_Nhatpsi_S3(dat = df_aggre, psi = pis.bc.hat))
```

Obtain Dirichlet-based credible interval

- `Dir_CI_phi_S3_beta` is a function to obtain Dirichlet-based credible
  intervals under constraints: (1)
  ![\psi = p\_{3\|12}/\phi_1](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C12%7D%2F%5Cphi_1 "\psi = p_{3|12}/\phi_1"), (2)
  ![\psi = p\_{3\|1\bar{2}}/\phi_2](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C1%5Cbar%7B2%7D%7D%2F%5Cphi_2 "\psi = p_{3|1\bar{2}}/\phi_2"),
  and (3)
  ![\psi = p\_{3\|\bar{1}2}/\phi_3](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C%5Cbar%7B1%7D2%7D%2F%5Cphi_3 "\psi = p_{3|\bar{1}2}/\phi_3").

In the function, arguments `a` and `b` are parameters in Beta priors.
Specifically, setting `a=1` and `b=0` corresponds to the bias-correction
based on Beta(1,0) priors. The argument `phi.type` specifies which
constraint is imposed, and the argument `phi` specifies the value of the
key parameter. For example, setting `phi.type="phi2"` and `phi=1`
indicates that the constraint
![\psi = p\_{3\|1\bar{2}}/\phi_2](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C1%5Cbar%7B2%7D%7D%2F%5Cphi_2 "\psi = p_{3|1\bar{2}}/\phi_2")
with
![\phi_2 = 1](https://latex.codecogs.com/png.latex?%5Cphi_2%20%3D%201 "\phi_2 = 1")
is imposed. Finally, the argument `n.post` specifies the number of
posterior samples generated.

``` r
# set the random seed
set.seed(1234)
re.CI <- Dir_CI_phi_S3_beta(df_aggre = df_aggre, n.post = 10000,
                              phi.type = paste0("phi", 2),
                              phi = 1, a = 1, b = 0)
# make sure the lower bound is larger than the nc
# (the number of uniquely identified cases)
re.CI$CI[1] <- max(nc.curr, re.CI$CI[1])
res.mat[1, c("lci", "uci")] <- ceiling(re.CI$CI)
```

Conduct the uncertainty analysis

Recall, we are focusing on the constraint
![\psi = p\_{3\|1\bar{2}}/\phi](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C1%5Cbar%7B2%7D%7D%2F%5Cphi "\psi = p_{3|1\bar{2}}/\phi").
In this example, the key parameter is
![\phi](https://latex.codecogs.com/png.latex?%5Cphi "\phi"), and we
assume ![\phi](https://latex.codecogs.com/png.latex?%5Cphi "\phi")
follows a Uniform distribution ranging from 0.8 to 1.2 to allow for
uncertainty in this assumption.

- `Dir_CI_phi_S3_unce` is a function to conduct the uncertainty analysis
  under constraints: (1)
  ![\psi = p\_{3\|12}/\phi_1](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C12%7D%2F%5Cphi_1 "\psi = p_{3|12}/\phi_1"), (2)
  ![\psi = p\_{3\|1\bar{2}}/\phi_2](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C1%5Cbar%7B2%7D%7D%2F%5Cphi_2 "\psi = p_{3|1\bar{2}}/\phi_2"),
  and (3)
  ![\psi = p\_{3\|\bar{1}2}/\phi_3](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C%5Cbar%7B1%7D2%7D%2F%5Cphi_3 "\psi = p_{3|\bar{1}2}/\phi_3").

In the function, arguments `a` and `b` are parameters in Beta priors.
Specifically, setting `a=1` and `b=0` corresponds to the bias-correction
based on Beta(1,0) priors. The argument `parm.name` specifies which
constraint is imposed, and the argument `assume.dist` specifies the
parametric distribution introduced for the key parameter. `parm.name`
can take values of `"p312"`, `"p312bar"`, and `"p31bar2"`, which
corresponds to constraints (1)
![\psi = p\_{3\|12}/\phi_1](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C12%7D%2F%5Cphi_1 "\psi = p_{3|12}/\phi_1"),
(2)
![\psi = p\_{3\|1\bar{2}}/\phi_2](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C1%5Cbar%7B2%7D%7D%2F%5Cphi_2 "\psi = p_{3|1\bar{2}}/\phi_2"),
and (3)
![\psi = p\_{3\|\bar{1}2}/\phi_3](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C%5Cbar%7B1%7D2%7D%2F%5Cphi_3 "\psi = p_{3|\bar{1}2}/\phi_3"),
respectively. `assume.dist` can take values of `"Uniform"` and
`"Normal"`. Arguments `a.dist` and `b.dist` specify parameters in the
assumed distribution for the key parameter. For example, setting
`parm.name="p312bar"`, `assume.dist="Uniform"`, `a.dist=0.8`, and
`b.dist=1.2` implies that the uncertainty analysis is conducted under
the constraint
![\psi = p\_{3\|1\bar{2}}/\phi](https://latex.codecogs.com/png.latex?%5Cpsi%20%3D%20p_%7B3%7C1%5Cbar%7B2%7D%7D%2F%5Cphi "\psi = p_{3|1\bar{2}}/\phi")
with
![\phi \sim \text{Uniform}(0.8, 1.2)](https://latex.codecogs.com/png.latex?%5Cphi%20%5Csim%20%5Ctext%7BUniform%7D%280.8%2C%201.2%29 "\phi \sim \text{Uniform}(0.8, 1.2)").
Finally, the argument `n.post` specifies the number of posterior samples
generated for a given value of the key parameter.

``` r
set.seed(1234)
res.CI.i <- Dir_CI_phi_S3_unce(df_aggre = df_aggre,
                                 n.post = 1000,
                                 a = 1, b = 0,
                                 parm.name = "p312bar",
                                 assume.dist = "Uniform",
                                 a.dist = 0.8,
                                 b.dist = 1.2)
res.CI.i$CI[1] <- max(nc.curr, res.CI.i$CI[1])
res.mat[1, c("lci_unce", "uci_unce")] <- ceiling(res.CI.i$CI)

pander(res.mat, style = "rmarkdown", split.table = Inf)
```

|  Est  |  lci  |  uci  | lci_unce | uci_unce |
|:-----:|:-----:|:-----:|:--------:|:--------:|
| 16530 | 11766 | 23883 |  11138   |  25012   |

### Notes

The implementation of the proposed modeling framework under other
possible constraints listed in Table 1 included in the manuscript based
on this 3-stream HIV CRC data can be found in the R script
`Application_3S.R`.

# Produce simulation results and real data analysis results included in the manuscript

## Data generation

Simulated data used in simulation studies can be generated by running
the R script `Simulations_data.R` we have provided.

## R scripts

- `FUNs.R` (self-defined functions)
- `Simulations_data.R`
- `SIM_ConsA.R`
- `SIM_Inde_ConsA.R`
- `SIM_RRunce.R`
- `SIM_OR.R`
- `SIM_referral.R`
- `Application_3S.R`
- `Application_4S.R`

## Bash shell scripts

- `sim_data.sh`
- `fit_consA.sh`
- `fit_inde_consA.sh`
- `fit_RRunce.sh`
- `fit_OR.sh`
- `fit_referral.sh`

## Instructions for conducting the simulation study

There are two options to run those provided R scripts: (1) RStudio or
(2) Bash Shell scripts. When running R codes in RStudio, please create a
new R project and make sure R scripts are stored at the same folder
where the created R project `*.Rproj` is; when running R codes using
provided bash shell scripts, please set the working directory to the
folder where R scripts and bash shell scripts are located at.

Run-time was approximated using a computer with Apple M2 Pro and 32 Gb
memory.

Steps to produce simulation results presented in Section 3 are given
below.

**Step 1**. Generate data

Run `Simulations_data.R` in RStudio; or execute the bash script
`sim_data.sh` using the command `bash sim_data.sh` in terminal.

All simulated data along with detailed simulation setups (saved into csv
files) will be saved to the folder named `Data` under the working
directory.

**Step 2**. Conduct simulation studies under various simulation
scenarios

Run `SIM_ConsA.R`, `SIM_Inde_ConsA.R`, `SIM_RRunce.R`, `SIM_OR.R`, and
`SIM_referral.R` in RStudio (arguments `sce_index` and `B_index` have to
be specified); or execute bash scripts `fit_consA.sh`,
`fit_inde_consA.sh`, `fit_RRunce.sh`, `fit_OR.sh`, and
`fit_referral.sh`. The computation time varies by simulation scenarios
and the number of simulations ran in the for-loop when the parallel
computation is implemented.

Simulation results will be saved to a sub-directory named `Results`.

## Instructions for conducting the real data application

Run `Application_3S.R` and `Application_4S.R` in RStudio. Figures 1 and
2 included in the manuscript will be generated and saved to a
sub-directory named `Results`.
