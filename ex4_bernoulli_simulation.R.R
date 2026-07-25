###############################################################
# Example 4: Bernoulli Testing
#
# Corresponds to:
#   Example 4 (Ex4), Table 1
#
# Problem:
#   Sequential large-scale multiple testing for
#   H0: X ~ Bernoulli(0.10)
#   H1: X ~ Bernoulli(0.15)
#
# This script performs the simulation study for Example 4
# using the Oracle Intersection (OI) and Data-driven
# Intersection (DI) procedures.
#
# Required code:
#   Bern.Hypo.alt.R
#
# Output:
#   A vector containing:
#     • Monte Carlo estimates of ASN, FDR and FNR
#       for the OI and DI procedures,
#     • Monte Carlo standard errors of these estimates.
#
###############################################################


library(parallel)

ncore <- detectCores() - 2

cl <- makeCluster(ncore)

clusterSetRNGStream(cl, 12345)

clusterEvalQ(cl, {
  
  source("Functions/EstNull.func.R.txt")
  source("Functions/epsest.func.R.txt")
  source("Functions/lin.itp.R.txt")
  source("Functions/adpt.cutz.R.txt")
  source("Functions/adaptZ.func.R.txt")
  source("Functions/Seq.adap.txt")
  source("Functions/lfdr.gen.R.txt")
  source("Functions/jin.cai.pi0.R.txt")
  
  source("Bern.Hypo.alt.R")
  
})

mat <- cbind(
  rep(c(100, 500, 1000, 5000), each = 3),
  rep(c(0.2, 0.5, 0.8), 4)
)

Nrep <- 200

bern.test <- function(x){
  
  r1 <- parSapply(
    cl,
    1:Nrep,
    function(i, ...){
      
      bern.hypo.or(
        m     = x[1],
        p     = x[2],
        m1    = 7,
        pi0   = 0.10,
        pi1   = 0.15,
        n     = 10,
        alpha = 0.05,
        beta  = 0.10
      )
      
    }
  )
  
  c(
    rowMeans(r1),
    apply(r1, 1, sd) / sqrt(Nrep)
  )
  
}

r <- apply(mat, 1, bern.test)

stopCluster(cl)
