###############################################################
# Example: Two-sample Student's t-test
#
# Problem:
#   Sequential large-scale multiple testing for
#   H0: μ = μ0
#   H1: μ ∈ {μ1, μ2}
#
# This script performs the simulation study for the
# two-sided Student's t-test using the adaptive
# local-FDR-based procedure.
#
# Required code:
#   Aadaptz.t.2side.R
#
# Output:
#   A vector containing:
#     • Monte Carlo estimates of ASN, FDR and FNR
#       for the two-sided procedure,
#     • Monte Carlo standard errors of these estimates.
###############################################################

library(parallel)

ncore = detectCores ( ) - 2

cl = makeCluster ( ncore )

clusterSetRNGStream(cl, 12345)

clusterEvalQ (cl,{
  
  source("Aadaptz.t.2side.R")
  
})


two.side.t.or.sc(100, c(0.2,0.6,0.2), 5, 0, 1, -1.5, 0.05)



matrix <- cbind ( rep ( c ( 100 , 500 , 1000 , 5000 ) , each = 3 ) ,
                  rbind ( c ( 0.2 , 0.6 , 0.2 ) ,
                          c ( 0.5 , 0.1 , 0.4 ) ,
                          c ( 0.8 , 0.15 , 0.05 ),
                          c ( 0.2 , 0.6 , 0.2 ) ,
                          c ( 0.5 , 0.1 , 0.4 ) ,
                          c ( 0.8 , 0.15 , 0.05 ),
                          c ( 0.2 , 0.6 , 0.2 ) ,
                          c ( 0.5 , 0.1 , 0.4 ) ,
                          c ( 0.8 , 0.15 , 0.05 ),
                          c ( 0.2 , 0.6 , 0.2 ) ,
                          c ( 0.5 , 0.1 , 0.4 ) ,
                          c ( 0.8 , 0.15 , 0.05 ) ) )



t.test.2 <- function ( X ) {
  
  a <- X[1]
  
  b <- X[-1]
  
  clusterExport(cl,c("a","b"))
  
  time.s <- Sys.time ( )
  
  r1 <- parSapply ( cl , 1 : 200 ,
                    function ( i , ... ) {
                      x <- two.side.t.or.sc(a, b, 5,
                                             0, 1, -1.5,
                                             0.05)
                    } )
  
  time.f <- Sys.time ( )
  
  tp <- time.f - time.s
  
  c ( rowMeans ( r1 ) ,
      apply ( r1 , 1 , sd ) / sqrt ( 200 ),
      tp )
  
}


r1 = apply ( matrix , 1 , t.test.2 )

r1
