###############################################################
# Example: Two-sided One-sample Student's t-test
#
# This script performs the simulation study for the
# sequential two-sided one-sample Student's t-test using the
# Data-driven Intersection (DI) procedure.
#
# Required code:
#   one.sample.t.di.R
#
# Output:
#   A vector containing:
#     • Monte Carlo estimates of sample size, FDR and FNR
#     • Monte Carlo standard errors of these estimates
#
###############################################################


library(parallel)

ncore = detectCores ( ) - 2

cl = makeCluster ( ncore )

clusterSetRNGStream ( cl , 12345 )

clusterEvalQ ( cl , {
  
  source ( "one.sample.t.di.R" )
  
})


one.sample.t.di ( 100 , c ( 0.2 , 0.6 , 0.2 ) ,
                 5 , 0 , 1 , -1.5 , 0.05 )



matrix <- cbind (
  
  rep ( c ( 100 , 500 , 1000 , 5000 ) , each = 3 ) ,
  
  rbind (
    
    c ( 0.2 , 0.6 , 0.2 ) ,
    c ( 0.5 , 0.1 , 0.4 ) ,
    c ( 0.8 , 0.15 , 0.05 ) ,
    
    c ( 0.2 , 0.6 , 0.2 ) ,
    c ( 0.5 , 0.1 , 0.4 ) ,
    c ( 0.8 , 0.15 , 0.05 ) ,
    
    c ( 0.2 , 0.6 , 0.2 ) ,
    c ( 0.5 , 0.1 , 0.4 ) ,
    c ( 0.8 , 0.15 , 0.05 ) ,
    
    c ( 0.2 , 0.6 , 0.2 ) ,
    c ( 0.5 , 0.1 , 0.4 ) ,
    c ( 0.8 , 0.15 , 0.05 )
    
  )
  
)



t.test.1 <- function ( X ) {
  
  m <- X[1]
  
  p <- X[-1]
  
  clusterExport ( cl , c ( "m" , "p" ) )
  
  time.s <- Sys.time ( )
  
  r1 <- parSapply (
    
    cl ,
    1 : 200 ,
    
    function ( i , ... ) {
      
      x <- one.sample.t.di (
        
        m , p , 5 ,
        0 , 1 , -1.5 ,
        0.05
        
      )
      
    }
    
  )
  
  time.f <- Sys.time ( )
  
  tp <- time.f - time.s
  
  c (
    
    rowMeans ( r1 ) ,
    
    apply ( r1 , 1 , sd ) / sqrt ( 200 ) ,
    
    tp
    
  )
  
}



r1 = apply ( matrix , 1 , t.test.1 )

r1
