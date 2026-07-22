
library(parallel)

cl <- makeCluster (14)
clusterEvalQ (cl,{
  
  setwd("C:\\Users\\royr2\\OneDrive\\Simulation")
  
  source("EstNull.func.R.txt")
  
  source("epsest.func.R.txt")
  
  source("lin.itp.R.txt")
  
  source("adpt.cutz.R.txt")
  
  source("adaptZ.func.R.txt")
  
  source("Seq.adap.txt")
  
  source("lfdr.gen.R.txt")
  
  source("jin.cai.pi0.R.txt")
  
  source("2.sample.f.or.R.txt")
  
})



two.side.t(x[1], x[2] , 10 , 0.5 , 0 , 1 , -1 , 0.05, 0.1)



matrix <- cbind ( rep ( c ( 100 , 500 , 1000 , 5000 ) , each = 3 ) ,rep(c(0.2,0.5,0.8),4) )



f.test.2 <- function ( X ) {
  
  clusterExport(cl,'X')
  
  time.s <- Sys.time ( )
  
  r1 <- parSapply ( cl , 1 : 200, function ( i , ... ) { x <- f.test.or(X[1], X[2] , 20 , 0.5, 0 , 0 , 1 , 1.5 , 0.05, 0.1) } )
  
  time.f <- Sys.time ( )
  
  
  tp <- time.f - time.s
  
  c ( rowMeans ( r1 ) , apply ( r1 , 1 , sd ) / sqrt ( 200 ), tp )
  
}

r1 = apply ( matrix , 1 , f.test.2 )

r1

X <- matrix[11,]

t.test.2(x)

clusterExport(cl,'X')
