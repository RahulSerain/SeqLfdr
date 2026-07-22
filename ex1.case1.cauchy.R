library ( parallel )

ncore = detectCores ( ) - 2

cl = makeCluster ( ncore )

clusterSetRNGStream(cl, 12345)

clusterEvalQ ( cl , {
  
  source("Functions/EstNull.func.R.txt")
  
  source("Functions/epsest.func.R.txt")
  
  source("Functions/lin.itp.R.txt")
  
  source("Functions/adpt.cutz.R.txt")
  
  source("Functions/adaptZ.func.R.txt")
  
  source("Functions/Seq.adap.txt")
  
  source("Functions/lfdr.gen.R.txt")
  
  source("Functions/jin.cai.pi0.R.txt")
  
  source("cauchy.hypo.R")
  
  
} )

clusterExport(cl,"W")

mat <- cbind ( rep ( c (2500,7500,10000) , each = 3 ) , rep ( c ( 0.2 , 0.5 , 0.8 ) , 3 ) )

Nrep = 200

norm.test <- function ( X ) {
  
  time.s <- Sys.time ( )
  
  r1 <- parSapply ( cl , 1 : Nrep , function ( i , ... ) { x <- cauchy.hypo ( X[1] , X[2], 100 , 0 , 0.25 , 0.05 , 0.1 , W ) } )
  
  time.f <- Sys.time ( )
  
  tp <- time.f - time.s
  
  c ( rowMeans ( r1 ) , apply (r1,1,sd)/ sqrt(Nrep))
  
}
x <- c(2500,0.2)
t1 <- Sys.time()
t1

r <- apply ( mat , 1 , norm.test )

write.csv(t(r),"cauchy.values.csv")

stopCluster(cl)


# [,1]       [,2]       [,3]       [,4]       [,5]     [,6]         [,7]        [,8]         [,9]       [,10]
# [1,] 209.78 0.05776212 0.09992321 0.05407508 0.10286542 3.991756 0.0022866210 0.002465617 0.0022948421 0.002456727
# [2,] 291.98 0.05303310 0.08275019 0.05190071 0.08456415 2.794089 0.0009859496 0.002292527 0.0010004404 0.002379786
# [3,] 261.22 0.04963576 0.10332964 0.04746005 0.10981123 2.491830 0.0007782613 0.003021629 0.0007977076 0.003089021
# [4,] 206.64 0.05659383 0.10120967 0.05274539 0.10413759 3.241914 0.0012554066 0.001923555 0.0013505359 0.001904405
# [5,] 293.34 0.05026870 0.08459350 0.04819898 0.08772998 1.810989 0.0005954655 0.001801310 0.0006970452 0.001889482
# [6,] 260.08 0.04910415 0.10818822 0.04613145 0.11567929 2.842511 0.0007307893 0.003299320 0.0007929264 0.003051457
# [7,] 207.56 0.05733922 0.09923614 0.05298581 0.10240240 2.471579 0.0013579698 0.001771987 0.0012478437 0.001685719
# [8,] 288.34 0.04991296 0.08938443 0.04809406 0.09205187 1.864953 0.0005867635 0.001670717 0.0006399139 0.001668241
# [9,] 275.06 0.05141086 0.08383169 0.04861982 0.09018323 3.209986 0.0005735680 0.003756527 0.0007807065 0.003608693
