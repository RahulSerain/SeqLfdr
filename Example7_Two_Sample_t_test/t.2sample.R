dat.gen = function ( theta , mu0 , mu1 , mu2 ) {
  
  theta1 <- theta[1,]
  
  theta2 <- theta[2,]
  
  theta3 <- theta[3,]
  
  m <- ncol ( theta )
  
  X <- theta1 * rnorm ( m , mu0 , 1 ) + theta2 * rnorm ( m , mu1 , 1 ) + theta3 * rnorm (m , mu2 , 1 )
  
  return ( X )
  
}

z.func <- function ( t , n ) {
  
  z <- qnorm ( pt (  t  , ( n - 2 ) ) )
  
  return (z)
  
}

t_val <- function ( x , n  ) {
  
  x1 <- x [ 1 : n ]
  
  x2 <- x [ ( n + 1 ) : length ( x ) ]
  
  n1 <- length ( x1 )
  
  n2 <- length ( x2 )
  
  t <- ( mean ( x2 ) - mean ( x1 ) ) / (sqrt ( ( ( n1 - 1 ) * ( sd ( x1 ) ) ^ 2 + ( n2 - 1 ) * ( sd ( x2 ) ) ^ 2 ) / ( n1 + n2 - 2 ) )*sqrt(1/n1+1/n2))
  
}

t.val.x = function ( X1 , X2  ) {
  
  X <- cbind ( X1 , X2 )
  
  n <- ncol (X1)
  
  t.v <- function(x){
    
    t_val(x,n)
    
  }
  
  apply ( X , 1 , t.v )
  
}


data.gen.n = function ( theta , mu0 , mu1 , mu2 ,  n) {
  
  replicate ( n , dat.gen ( theta , mu0 , mu1 , mu2 ) ) 
  
}


data.gen.x = function ( theta , mu0 , mu1 , mu2 , X ) {
  
  u <-  dat.gen ( theta , mu0 , mu1 , mu2 ) 
  
  cbind ( X , u )
  
}




or.t.lfdr <- function ( t , n1 , n2  , mu0 , mu1 , mu2 , p ){
  
  n <- n1 + n2
  
  num <- 1
  
  denom <- (p [ 2 ] * dt ( t , n - 2 , ( mu0 - mu1 ) / sqrt ( 1 / n1 + 1 / n2 )   ) + p [ 3 ] * dt ( t , n - 2 , ( mu0 - mu2 ) / sqrt ( 1 / n1 + 1 / n2 ) ))/(p[1]* dt ( t , n - 2 ))
  
  num / ( num + denom )
  
}


#m <- 100; p <- c(0.2,0.6,0.2); k <- 5; mu0 <- 0; mu1 <- 1; mu2 <- -1.5; alpha <- 0.05 

two.side.t.or.sc <- function (m , p , k , mu0 , mu1 , mu2 , alpha  ) {
  
  theta <- rmultinom ( m , 1 , p ) 
  
  n <- 2 * k
  
  X1 <- replicate ( k,rnorm (m,mu0,1))
    
  X2 <- data.gen.n ( theta , mu0 , mu1 , mu2 , k) 
  
  t <- t.val.x ( X1 , X2 )
  
  n1 <- ncol ( X1 )
  
  n2 <- ncol ( X2 )
  
  z1 = z.func ( t , n )
  
  # if(length(which(z1==Inf))>0){
  #   z1[which(z1==Inf)] <- rep(20,length(which(z1==Inf)))
  # }
  # 
  # if(length(which(z1==-Inf))>0){
  #   z1[which(z1==-Inf)] <- rep(-20,length(which(z1==-Inf)))
  # }
  
  lfdr1 <- lfdr.gen ( z1 )
  
  adpz <- adpt.cutz ( lfdr1 , alpha )
  
  theta1 <- 1 - theta [ 1 , ]
  
  D <- rep( 0 , m )
  
  D [ adpz $ re ] = rep ( 1 , length ( adpz $ re ) )
  
  fdr <- sum ( D > theta1 ) / max ( sum ( D ) , 1 )
  
  fnr <- sum ( D < theta1 ) / max ( sum ( 1 - D ) , 1 )
  
  c ( n , fdr , fnr )
  
}
