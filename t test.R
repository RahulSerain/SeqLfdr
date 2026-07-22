
dat.gen = function ( theta , mu0 , mu1 , mu2 ) {
  
  theta1 <- theta[1,]
  
  theta2 <- theta[2,]
  
  theta3 <- theta[3,]
  
  m <- ncol ( theta )
  
  X <- theta1 * rnorm ( m , mu0 , 1 ) + theta2 * rnorm ( m , mu1 , 1 ) + theta3 * rnorm (m , mu2 , 1 )
  
  return ( X )

}


z_val <- function ( x , n  ) {
  
  x1 <- x [ 1 : n ]
  
  x2 <- x [ ( n + 1 ) : length ( x ) ]
  
  #t <- ( mean ( x2 ) - mean ( x1 ) ) / (sqrt ( ( ( n1 - 1 ) * ( sd ( x1 ) ) ^ 2 + ( n2 - 1 ) * ( sd ( x2 ) ) ^ 2 ) / ( n1 + n2 - 2 ) )*sqrt(1/n1+1/n2))
  
  t.t <- t.test(x1, x2, alternative="two.sided", var.equal=TRUE)
  
  z <- qnorm ( pt ( as.numeric ( t.t$statistic ) , as.numeric ( t.t$parameter ) ) )
  
  #p <- as.numeric ( t.t$p.value )
  
  return ( z )

}



data.gen.n = function ( theta , mu0 , mu1 , mu2 ,  n) {
  
  replicate ( n , dat.gen ( theta , mu0 , mu1 , mu2 ) ) 
  
  }

  
data.gen.x = function ( theta , mu0 , mu1 , mu2 , X ) {
  
  u <-  dat.gen ( theta , mu0 , mu1 , mu2 ) 
  
  cbind ( X , u )
  
}

p.func <- function ( t , n ) {
  
  p <- 2 * min ( 1 - pt (  t  , n-2) , pt ( t , ( n - 2 ) ) )
    
  return(p)
  
}

z.func <- function ( t , n ) {
  
  qnorm ( pt ( t , ( n - 2 ) ) )
  
}


two.side.t <- function (m , p , k , mu0 , mu1 , mu2 , alpha , beta ) {
  
  theta <- rmultinom ( m , 1 , p )

  n <- 2 * k
  
  X1 <- replicate ( k,rnorm (m,0,1))
    
  X2 <- data.gen.n ( theta , mu0 , mu1 , mu2 , k) 
  
  X <- cbind ( X1 , X2 )
  
  t.v <- function(x) z_val( x , n )
  
  z1 <- apply ( X , 1 , t.v )

  lfdr1 <- lfdr.gen ( z1 )

  lf <- LfdrI ( lfdr1 , alpha , beta )

  del <- lf $ cutoff [ 1 ] - lf $ cutoff [ 2 ]

  while ( del < 0 ) {
  
    k = k + 1
  
    n <- 2 * k
    
    X1 <- cbind ( X1 , rnorm ( m , 0 , 1 ) )
    
    X2 <- data.gen.x ( theta , mu0 , mu1 , mu2 , X2 ) 
  
    X <- cbind ( X1 , X2 )
    
    t.v <- function(x) z_val( x , n )
    
    z1 <- apply ( X , 1 , t.v )
  
    lfdr1 <- lfdr.gen ( z1 )
  
    lf <- LfdrI ( lfdr1 , alpha , beta )
  
    del <- lf $ cutoff [ 1 ] - lf $ cutoff [ 2 ] 
    
    #print(c ( n , del , jin.cai.pi0(z1,0,1)))
  }

  theta1 <- 1 - theta [ 1 , ]

  D1 = rep ( 0 , m )
  
  if ( is.na(length ( lf $ rej.hypo )) == 0  ) {
    
    D1 [ lf $ rej.hypo ] <-  rep ( 1 , length ( lf $ rej.hypo ) ) 
    
  }

  fr1 <- sum ( D1 > theta1) / max (sum ( D1 ) , 1 )

  fn1 <- sum (D1 < theta1) / max ( sum ( 1 - D1 ), 1 )

  D2 = rep ( 1 , m )
  
  if ( is.na(length ( lf $ acc.hypo )) == 0 ) {
    
    D2 [ lf $ acc.hypo ] <-  rep ( 0 , length ( lf $ acc.hypo ) )
    
  }

  fr2 <- sum ( D2 > theta1) / max (sum ( D2 ) , 1 )

  fn2 <- sum (D2 < theta1) / max ( sum ( 1 - D2 ), 1 )

 c( n , fr1 , fn1  , fr2 , fn2 ) 

}




