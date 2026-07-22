
dat.gen = function ( theta , mu0 , mu1 , sigma0 , sigma1) {
  
  m <- length ( theta )
  
  X <- (1-theta) * rnorm ( m , mu0 , sigma0 ) + theta * rnorm ( m , mu1 , sigma1 ) 
  
  return ( X )
  
}


f_val <- function ( x , n  ) {
  
  x1 <- x [ 1 : n ]
  
  x2 <- x [ ( n + 1 ) : length ( x ) ]
  
  n1 <- length ( x1 )
  
  n2 <- length ( x2 )
  
  #t <- ( mean ( x2 ) - mean ( x1 ) ) / (sqrt ( ( ( n1 - 1 ) * ( sd ( x1 ) ) ^ 2 + ( n2 - 1 ) * ( sd ( x2 ) ) ^ 2 ) / ( n1 + n2 - 2 ) )*sqrt(1/n1+1/n2))
  
  f <- var ( x1 ) / var ( x2 )
  
  f
  
}

f.val.x = function ( X1 , X2  ) {
  
  X <- cbind ( X1 , X2 )
  
  n <- ncol (X1)
  
  f.v <- function(x){
    
    f_val(x,n)
    
  }
  
  apply ( X , 1 , f.v )
  
}


data.gen.n = function ( theta , mu0 , mu1 , sigma0 , sigma1 ,  n) {
  
  replicate ( n , dat.gen ( theta , mu0 , mu1 , sigma0 , sigma1 ) ) 
  
}


data.gen.x = function ( theta , mu0 , mu1 , sigma0 , sigma1 , X ) {
  
  u <-  dat.gen ( theta , mu0 , mu1 , sigma0 , sigma1 ) 
  
  cbind ( X , u )
  
}

p.func <- function ( t , n1 , n2 , sigma0 , sigma1 ) {
  
   if ( sigma1 > sigma0 ) {
    
    p <- 1 - pt ( t , ( n - 2 ) )
    
  } else if ( mu1 < mu0 ) {
    
    p <- pt ( t , ( n - 2 ) )
    
  } 
  
  return(p)
  
}

z.func <- function ( f , n1 , n2 ) {
  
  qnorm ( pf ( f , n1-1 , n2-1 ) )
  
}


lfdr.f <- function (f , n1, n2, sigma0, sigma1, p) {

  g <- ( sigma1 / sigma0 )^2*f
  
  num <- (1-p)* df ( f , n1-1 , n2-1 )
  
  den <- p * df ( g , n1-1 , n2-1 )
  
  l <- num / ( num + den )
  
  l
  
}

f.test.or <- function (m , p ,n , p1 , mu0 , mu1 , sigma0 , sigma1 , alpha , beta ) {
  
  theta <- rbinom ( m , 1 , p ) 
  
  r = 0
  
  while(r <= 1  || r >= ( n - 1 ) ) {
    
    r <- rbinom ( 1 , n , p1 )
    
    n = n + 1
    
  }
  
  n = n - 1
  
  X1 <- replicate ( n - r , rnorm ( m , mu0 , sigma0 ) )
    
  X2 <- data.gen.n ( theta , mu0 , mu1 , sigma0 , sigma1 , r) 
    
  f <- f.val.x ( X1 , X2 )
  
  n1 <- ncol (X1)
  
  n2 <- ncol (X2)
  
  #z1 = z.func ( f , n1 , n2  )
  
  lfdr1 <- lfdr.f (f , n1, n2, sigma0, sigma1, p)
  
  lf <- LfdrI ( lfdr1 , alpha , beta )
  
  del <- lf $ cutoff [ 1 ] - lf $ cutoff [ 2 ]
  
  while ( del <= 0 ) {
    
    n = n+1
    
    r <- rbinom(1,1,p1)
    
    if ( r == 0 ) {
      
      X1 <- cbind ( X1 , rnorm ( m , mu0 , sigma0 ) )
      
    } else {
      
      X2 <- data.gen.x ( theta , mu0 , mu1 ,  sigma0 ,sigma1 , X2 ) 
      
    }
    
    f <- f.val.x ( X1 , X2 )
    
    n1 <- ncol (X1)
    
    n2 <- ncol (X2)
    
    #z1 <- z.func ( f , n1 , n2 )
    
    lfdr1 <- lfdr.f (f , n1, n2, sigma0, sigma1, p)
    
    lf <- LfdrI ( lfdr1 , alpha , beta )
    
    del <- lf $ cutoff [ 1 ] - lf $ cutoff [ 2 ]
    
    print(c ( n , del ))#, jin.cai.pi0(z1,0,1)))
  }
  
  #theta1 <- 1 - theta [ 1 , ]
  
  D1 = rep ( 0 , m )
  
  D1 [ lf $ rej.hypo ] <-  rep ( 1 , length ( lf $ rej.hypo ) )
  
  fr1 <- sum ( D1 > theta) / max (sum ( D1 ) , 1 )
  
  fn1 <- sum (D1 < theta) / max ( sum ( 1 - D1 ), 1 )
  
  D2 = rep ( 1 , m )
  
  D2 [ lf $ acc.hypo ] <-  rep ( 0 , length ( lf $ acc.hypo ) )
  
  fr2 <- sum ( D2 > theta) / max (sum ( D2 ) , 1 )
  
  fn2 <- sum (D2 < theta) / max ( sum ( 1 - D2 ), 1 )
  
  c( n , fr1 , fn1  , fr2 , fn2 ) 
  
}



