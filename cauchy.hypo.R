#m <- 2000;  p <- 0.2 ; n<-20 ; mu0 = 0 ; mu1 = 0.25 ; alpha = 0.05 ; beta = 0.1
data.generate <- function ( theta , mu1 , mu0 ) {
  
  m <- length ( theta )  
  
  X <- theta * rcauchy ( m , mu1 , 1 ) + ( 1 - theta ) * rcauchy ( m , mu0 , 1 )
  
  X
  
}

data.generate.n <- function ( theta , mu1 , mu0 , n ) {
  
  replicate ( n , data.generate ( theta , mu1 , mu0 ) ) 
  
}

data.generate.x <- function ( theta , mu1 , mu0 , X ) {
  
  u <-  data.generate ( theta , mu1 , mu0 ) 
  
  cbind ( X , u )
  
}


test.stat <- function(x , mu0, mu1){
  
  t <- sum(log(dcauchy(x,mu1,1)/dcauchy(x,mu0,1)))
  
  t
  
}


z_score = function ( t , w )   {
  #qnorm(pt(sqrt(length(x))*mean(x)/sd(x),df=(length(x)-1)))
  
  emp.prb <-function (x){
    
    sum(w <= x)/ length(w)
    
  }
  
  p <- sapply(t,emp.prb)
  
  z <- qnorm(p)
  
  return (z)
  
}

or.lfr <- function (z , mu0 , mu1 , p , n ) {
  
  z.n <- ( 1 - p ) * dnorm ( z , 0 , 1 )
  
  z.a <- p * dnorm( z , sqrt(n)*(mu1-mu0) , 1 )
  
  lfr <- z.n / ( z.n + z.a )
  
}



cauchy.hypo = function ( m , p , n , mu0 , mu1 , alpha , beta , W ) {

  theta <- rbinom ( m , 1 , p )

  X <- data.generate.n ( theta , mu1 , mu0 , n)

  t.s <- function ( x ) {
  
    test.stat ( x , mu0, mu1 )
  
  }

  t <- apply ( X , 1 , t.s )
  
  w = W[,n]
  
  z <- z_score ( t , w )
  
  z[which(z==Inf)]= rep (20, length(which(z==Inf)))
  
  z[which(z==-Inf)]= rep (-20, length(which(z==-Inf)))
  
  lfr1 <- lfdr.gen (z , 0 )
 
  #lfr2 <- or.lfr ( Z , mu0 , mu1 , p , n )
 
  c1 <- LfdrI ( lfr1 , alpha , beta )
 
  #c2 <- LfdrI ( lfr2 , alpha , beta )
 
  l1 <- c1 $ cutoff [ 1 ] - c1 $ cutoff [ 2 ] 
 
  #l2 <- c2 $ cutoff [ 1 ] - c2 $ cutoff [ 2 ] 
 
  # i = 0
  # 
  # j = 0
  # 
  # n1 = 0
  # 
  # n2 = 0 
 
  while ( l1 < 0  ) {
   
     n = n + 1
   
     X <- data.generate.x ( theta , mu1 , mu0 , X)
   
     t.s <- function ( x ) {
       
       test.stat ( x , mu0, mu1 )
       
     }
     
     t <- apply ( X , 1 , t.s )
     
     w = W[,n]
     
     z <- z_score ( t , w )
     
     z[which(z==Inf)]= rep (20, length(which(z==Inf)))
     
     z[which(z==-Inf)]= rep (-20, length(which(z==-Inf)))
     
     lfr1 <- lfdr.gen (z , 0 )
     
     c1 <- LfdrI ( lfr1 , alpha , beta )
   
     l1 <- c1 $ cutoff [ 1 ] - c1 $ cutoff [ 2 ] 
   
     #l2 <- c2 $ cutoff [ 1 ] - c2 $ cutoff [ 2 ]
   
   #   if ( l1 > 0 && i==0 ) {
   #   
   #     i <- 1
   #   
   #     n1 <- n
   #   
   #     a1 <- c1 $ acc.hypo
   #   
   #     r1 <- c1 $ rej.hypo
   #   
   #   }
   # 
   #   if ( l2 > 0 && j == 0) {
   #   
   #     j <- 1
   #   
   #     n2 <- n
   #   
   #     a2 <- c2 $ acc.hypo
   #   
   #     r2 <- c2 $ rej.hypo
   #   
   #   }
   # print(c(n,l1,l2))
  }
 
  # if ( i == 0 ) {
  #  
  #   n1 <- n
  #  
     a1 <- c1 $ acc.hypo
  #  
     r1 <- c1 $ rej.hypo
  #  
  # }
  # 
  # if ( j == 0 ) {
  #  
  #   n2 <- n
  #  
  #   a2 <- c2 $ acc.hypo
  #  
  #   r2 <- c2 $ rej.hypo
  #  
  #}
 
  D1 <- rep ( 0 , m )

  D1 [ r1 ] <- rep ( 1 , length ( r1 ) )

  fd1 <- length ( which ( D1 > theta ) ) / max ( 1 , sum ( D1 ) )

  fn1 <- length ( which ( D1 < theta ) ) / max ( 1 , sum ( 1 - D1 ) )
  
  #msce1 <- sum(D1 * ( 1 - theta ) + ( 1 - D1 ) * theta)/m

  D2 <- rep ( 1 , m )

  D2 [ a1 ] <- rep ( 0 , length ( a1 ) )

  fd2 <- length ( which ( D2 > theta ) ) / max ( 1 , sum ( D2 ) )

  fn2 <- length ( which ( D2 < theta ) ) / max ( 1 , sum ( 1 - D2 ) )
  
  #msce2 <- sum(D2 * ( 1 - theta ) + ( 1 - D2 ) * theta)/m

  # D3 <- rep ( 0 , m )
  # 
  # D3 [ r2 ] <- rep ( 1 , length ( r2 ) )
  # 
  # fd3 <- length ( which ( D3 > theta ) ) / max ( 1 , sum ( D3 ) )
  # 
  # fn3 <- length ( which ( D3 < theta ) ) / max ( 1 , sum ( 1 - D3 ) )
  # 
  # #msce3 <- sum(D3 * ( 1 - theta ) + ( 1 - D3 ) * theta)/m
  # 
  # D4 <- rep ( 1 , m )
  # 
  # D4 [ a2 ] <- rep ( 0 , length ( a2 ) )
  # 
  # fd4 <- length ( which ( D4 > theta ) ) / max ( 1 , sum ( D4 ) )
  # 
  # fn4 <- length ( which ( D4 < theta ) ) / max ( 1 , sum ( 1 - D4 ) )
  # 
  #msce4 <- sum(D4 * ( 1 - theta ) + ( 1 - D4 ) * theta)/m
  
  c ( n , fd1 , fn1 , fd2 , fn2 )

  #c ( n1 , fd1 , fn1 , fd2 , fn2 , n2 , fd3 , fn3 , fd4 , fn4 )
  
  #c ( n1 , msce1 , msce2 , n2 , msce3 , msce4 )
  
}
#norm.hypo ( m , p , n , mu0 , mu1 , alpha , beta )
