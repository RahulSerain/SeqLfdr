
###############################################################
# Example 3: Normal Variance Testing
#
# Corresponds to:
#   Example 3 (Ex3), Table 1
#
# Problem:
#   Sequential large-scale multiple testing for
#   H0: σ = σ0
#   H1: σ = σ1
#
# This file contains the data-generation functions and the
# main simulation function for the Normal Variance testing
# problem.
#
# Main function:
#   norm.sigma.hypo(m, p, n, mu0, mu1, sigma0, sigma1,
#                   alpha, beta)
#
# Parameters:
#   m     - number of hypotheses
#   p     - proportion of hypotheses under the alternative
#   n     - initial number of observations per hypothesis
#   mu0   - mean under the null hypothesis
#   mu1   - mean under the alternative hypothesis
#   sigma0 - standard deviation under the null hypothesis
#   sigma1 - standard deviation under the alternative hypothesis
#   alpha - target false discovery rate
#   beta  - target false non-discovery rate
#
# The simulation compares:
#   • the data-driven local-FDR procedure, and
#   • the oracle local-FDR procedure.
#
# Functions:
#   dat.gen()           - Generates one observation for each
#                        hypothesis.
#
#   data.generate.n()  - Generates n observations for each
#                        hypothesis.
#
#   data.generate.x()  - Adds one new observation to the
#                        existing observations.
#
#   test.stat()         - Calculates the test statistic based
#                        on the sample variance.
#
#   or.lfr()            - Calculates the oracle local false
#                        discovery rate.
#
#   z_score()           - Transforms the test statistic to a
#                        normal score.
#
#   norm.sigma.hypo()   - Performs the sequential simulation
#                        and calculates the operating
#                        characteristics of the procedures.
#
# Output of norm.sigma.hypo():
#   • sample size, FDR and FNR for the data-driven procedure,
#   • sample size, FDR and FNR for the oracle procedure.
#
# Required functions:
#   lfdr.gen()
#   LfdrI()
#
# These functions are sourced from the common Functions/
# directory by the simulation driver.
###############################################################

dat.gen = function ( theta , mu0 , mu1 , sigma0 , sigma1) {
  
  m <- length ( theta )
  
  X <- (1-theta) * rnorm ( m , mu0 , sigma0 ) + theta * rnorm ( m , mu1 , sigma1 ) 
  
  return ( X )
  
}

data.generate.n <- function ( theta , mu0 , mu1 , sigma0 , sigma1, n ) {
  
  replicate ( n , dat.gen (  theta , mu0 , mu1 , sigma0 , sigma1 ) ) 
  
}

data.generate.x <- function (  theta , mu0 , mu1 , sigma0 , sigma1, X ) {
  
  u <-  dat.gen (  theta , mu0 , mu1 , sigma0 , sigma1 ) 
  
  cbind ( X , u )
  
}

test.stat = function ( x , mu0 , mu1, sigma0 )   {
  #qnorm(pt(sqrt(length(x))*mean(x)/sd(x),df=(length(x)-1)))
  
  n <- length(x)
  
  mu <- c(mu0 , mu1)
  
  biyog <- abs(mean(x) - mu)
  
  mu <- mu[which(biyog==min(biyog))]
  
  zi <- x-mu
  
  v <- sum (zi^2)/sigma0^2
  
  #p <- pcauchy ( z1 , n*mu0)
  
  #z2 <- qnorm ( p )
  
  #z2 <- qnorm(p)
  #p <- pnorm(z)
  
  #ret <- c ( z2 )
  
  #y <- cbind(z1,z2)
  
  return ( v)
  
}

or.lfr <- function (z , sigma0 , sigma1 , p , n ) {
  
  z.n <- ( 1 - p ) * dchisq ( z , df = n )/sigma0^2
  
  z.a <- p * dchisq( z * sigma0^2 / sigma1^2 ,  df = n )/sigma1^2
  
  lfr <- z.n / ( z.n + z.a )
  
}

z_score <- function( V, n ){
  
  p.val <- pchisq ( V , df = n , lower.tail = FALSE )
  
  z.sc <- qnorm ( p.val )
  
  return (z.sc)
  
}

norm.sigma.hypo = function ( m , p , n , mu0 , mu1 , sigma0, sigma1 , alpha , beta ) {
  
  theta <- rbinom ( m , 1 , p )
  
  X <- data.generate.n ( theta , mu0 , mu1 , sigma0 , sigma1, n )
  
  t.s <- function ( x ) {
    
    test.stat  ( x , mu0 , mu1, sigma0 )
    
  }
  
  V <- apply ( X , 1 , t.s )
  
  Z <- z_score ( V , n )
  
  lfr1 <- or.lfr ( V , sigma0 , sigma1 , p , n )
  
  lfr2 <- lfdr.gen ( Z )
  
  c1 <- LfdrI ( lfr1 , alpha , beta )
  
  c2 <- LfdrI ( lfr2 , alpha , beta )
  
  l1 <- c1 $ cutoff [ 1 ] - c1 $ cutoff [ 2 ] 
  
  l2 <- c2 $ cutoff [ 1 ] - c2 $ cutoff [ 2 ] 
  
  i = 0
  
  j = 0
  
  n1 = 0
  
  n2 = 0 
  
  #print(c(n,l1,l2))
  
  while ( l1 < 0 || l2 < 0 ) {
    
    n = n + 1
    
    X <- data.generate.x ( theta , mu0 , mu1 , sigma0 , sigma1 , X)
    
    t.s <- function ( x ) {
      
      test.stat  ( x , mu0 , mu1, sigma0 )
      
    }
    
    V <- apply ( X , 1 , t.s )
    
    Z <- z_score ( V , n )
    
    lfr1 <- or.lfr ( V , sigma0 , sigma1 , p , n )
    
    lfr2 <- lfdr.gen ( Z )
    
    c1 <- LfdrI ( lfr1 , alpha , beta )
    
    c2 <- LfdrI ( lfr2 , alpha , beta )
    
    l1 <- c1 $ cutoff [ 1 ] - c1 $ cutoff [ 2 ] 
    
    l2 <- c2 $ cutoff [ 1 ] - c2 $ cutoff [ 2 ]
    
    if ( l1 > 0 && i==0 ) {
      
      i <- 1
      
      n1 <- n
      
      a1 <- c1 $ acc.hypo
      
      r1 <- c1 $ rej.hypo
      
    }
    
    if ( l2 > 0 && j == 0) {
      
      j <- 1
      
      n2 <- n
      
      a2 <- c2 $ acc.hypo
      
      r2 <- c2 $ rej.hypo
      
    }
    #print(c(n,l1,l2))
  }
  
  if ( i == 0 ) {
    
    n1 <- n
    
    a1 <- c1 $ acc.hypo
    
    r1 <- c1 $ rej.hypo
    
  }
  
  if ( j == 0 ) {
    
    n2 <- n
    
    a2 <- c2 $ acc.hypo
    
    r2 <- c2 $ rej.hypo
    
  }
  
  D1 <- rep ( 0 , m )
  
  D1 [ r1 ] <- rep ( 1 , length ( r1 ) )
  
  fd1 <- length ( which ( D1 > theta ) ) / max ( 1 , sum ( D1 ) )
  
  fn1 <- length ( which ( D1 < theta ) ) / max ( 1 , sum ( 1 - D1 ) )
  
  D2 <- rep ( 1 , m )
  
  D2 [ a1 ] <- rep ( 0 , length ( a1 ) )
  
  fd2 <- length ( which ( D2 > theta ) ) / max ( 1 , sum ( D2 ) )
  
  fn2 <- length ( which ( D2 < theta ) ) / max ( 1 , sum ( 1 - D2 ) )
  
  D3 <- rep ( 0 , m )
  
  D3 [ r2 ] <- rep ( 1 , length ( r2 ) )
  
  fd3 <- length ( which ( D3 > theta ) ) / max ( 1 , sum ( D3 ) )
  
  fn3 <- length ( which ( D3 < theta ) ) / max ( 1 , sum ( 1 - D3 ) )
  
  D4 <- rep ( 1 , m )
  
  D4 [ a2 ] <- rep ( 0 , length ( a2 ) )
  
  fd4 <- length ( which ( D4 > theta ) ) / max ( 1 , sum ( D4 ) )
  
  fn4 <- length ( which ( D4 < theta ) ) / max ( 1 , sum ( 1 - D4 ) )
  
  c ( n1 , fd1 , fn1 , fd2 , fn2 , n2 , fd3 , fn3 , fd4 , fn4 )
  
}
