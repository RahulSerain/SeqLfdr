###############################################################
# Example 2: Exponential Scale Testing
#
# Corresponds to:
#   Example 2 (Ex2), Table 1
#
# Problem:
#   Sequential large-scale multiple testing for
#   H0: μ = μ0
#   H1: μ = μ1
#
# This file contains the data-generation functions and the
# main simulation function for the Exponential Scale testing
# problem.
#
# Main function:
#   exp.hypo(m, p, n, mu0, mu1, alpha, beta)
#
# Parameters:
#   m     - number of hypotheses
#   p     - proportion of hypotheses under the alternative
#   n     - initial number of observations per hypothesis
#   mu0   - scale parameter under the null hypothesis
#   mu1   - scale parameter under the alternative hypothesis
#   alpha - target false discovery rate
#   beta  - target false non-discovery rate
#
# The simulation compares:
#   • the data-driven local-FDR procedure, and
#   • the oracle local-FDR procedure.
#
# Functions:
#   data.generate()    - Generates one observation for each
#                        hypothesis from the corresponding
#                        exponential distribution.
#
#   data.generate.n()  - Generates n observations for each
#                        hypothesis.
#
#   data.generate.x()  - Adds one new observation to the
#                        existing observations.
#
#   z_score()          - Calculates the test statistic and
#                        its normal transformation.
#
#   or.lfr()           - Calculates the oracle local false
#                        discovery rate.
#
#   exp.hypo()         - Performs the sequential simulation
#                        and calculates the operating
#                        characteristics of the procedures.
#
# Output of exp.hypo():
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
data.generate <- function ( theta , mu1 , mu0 ) {
  
  m <- length ( theta )  
  
  X <- theta * rexp ( m , rate = 1 / mu1  ) + ( 1 - theta ) * rexp ( m , rate = 1/ mu0  )
  
  X
  
}

data.generate.n <- function ( theta , mu1 , mu0 , n ) {
  
  replicate ( n , data.generate ( theta , mu1 , mu0 ) ) 
  
}

data.generate.x <- function ( theta , mu1 , mu0 , X ) {
  
  u <-  data.generate ( theta , mu1 , mu0 ) 
  
  cbind ( X , u )
  
}


z_score = function ( x , mu0 )   {
  #qnorm(pt(sqrt(length(x))*mean(x)/sd(x),df=(length(x)-1)))
  
  z1 <-  sum ( x ) 
  
  p <- pgamma ( z1 , shape = length ( x ) , scale = mu0 )
  
  z2 <- qnorm ( p )
  
  #z2 <- qnorm(p)
  #p <- pnorm(z)
  
  ret <- c ( z1 , z2 )
  
  return ( ret )
  
}

or.lfr <- function (z , mu0 , mu1 , p , n ) {
  
  z.n <- ( 1 - p ) * dgamma ( z , shape = n , scale = mu0 )
  
  z.a <- p * dgamma ( z , shape = n , scale = mu1 )
  
  lfr <- z.n / ( z.n + z.a )
  
}

exp.hypo = function ( m , p , n , mu0 , mu1 , alpha , beta ) {
  
  theta <- rbinom ( m , 1 , p )
  
  X <- data.generate.n ( theta , mu1 , mu0 , n)
  
  z.s <- function ( x ) {
    
    z_score ( x , mu0 )
    
  }
  
  z <- apply ( X , 1 , z.s )
  
  z1 <- z[1,]
  
  z2 <- z[2,]
  
  lfr1 <- lfdr.gen ( z2 )
  
  lfr2 <- or.lfr ( z1 , mu0 , mu1 , p , n )
  
  c1 <- LfdrI ( lfr1 , alpha , beta )
  
  c2 <- LfdrI ( lfr2 , alpha , beta )
  
  l1 <- c1 $ cutoff [ 1 ] - c1 $ cutoff [ 2 ] 
  
  l2 <- c2 $ cutoff [ 1 ] - c2 $ cutoff [ 2 ] 
  
  i = 0
  
  j = 0
  
  n1 = 0
  
  n2 = 0 
  
  #print( c(l1,l2))
  while ( l1 < 0 || l2 < 0 ) {
    
    n = n + 1
    
    X <- data.generate.x ( theta , mu1 , mu0 , X)
    
    z.s <- function ( x ) {
      
      z_score ( x , mu0 )
      
    }
    
    z <- apply ( X , 1 , z.s )
    
    z1 <- z[1,]
    
    z2 <- z[2,]
    
    lfr1 <- lfdr.gen ( z2 )
    
    lfr2 <- or.lfr ( z1 , mu0 , mu1 , p , n )
    
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
   # print( c(l1,l2))
    
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

#exp.hypo ( 2000 , 0.5 , 20 , 1 , 1.2 , 0.05 , 0.1 )


#m <- 2000 ; p <- 0.5 ; n <- 50 ; mu0 <- 1 ; mu1 <- 1.2 ; alpha <- 0.05 ; beta <- 0.1
