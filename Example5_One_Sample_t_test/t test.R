###############################################################
# Example: Two-sided One-sample Student's t-test
#
# Problem:
#   Sequential large-scale multiple testing for
#
#       H0: μ = μ0
#       H1: μ ∈ {μ1, μ2}
#
# This script implements the sequential two-sided one-sample
# Student's t-test using the Data-driven Intersection (DI)
# procedure.
#
# For each hypothesis, observations are generated from one of
# the distributions
#
#       N(μ0, 1), N(μ1, 1), or N(μ2, 1),
#
# according to the mixing proportions specified by p.
#
# At each stage, the one-sample Student's t-statistic is
# calculated and transformed to a z-score using its
# Student's t distribution. The local false discovery rates
# are estimated from the observed z-scores using lfdr.gen().
# The adaptive cutoff is then obtained using adpt.cutz().
#
# Sampling continues until the DI procedure produces a valid
# rejection/acceptance region.
#
# Procedure:
#   Data-driven Intersection (DI) procedure
#
# Required code:
#   lfdr.gen.R.txt
#   adpt.cutz.R.txt
#
# Main function:
#   one.sample.t.di()
#
# Parameters:
#   m       : Number of hypotheses
#   p       : Vector of proportions corresponding to
#             μ0, μ1 and μ2
#   k       : Initial sample size for each hypothesis
#   mu0     : Mean under the null hypothesis
#   mu1     : First alternative mean
#   mu2     : Second alternative mean
#   alpha   : Target false discovery rate level
#
# Output:
#   A vector containing:
#     • Sample size at stopping
#     • False Discovery Rate (FDR)
#     • False Non-discovery Rate (FNR)
#
###############################################################


dat.gen = function ( theta , mu0 , mu1 , mu2 ) {
  
  theta1 <- theta[1,]
  
  theta2 <- theta[2,]
  
  theta3 <- theta[3,]
  
  m <- ncol ( theta )
  
  X <- theta1 * rnorm ( m , mu0 , 1 ) +
       theta2 * rnorm ( m , mu1 , 1 ) +
       theta3 * rnorm ( m , mu2 , 1 )
  
  return ( X )
  
}


data.gen.n = function ( theta , mu0 , mu1 , mu2 , n ) {
  
  replicate ( n , dat.gen ( theta , mu0 , mu1 , mu2 ) )
  
}


data.gen.x = function ( theta , mu0 , mu1 , mu2 , X ) {
  
  u <- dat.gen ( theta , mu0 , mu1 , mu2 )
  
  cbind ( X , u )
  
}


t_val = function ( x , mu0 ) {
  
  n <- length ( x )
  
  t <- ( mean ( x ) - mu0 ) /
       ( sd ( x ) / sqrt ( n ) )
  
  return ( t )
  
}


t.val.x = function ( X , mu0 ) {
  
  t.v <- function ( x ) {
    
    t_val ( x , mu0 )
    
  }
  
  apply ( X , 1 , t.v )
  
}


z.func = function ( t , n ) {
  
  z <- qnorm ( pt ( t , n - 1 ) )
  
  return ( z )
  
}


one.sample.t.di = function ( m , p , k , mu0 , mu1 , mu2 ,
                            alpha ) {
  
  theta <- rmultinom ( m , 1 , p )
  
  X <- data.gen.n ( theta , mu0 , mu1 , mu2 , k )
  
  t <- t.val.x ( X , mu0 )
  
  z1 <- z.func ( t , k )
  
  lfdr1 <- lfdr.gen ( z1 )
  
  adpz <- adpt.cutz ( lfdr1 , alpha )
  
  del <- adpz$cutoff[1] - adpz$cutoff[2]
  
  
  while ( del < 0 ) {
    
    k <- k + 1
    
    X <- data.gen.x ( theta , mu0 , mu1 , mu2 , X )
    
    t <- t.val.x ( X , mu0 )
    
    z1 <- z.func ( t , k )
    
    lfdr1 <- lfdr.gen ( z1 )
    
    adpz <- adpt.cutz ( lfdr1 , alpha )
    
    del <- adpz$cutoff[1] - adpz$cutoff[2]
    
  }
  
  
  theta1 <- 1 - theta[1,]
  
  D <- rep ( 0 , m )
  
  if ( length ( adpz$re ) > 0 ) {
    
    D[adpz$re] <- rep ( 1 , length ( adpz$re ) )
    
  }
  
  
  fdr <- sum ( D > theta1 ) /
         max ( sum ( D ) , 1 )
  
  fnr <- sum ( D < theta1 ) /
         max ( sum ( 1 - D ) , 1 )
  
  
  c ( k , fdr , fnr )
  
}
