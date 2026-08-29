

dat.gen <- function(theta, mu0, mu1, sigma0, sigma1) {
  
  m <- length(theta)
  
  X <- (1 - theta) * rnorm(m, mu0, sigma0) +
       theta * rnorm(m, mu1, sigma1)
  
  return(X)
}

data.gen.n <- function(theta, mu0, mu1, sigma0, sigma1, n) {
  
  replicate(
    n,
    dat.gen(theta, mu0, mu1, sigma0, sigma1)
  )
}



data.gen.x <- function(theta, mu0, mu1, sigma0, sigma1, X) {
  
  u <- dat.gen(theta, mu0, mu1, sigma0, sigma1)
  
  cbind(X, u)
}


test.stat <- function(x, n1) {
  
  x1 <- x[1:n1]
  x2 <- x[(n1 + 1):length(x)]
  
  f <- var(x1) / var(x2)
  
  return(f)
}


f.val.x <- function(X1, X2) {
  
  X <- cbind(X1, X2)
  n1 <- ncol(X1)
  
  f.s <- function(x) {
    test.stat(x, n1)
  }
  
  apply(X, 1, f.s)
}



z_score <- function(f, n1, n2) {
  
  qnorm(pf(f, n1 - 1, n2 - 1))
}



lfdr.f <- function(f, n1, n2, sigma0, sigma1, p) {
  
  g <- (sigma1 / sigma0)^2 * f
  
  num <- (1 - p) * df(f, n1 - 1, n2 - 1)
  
  den <- p * df(g, n1 - 1, n2 - 1)
  
  lfr <- num / (num + den)
  
  return(lfr)
}



f.test.or <- function(
  m, p, n, p1,
  mu0, mu1,
  sigma0, sigma1,
  alpha, beta
) {

  
  theta <- rbinom(m, 1, p)
 
  r <- 0
  
  while(r <= 1 || r >= (n - 1)) {
    
    r <- rbinom(1, n, p1)
    n <- n + 1
  }
  
  n <- n - 1
  
  X1 <- replicate(
    n - r,
    rnorm(m, mu0, sigma0)
  )
  
  X2 <- data.gen.n(
    theta,
    mu0, mu1,
    sigma0, sigma1,
    r
  )
  
  f <- f.val.x(X1, X2)
  
  n1 <- ncol(X1)
  n2 <- ncol(X2)
  lfdr1 <- lfdr.f(
    f,
    n1, n2,
    sigma0, sigma1,
    p
  )
  
  lf <- LfdrI(lfdr1, alpha, beta)
  
  del <- lf$cutoff[1] - lf$cutoff[2]
  
  while(del <= 0) {
    
    n <- n + 1
    
    r <- rbinom(1, 1, p1)
    
    if(r == 0) {
      
      X1 <- cbind(
        X1,
        rnorm(m, mu0, sigma0)
      )
      
    } else {
      
      X2 <- data.gen.x(
        theta,
        mu0, mu1,
        sigma0, sigma1,
        X2
      )
    }
    
    
    # Recalculate F statistic
    
    f <- f.val.x(X1, X2)
    
    n1 <- ncol(X1)
    n2 <- ncol(X2)
    
    
    # Recalculate LFDR
    
    lfdr1 <- lfdr.f(
      f,
      n1, n2,
      sigma0, sigma1,
      p
    )
    
    
    # Reapply LfdrI
    
    lf <- LfdrI(
      lfdr1,
      alpha, beta
    )
    
    del <- lf$cutoff[1] - lf$cutoff[2]
  }
  
  D1 <- rep(0, m)
  
  D1[lf$rej.hypo] <- 1
  
  fdr1 <- sum(D1 > theta) /
          max(sum(D1), 1)
  
  fnr1 <- sum(D1 < theta) /
          max(sum(1 - D1), 1)
  
  
  D2 <- rep(1, m)
  
  D2[lf$acc.hypo] <- 0
  
  fdr2 <- sum(D2 > theta) /
          max(sum(D2), 1)
  
  fnr2 <- sum(D2 < theta) /
          max(sum(1 - D2), 1)
  
  
  return(
    c(n, fdr1, fnr1, fdr2, fnr2)
  )
}
