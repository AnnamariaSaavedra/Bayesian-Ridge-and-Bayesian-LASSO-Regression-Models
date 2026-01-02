# 0. Set seed

set.seed(123)

# Full conditional distribution of beta

sample_beta_lasso <- function(sigma2, tau, X, X_y){
  Z <- solve(X + diag(sigma2/tau, p)) # Compute covariance matrix
  
  mu <- Z%*%X_y # Compute mean vector
  
  beta <- c(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = sigma2*Z)) # Sample beta
  return(beta)
}

# Full conditional distribution of tau2

sample_tau <- function(beta, lambda){
  tau <- GIGrvg::rgig(n = p, lambda = 0.5, chi = beta^2, psi = lambda) # Sample tau2
  return(tau)
}

# Full conditional distribution of sigma2

sample_sigma2_lasso <- function(beta, y, x, e, f){
  residuals <- (y - x%*%beta)
  RSS <- sum(residuals^2) # Compute residual sum of squares
  
  shape <- (0.5*n) + e # Shape parameter
  rate <- (0.5*RSS) + f # Rate parameter
  
  sigma2 <- 1/rgamma(n = 1, shape = shape, rate = rate)
  return(sigma2)
}

# Full conditional distribution of lambda

sample_lambda_lasso <- function(tau, g, h){
  shape <- p + g # Shape parameter
  rate <- (0.5*sum(tau)) + h # Rate parameter
  
  lambda <- rgamma(n = 1, shape = shape, rate = rate)
  return(lambda)
}

# Gibbs sampling algorithm

Gibbs_lasso <- function(y, x, e, f, g, h, n_skip, n_sams, n_burn, verbose = TRUE)
{
  X <- t(x)%*%x # Compute X^{\top} X
  X_y <- t(x)%*%y # Compute X^{\top} y
  
  # Number of iterations of the Gibbs sampling algorithm
  B <- n_burn + n_sams*n_skip
  
  # Objects where the samples of beta, tau, sigma2, and lambda, and log-likelihood will be stored
  BETA <- matrix(data = NA, nrow = n_sams, ncol = p)
  TAU <- matrix(data = NA, nrow = n_sams, ncol = p)
  SIGMA <- matrix(data = NA, nrow = n_sams, ncol = 1)
  LAMBDA <- matrix(data = NA, nrow = n_sams, ncol = 1)
  LL <- matrix(data = NA, nrow = n_sams, ncol = 1)
  
  # Initialize beta, tau, sigma2, and lambda values
  beta <- rep(0, p)
  tau <- rep(1, p)
  sigma2 <- 1
  lambda <- 1
  
    # Gibbs sampling algorithm
    for (i in 1:B) {
      beta <- sample_beta_lasso(sigma2, tau, X, X_y) # Update beta
      tau <- sample_tau(beta, lambda) # Update tau2
      sigma2 <- sample_sigma2_lasso(beta, y, x, e, f) # Update sigma2
      lambda <- sample_lambda_lasso(tau, g, h) # Update lambda
      
      # Save effective samples
      if (i > n_burn && (i - n_burn) %% n_skip == 0) {
        t <- (i - n_burn) / n_skip 
        BETA[t,] <- beta
        TAU[t,] <- tau
        SIGMA[t] <- sigma2
        LAMBDA[t] <- lambda
        LL[t] <- sum(dnorm(x = y, mean = c(x%*%beta), sd = sqrt(sigma2), log = TRUE))
        
      }
  }
  
  return(list(BETA = BETA, SIGMA = SIGMA, TAU = TAU, LAMBDA = LAMBDA, LL = LL))
}
