# 0. Set seed

set.seed(123)

# Full conditional distribution of beta

sample_beta_ridge <- function(sigma2, lambda, X, X_y){
  Z <- solve(X + diag(sigma2*lambda, p)) # Compute covariance matrix
  
  mu <- Z%*%X_y # Compute mean vector
  
  beta <- c(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = sigma2*Z)) # Sample beta
  return(beta)
}

# Full conditional distribution of sigma2

sample_sigma2_ridge <- function(beta, y, x, a, b){
  residuals <- (y - x%*%beta)
  RSS <- sum(residuals^2) # Compute residual sum of squares
  
  shape <- (0.5*n) + a # Shape parameter
  rate <- (0.5*RSS) + b # Rate parameter
  
  sigma2 <- 1/rgamma(n = 1, shape = shape, rate = rate) # Sample sigma2
  return(sigma2)
}

# Full conditional distribution of lambda

sample_lambda_ridge <- function(beta, c, d){
  shape <- (0.5*p) + c # Shape parameter
  rate <- (0.5*(sum(beta^2))) + d # Rate parameter
  
  lambda <- rgamma(n = 1, shape = shape, rate = rate) # Sample lambda
  return(lambda)
}

# Gibbs sampling algorithm

Gibbs_ridge <- function(y, x, a, b, c, d, n_skip, n_sams, n_burn, verbose = TRUE)
{
  X <- t(x)%*%x # Compute X^{\top} X
  X_y <- t(x)%*%y # Compute X^{\top} y
  
  # Number of iterations of the Gibbs sampling algorithm
  B <- n_burn + n_sams*n_skip
  
  # Objects where the samples of beta, sigma2, and lambda, and log-likelihood will be stored
  BETA <- matrix(data = NA, nrow = n_sams, ncol = p)
  SIGMA <- matrix(data = NA, nrow = n_sams, ncol = 1)
  LAMBDA <- matrix(data = NA, nrow = n_sams, ncol = 1)
  LL <- matrix(data = NA, nrow = n_sams, ncol = 1)
  
  # Initialize beta, sigma2, and lambda values
  beta <- rep(0, p)
  sigma2 <- 1
  lambda <- 1
    
    # Gibbs sampling algorithm
    for (i in 1:B) {
      beta <- sample_beta_ridge(sigma2, lambda, X, X_y) # Update beta
      sigma2 <- sample_sigma2_ridge(beta, y, x, a, b) # Update sigma2
      lambda <- sample_lambda_ridge(beta, c, d) # Update lambda
      
      # Save effective samples
      if (i > n_burn && (i - n_burn) %% n_skip == 0) {
        t <- (i - n_burn) / n_skip 
        BETA[t,] <- beta
        SIGMA[t] <- sigma2
        LAMBDA[t] <- lambda
        LL[t] <- sum(dnorm(x = y, mean = c(x%*%beta), sd = sqrt(sigma2), log = TRUE))
      }
      # Algorithm progress
      ncat <- floor(B / 10)
      if (b %% ncat == 0) {
        cat(100 * round(b / B, 1), "% completado ... \n", sep = "")
    }
  }
  return(list(BETA = BETA, SIGMA = SIGMA, LL = LL, LAMBDA = LAMBDA))
}
