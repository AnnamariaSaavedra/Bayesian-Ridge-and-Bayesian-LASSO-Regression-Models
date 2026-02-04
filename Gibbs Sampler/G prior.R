# 0. Set seed

set.seed(123)

# Monte Carlo algorithm

G_prior <- function(y, x, sigma2_0, nu_0, g, n, p, n_sams, n_skip, n_burn, verbose = TRUE)
{
  # Marginal posterior distribution of sigma2
  shape <- 0.5*(nu_0 + n) # Shape parameter
  
  # Rate parameter
  Matrix <- x%*%solve(t(x)%*%x)%*%t(x)
  RSS <- t(y)%*%(diag(x = 1, n) - (g / (g + 1))*Matrix)%*%y # Residual sum of squares
  rate <- 0.5*((nu_0*sigma2_0) + RSS)
  
  # Full conditional distribution of beta
  Sigma <- (g / (g + 1))*(solve(t(x)%*%x)) # Compute covariance matrix
  mu <- Sigma%*%t(x)%*%y # Compute mean vector
  
  # Number of iterations of the Monte Carlo algorithm
  B <- n_burn + n_sams*n_skip
  ncat <- floor(0.01*B)
  
  # Objects where the samples of beta, and sigma2, and log-likelihood will be stored
  SIGMA <- numeric(n_sams)
  BETA <- matrix(data = NA, nrow = n_sams, ncol = p)
  LL <- numeric(n_sams)
  
  # Monte Carlo algorithm
  for (b in 1:B) {
    sigma2 <- 1/rgamma(n = 1, shape = shape, rate = rate) # Update sigma2
    beta <- c(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = sigma2*Sigma)) # Update beta
    ll <- sum(dnorm(x = y, mean = c(x%*%beta), sd = sqrt(sigma2), log = TRUE)) # Compute log-likelihood
    
    # Save effective samples
    if (b > n_burn && (b - n_burn) %% n_skip == 0) {
      i <- (b - n_burn) / n_skip
      SIGMA[i] <- sigma2
      BETA[i,] <- beta
      LL[i] <- ll
    }
    # Algorithm progress
    if (verbose && b %% ncat == 0)
      cat(sprintf("%.1f%% completado\n", 100*b/B))
  }
  
  return(list(BETA = BETA, SIGMA = SIGMA, LL = LL))
}
