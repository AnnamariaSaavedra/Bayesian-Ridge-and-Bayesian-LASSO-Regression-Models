# 0. Set seed

set.seed(123)

# Full conditional distribution of xi

sample_xi_ridge <- function(y, x, xi, beta, sigma2, alpha, a, b, c, d) {
  n <- length(y) # Number of observations
  
  for (i in 1:n) {
    xi_new <- xi
    xi_new[i] <- 0 # Remove cluster assignment of the i-th observation
    
    # Relabel clusters
    xi_unique <- unique(xi_new[-i]) # Number of clusters excluding the i-th observation
    xi_unique <- xi_unique[order(xi_unique)]
    xi_new[-i] <- as.numeric(factor(xi_new[-i], levels = xi_unique, labels = seq_along(xi_unique)))
    xi_unique <- as.numeric(factor(xi_unique))
    
    K_new <- max(xi_unique)  # Number of clusters excluding the i-th observation if all elements are equal to zero
    
    # Active cluster parameters
    beta <- matrix(data = beta[xi_unique, ], nrow = length(xi_unique), ncol = p, byrow = FALSE) # Mean vector
    sigma2 <- sigma2[xi_unique] # Variance parameter
    
    # Vector of explanatory variables, and response variable of the i-th observation
    x_i <- matrix(data = x[i,], nrow = 1, ncol = p, byrow = FALSE)
    y_i <- y[i]
    
    # Object where the probability for the i-th observation to be assigned to a cluster will be stored
    log_probs <- numeric(K_new + 1)
    
    # Probability for the i-th observation to be assigned to an existing cluster k
    for (k in 1:K_new) {
      n_k <- sum(xi_new == k) # Number of observations in the k-th cluster
      log_probs[k] <- log(n_k) + dnorm(y_i, mean = x_i%*%beta[k,], sd = sqrt(sigma2[k]), log = TRUE)
    }
    
    # Sample beta from its prior distribution
    lambda_prior <- rgamma(n = 1, shape = c, rate = d)
    beta_prior <- c(mvtnorm::rmvnorm(n = 1, mean = rep(0, p), sigma = (1/lambda_prior)*(diag(x = 1, p))))
    
    # Sample sigma2 from its prior distribution
    sigma2_prior <- 1/rgamma(n = 1, shape = a, rate = b)
    
    # Probability for the i-th observation to be assigned to a new cluster K + 1
    log_probs[K_new + 1] <- log(alpha) + dnorm(y_i, mean = x_i%*%beta_prior, sd = sqrt(sigma2_prior), log = TRUE)
    
    # Sample xi
    new_cluster <- sample(1:(K_new + 1), size = 1, prob = exp(log_probs - max(log_probs)))
    xi_new[i] <- new_cluster

    # Cluster parameters for the new cluster
    if (new_cluster == (K_new + 1)) {
      
      # Full conditional distribution of sigma2
      shape <- a + 0.5 # Shape parameter
      rate <- b + 0.5*(((y[i] - x_i%*%beta_prior)^2)) # Rate parameter
      
      new_sigma2 <- 1 / rgamma(1, shape = shape, rate = rate)
      
      # Full conditional distribution of beta
      Sigma <- solve((1/sigma2_prior)*(t(x_i)%*%x_i) + (lambda_prior*diag(1, p))) # Covariance matrix
      mu <- Sigma%*%((1/sigma2_prior)*(t(x_i)%*%y[i])) # Mean vector
      
      new_beta <- c(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma))
      
      beta <- rbind(beta, new_beta)
      sigma2 <- c(sigma2, new_sigma2)
      }
    xi <- xi_new
  }
  
  return(list(xi = xi, beta = beta, sigma2 = sigma2))
}

# Full conditional distribution of beta

sample_beta_ridgenp <- function(y, x, xi, sigma2, lambda) 
{
  K <- max(xi) # Number of clusters
  
  # Object where the samples of beta for each cluster will be stored
  beta <- matrix(data = NA, nrow = K, ncol = p)
  
  for (k in 1:K) {
    y_k <- y[xi == k] # Response variables of the k-th cluster
    n_k <- length(y_k) # Number of observations in the k-th cluster
    x_k <- matrix(data = x[xi == k,], nrow = n_k, ncol = p, byrow = FALSE) # Explanatory variables of the k-th cluster
    
    Z <- solve((1/sigma2[k])*(t(x_k)%*%x_k) + (lambda*diag(1, p))) # Compute covariance matrix
    mu <- Z%*%((1/sigma2[k])*(t(x_k)%*%y_k)) # Compute mean vector
    
    beta[k,] <- c(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Z)) # Sample beta
  }
  
  return(beta)
}

# Full conditional distribution of sigma2

sample_sigma2_ridgenp <- function(y, x, xi, beta, a, b)
{
  K <- max(xi) # Number of clusters
  
  # Object where the samples of sigma2 for each cluster will be stored
  sigma2 <- numeric(K)
  
  for (k in 1:K) {
    y_k  <- y[xi == k] # Response variables of the k-th cluster
    n_k  <- length(y_k) # Number of observations in the k-th cluster
    x_k <- matrix(data = x[xi == k,], nrow = n_k, ncol = p, byrow = FALSE) # Explanatory variables of the k-th cluster
    
    beta_k <- beta[k,] # Mean parameter of the k-th cluster
    
    shape <- (0.5*n_k) + a # Shape parameter
    rate <- (0.5*(sum((y_k - x_k%*%beta_k)^2))) + b # Rate parameter
    
    sigma2[k] <- 1/rgamma(1, shape = shape, rate = rate) # Sample sigma2
  }
  
  return(sigma2)
}

# Full conditional distribution of lambda

sample_lambda_ridgenp <- function(xi, beta, c, d) 
{
  K <- length(unique(xi)) # Number of clusters
  
  shape <- (0.5*K*p) + c # Shape parameter
  rate <- (0.5*sum(beta^(2))) + d # Rate parameter
  
  lambda <- rgamma(n = 1, shape = shape, rate = rate) # Sample lambda
  return(lambda)
}

# Full conditional distribution of alpha

sample_alpha <- function(alpha, xi, n, l, m) 
{
  K <- length(unique(xi)) # Number of clusters
  phi <- rbeta(1, shape1 = alpha + 1, shape2 = n)
  pi <- (l + K - 1)/(l + K - 1 + n*(m - log(phi)))
  
  if (runif(1) < pi) {
    return(rgamma(1, shape = l + K, rate = m - log(phi)))
  } else {
    return(rgamma(1, shape = l + K - 1, rate = m - log(phi)))
  }
}

# Gibbs sampling algorithm

Gibbs_ridgenp <- function(y, x, n_burn, n_sams, n_skip, a, b, c, d, l, m, verbose = TRUE) {
  
  max_K <- floor(n / 2) # Maximum number of clusters
  
  # Number of iterations of the Gibbs sampling algorithm
  B <- n_burn + n_sams*n_skip
  ncat <- floor(0.01*B)
  
  # Objects where the samples of xi, beta, sigma2, lambda, and alpha and log-likelihood will be stored
  XI <- matrix(NA, nrow = n_sams, ncol = n)
  BETA <- list(n_sams)
  SIGMA <- vector("list", n_sams)
  LAMBDA <- numeric(n_sams)
  ALPHA <- numeric(n_sams)
  LL <- numeric(n_sams)     
  
  # Initialize xi, beta, sigma2, lambda, and alpha values
  alpha <- rgamma(1, l, m)
  lambda <- rgamma(1, c, d)
  xi <- sample(1:max_K, n, replace = TRUE)
  xi <- as.numeric(factor(xi, levels = unique(xi), labels = seq_along(unique(xi))))
  K <- length(unique(xi)) # Initialize the number of clusters
  sigma2 <- 1/rgamma(K, shape = a, rate = b)
  beta  <- matrix(data = 0, nrow = K, ncol = p)
  
  for (t in 1:B) {
    # Update xi
    xi_b <- sample_xi_ridge(y, x, xi, beta, sigma2, alpha, a, b, c, d)
    xi <- xi_b$xi
    beta <- xi_b$beta
    sigma2 <- xi_b$sigma2
    
    # Update cluster parameters
    beta <- sample_beta_ridgenp(y, x, xi, sigma2, lambda) # Update beta
    sigma2 <- sample_sigma2_ridgenp(y, x, xi, beta, a, b) # Update sigma2
    
    # Update hyperparameters
    lambda <- sample_lambda_ridgenp(xi, beta, c, d) # Update lambda
    alpha <- sample_alpha(alpha, xi, n, l, m) # Update alpha
    
    # Compute log-likelihood
    ll <- 0
    for (j in 1:n) {
      xi_j <- xi[j]
      ll_j <- dnorm(y[j], mean = x[j, ]%*%beta[xi_j,], sd = sqrt(sigma2[xi_j]), log = TRUE)
      ll <- ll + ll_j
    }
    
    # Save effective samples
    if (t > n_burn && (t - n_burn) %% n_skip == 0) {
      i <- (t - n_burn) / n_skip 
      XI[i, ] <- xi
      BETA[[i]] <- beta
      SIGMA[[i]] <- sigma2
      LAMBDA[i] <- lambda
      ALPHA[i] <- alpha
      LL[i] <- ll
     }
    # Algorithm progress
    if (verbose && t %% ncat == 0)
      cat(sprintf("%.1f%% completado\n", 100*t/B))
  }
  
  return(list(XI = XI, BETA = BETA, SIGMA = SIGMA, LAMBDA = LAMBDA, ALPHA = ALPHA, LL = LL))
}
