# 0. Set seed

set.seed(123)

# 1. Simulation framework to evaluate Bayesian regression models

n <- c(rep(50, 6), rep(100, 6), rep(200, 6), rep(500, 6)) # Number of observations

K <- c(rep(c(rep(2, 2), rep(3, 2), rep(4, 2)), 4)) # Number of clusters

p <- c(rep(c(10, 20), 12)) # Number of parameters

p_0 <- c(rep(5, 24)) # Number of non-zero elements

Scenery <- data.frame(Escenario = 1:length(n), n = n, K = K, p = p, s = p_0)

sigma2_sim <- list(sigma2_k2 = c(1, 2), sigma2_k3 = c(1, 2, 1), sigma2_k4 = c(1, 2, 1, 3)) # Variance parameter

beta_sim <- matrix(data = c(1.5, 2.5, 3.5, 4.5, 5.5,
                            -1, -2, -3, -4, -5,
                            6, 5, 4, 3, 2,
                            -6.5, -5.5, -4.5, -3.5, -2.5),
                   nrow = max(K), ncol = 5, byrow = TRUE) # Mean parameter

# 2. Data simulation

# Auxiliary functions

# Set the mean vector for the i-th scenery

mean_parameter <- function(beta, beta_sim, K, p, p_0){
  if (K == 2) {
    for (i in 1:2) {
      for (j in 1:p_0) {
        beta[i, j] <- beta_sim[i, j]
      }
    }
  } else if (K == 3) {
    for (i in 1:2) {
      for (j in 1:p_0) {
        beta[i, j] <- beta_sim[i, j]
      }
    }
    for (j in 1:p_0) {
      beta[3, p_0 + j] <- beta_sim[3,j]
    }
  } else{
    for (i in 1:2) {
      for (j in 1:p_0) {
        beta[i, j] <- beta_sim[i, j]
      }
    }
    for (i in 3:4) {
      for (j in 1:p_0) {
        beta[i, p_0 + j] <- beta_sim[i,j]
      }
    }
  }
  return(beta)
}

# Set the variance parameter for the i-th scenery

variance_parameter <- function(sigma2_sim, K){
  if (K == 2) {
    sigma2 <- sigma2_sim[[1]]
  } else if (K == 3) {
    sigma2 <- sigma2_sim[[2]]
  } else{
    sigma2 <- sigma2_sim[[3]]
  }
  return(sigma2)
}

# Compute Deviance Information Criterion

compute_DIC <- function(model, y, x, beta, sigma2){
  LL_HAT <- sum(dnorm(x = y, mean = x%*%beta, sd = sqrt(sigma2), log = TRUE))
  LL_B <- model$LL
  
  pDIC <- 2*(LL_HAT - mean(LL_B))
  DIC <- -2*LL_HAT + 2*pDIC
  
  return(DIC)
}

# Compute Watanabe-Akaike Information Criterion

compute_WAIC <- function(model, y, x, n){
  LPPD <- 0
  pWAIC <- 0
  
  for (i in 1:n) {
    # LPPD
    TMP <- dnorm(x = y[i], mean = x[i,]%*%t(model$BETA), sd = sqrt(model$SIGMA))
    LPPD <- LPPD + log(mean(TMP))
    # pWAIC
    TMP_2 <- dnorm(x = y[i], mean = x[i,]%*%t(model$BETA), sd = sqrt(model$SIGMA), log = TRUE)
    pWAIC <- pWAIC + 2*(log(mean(TMP)) - mean(TMP_2))
  }
  WAIC <- -2*LPPD + 2*pWAIC
  
  return(WAIC)
}

# Bayesian inference in finite mixture models

function(model, B, K, p){
  permu <- gtools::permutations(n = K, r = K) # Permutations
  
  # Average over the permuted spaces
  # Objects where the samples of beta, and sigma2 will be stored
  BETA_corrected <- vector("list", p)
  SIGMA2_corrected <- vector("list")
  
  # Posterior distribution of K
  k <- apply(X = model$XI, MARGIN = 1, function(x) length(unique(x)))
  k_tab <- table(factor(x = k, levels = k, labels = k))
  
  beta_pos <- matrix(0, nrow = K, ncol = p)
  sigma2_pos <- 0
  
  for (b in 1:B) {
    if (length(unique(model$XI[b,])) == K) {
      beta_pos <- beta_pos +
        model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE] / K_tab[1] # Posterior mean
      
      sigma_pos <- sigma_pos + model$SIGMA[[b]] / K_tab[1] # Posterior mean
    }
  }
  
  for (b in 1:B) {
    if (length(table(model$XI[b,])) == K) {
      for (j in 1:p) {
        beta_current <- model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE][,j]
        # Reorder according to the permutations, and compute the distance of each sample to its posterior mean
        dist <- apply(X = permu, 
                      MARGIN = 1, 
                      FUN = function(p) {
                        permuted_beta <- beta_current[p]
                        sum((permuted_beta - beta_pos[,j])^2)
                      }
        )
        # Select the optimum permutation
        best_permu <- permu[which.min(dist),]
        BETA_corrected[[j]] <- rbind(BETA_corrected[[j]], beta_current[best_permu]) 
      }
      sigma2_current <- model$SIGMA[[b]]
      # Reorder according to the permutations, and compute the distance of each sample to its posterior mean
      dist <- apply(X = permu, 
                    MARGIN = 1, 
                    FUN = function(p) {
                      permuted_sigma2 <- sigma2_current[p]
                      sum((permuted_sigma2 - sigma2_pos)^2)
                    }
      )
      # Select the optimum permutation
      best_permu <- permu[which.min(dist),]
      SIGMA2_corrected <- rbind(SIGMA2_corrected, sigma2_current[best_permu])
    }
  }
  # Posterior mean, and 95% credible interval
  beta <- matrix(data = NA, nrow = (p + 1), ncol = 3*K)
  if (K == 2) {
    for (i in 1:p) {
      beta[i,1] <- colMeans(BETA_corrected[[i]])[1]
      beta[i,2] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,1]
      beta[i,3] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,1]
      beta[i,4] <- colMeans(BETA_corrected[[i]])[2]
      beta[i,5] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,2]
      beta[i,6] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,2]
    }
    beta[(p + 1),1] <- mean(SIGMA2_corrected[[1]])
    beta[(p + 1),2] <- quantile(SIGMA2_corrected[[1]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),3] <- quantile(SIGMA2_corrected[[1]], probs = c(0.025, 0.975))[2]
    beta[(p + 1),4] <- mean(SIGMA2_corrected[[2]])
    beta[(p + 1),5] <- quantile(SIGMA2_corrected[[2]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),6] <- quantile(SIGMA2_corrected[[2]], probs = c(0.025, 0.975))[2]  
  } else if (K == 3){
    for (i in 1:p) {
      beta[i,1] <- colMeans(BETA_corrected[[i]])[1]
      beta[i,2] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,1]
      beta[i,3] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,1]
      beta[i,4] <- colMeans(BETA_corrected[[i]])[2]
      beta[i,5] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,2]
      beta[i,6] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,2]
      beta[i,7] <- colMeans(BETA_corrected[[i]])[3]
      beta[i,8] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,3]
      beta[i,9] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,3]
    }
    beta[(p + 1),1] <- mean(SIGMA2_corrected[[1]])
    beta[(p + 1),2] <- quantile(SIGMA2_corrected[[1]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),3] <- quantile(SIGMA2_corrected[[1]], probs = c(0.025, 0.975))[2]
    beta[(p + 1),4] <- mean(SIGMA2_corrected[[2]])
    beta[(p + 1),5] <- quantile(SIGMA2_corrected[[2]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),6] <- quantile(SIGMA2_corrected[[2]], probs = c(0.025, 0.975))[2]
    beta[(p + 1),7] <- mean(SIGMA2_corrected[[3]])
    beta[(p + 1),8] <- quantile(SIGMA2_corrected[[3]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),9] <- quantile(SIGMA2_corrected[[3]], probs = c(0.025, 0.975))[2]
  } else{
    for (i in 1:p) {
      beta[i,1] <- colMeans(BETA_corrected[[i]])[1]
      beta[i,2] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,1]
      beta[i,3] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,1]
      beta[i,4] <- colMeans(BETA_corrected[[i]])[2]
      beta[i,5] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,2]
      beta[i,6] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,2]
      beta[i,7] <- colMeans(BETA_corrected[[i]])[3]
      beta[i,8] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,3]
      beta[i,9] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,3]
      beta[i,10] <- colMeans(BETA_corrected[[i]])[4]
      beta[i,11] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[1,4]
      beta[i,12] <- apply(BETA_corrected[[i]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))[2,4]
    }
    beta[(p + 1),1] <- mean(SIGMA2_corrected[[1]])
    beta[(p + 1),2] <- quantile(SIGMA2_corrected[[1]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),3] <- quantile(SIGMA2_corrected[[1]], probs = c(0.025, 0.975))[2]
    beta[(p + 1),4] <- mean(SIGMA2_corrected[[2]])
    beta[(p + 1),5] <- quantile(SIGMA2_corrected[[2]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),6] <- quantile(SIGMA2_corrected[[2]], probs = c(0.025, 0.975))[2]
    beta[(p + 1),7] <- mean(SIGMA2_corrected[[3]])
    beta[(p + 1),8] <- quantile(SIGMA2_corrected[[3]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),9] <- quantile(SIGMA2_corrected[[3]], probs = c(0.025, 0.975))[2]
    beta[(p + 1),10] <- mean(SIGMA2_corrected[[3]])
    beta[(p + 1),11] <- quantile(SIGMA2_corrected[[3]], probs = c(0.025, 0.975))[1]
    beta[(p + 1),12] <- quantile(SIGMA2_corrected[[3]], probs = c(0.025, 0.975))[2]
  }
}

# Deviance Information Criterion

beta <- beta[,c(1,4)] # Estimacion puntual corregida

LL_HAT <- NA

for (i in 1:n) {
  LL_HAT[i] <- dnorm(x = y[i], mean = x[i,]%*%beta[,data_2$xi[i]], sd = sqrt(sigma2[data_2$xi[i]]), log = TRUE) 
}

LL_B <- model$LL

pDIC <- 2*(sum(LL_HAT) - mean(LL_B))

DIC <- -2*sum(LL_HAT) + 2*pDIC

# Watanabe-Akaike Information Criterion

B <- 10000; n <- length(data_1$y)

TMP <- matrix(data = NA, nrow = B, ncol = n)

TMP_2 <- matrix(data = NA, nrow = B, ncol = n)

for (i in 1:n) {
  for (b in 1:B) {
    xi <- model$XI
    beta <- model$BETA[[b]][xi[b,i], , drop = FALSE]
    sigma2 <- model$SIGMA[[b]][xi[b,i]]
    
    # LPPD
    TMP[b,i] <- dnorm(x = y[i], mean = x[i,]%*%t(beta), sd = sqrt(sigma2))
    
    #pWAIC
    TMP_2[b,i] <- dnorm(x = y[i], mean = x[i,]%*%t(beta), sd = sqrt(sigma2), log = TRUE) 
  }
}

LPPD <- sum(log(colMeans(TMP)))

pWAIC <- 2*sum((log(colMeans(TMP)) - colMeans(TMP_2)))

WAIC <- -2*LPPD + 2*pWAIC

# 3. Create dataset for each scenary

create_dataset <- function(Scenary) {
  # Objects where the response variable, explanatory variables, density function, and cluster indicator variable will be stored
  Y <- list()
  X <- list()
  F_TRUE <- list()
  XI <- list()
  
  BETA <- list()
  SIGMA2 <- list()
  
  # Simulated dataset for each scenery
  for (i in 1:nrow(Scenery)) {
    n <- Scenery$n[i] # Number of observations
    K <- Scenery$K[i] # Number of clusters
    p <- Scenery$p[i] # Number of parameters
    p_0 <- Scenery$s[i] # Number of non-zero elements
    
    omega <- rep(1/K, K) # Vector of mixing proportions
    xi <- sample(x = 1:K, size = n, replace = TRUE, prob = omega) # Vector of cluster assignment
    
    n_k <- table(xi) # Size of cluster k
    
    x <- matrix(data = NA, nrow = n, ncol = p) # Matrix of explanatory variables
    
    x[,1] <- rep(1, n) # Intercept column
    x[,2] <- sample(x = 0:1, size = n, replace = TRUE, prob = c(0.5, 0.5)) # Dummy variable
    
    lambda <- c(1, 3, 9) # Rate parameter of the Poisson distribution
    
    for (j in 3:5) {
      x[,j] <- rpois(n = n, lambda = lambda[j - 2]) # Poisson-distributed variables
    }
    
    if (p == 10) {
      mu <- seq(from = 1, to = 5, by = 1) # Mean parameter
      sigma2 <- seq(from = 1, to = 5, by = 1) # Variance parameter
      
      for (j in 6:p) {
        x[,j] <- rnorm(n = n, mean = mu[j - 5], sd = sqrt(sigma2[j - 5])) # Normal-distributed variables
      }
    } else if (p == 20) {
      mu <- seq(from = -9, to = 5, by = 1) # Mean parameter
      sigma2 <- rep(c(seq(from = 1, to = 3, by = 1)), 5) # Variance parameter
      
      for (j in 6:p) {
        x[,j] <- rnorm(n = n, mean = mu[j - 5], sd = sqrt(sigma2[j - 5])) # Normal-distributed variables
      }
    } else {
      mu <- seq(from = -9, to = 10, by = 1) # Mean parameter
      sigma2 <- rep(c(seq(from = 1, to = 5, by = 1)), 4) # Variance parameter
      
      for (j in 6:p) {
        x[,j] <- rnorm(n = n, mean = mu[j - 5], sd = sqrt(sigma2[j - 5])) # Normal-distributed variables
      }}
    
    # Save the matrix of explanatory variables
    X[[i]] <- x
    
    beta <- matrix(data = rep(0, K*p), nrow = K, ncol = p) # Mean vector
    
    beta <- mean_parameter(beta, beta_sim, K, p, p_0) # Mean vector
    sigma2 <- variance_parameter(sigma2_sim, K) # Variance parameter
    y <- numeric() # Response variable
    f_true <- numeric() # Density function
    
    # Sample response variable for each cluster
    for (k in 1:K) {
    y[xi == k] <- rnorm(n = n_k[k], mean = x[xi == k,]%*%beta[k,], sd = sqrt(sigma2[k]))
    f_true[xi == k] <- dnorm(x = y[xi == k], mean = x[xi == k,]%*%beta[k,], sd = sqrt(sigma2[k]))
    }
    
    # Save the response variable, the density function, and the cluster indicator variables
    Y[[i]] <- y
    F_TRUE[[i]] <- f_true
    XI[[i]] <- xi
    
    BETA[[i]] <- beta
    SIGMA2[[i]] <- sigma2
  }
  return(list(X = X, Y = Y, F_TRUE = F_TRUE, XI = XI, BETA = BETA, SIGMA2 = SIGMA2))
}

data <- create_dataset(Scenary)

save(data, file = "~/Downloads/Dataset - Simulation framework")

