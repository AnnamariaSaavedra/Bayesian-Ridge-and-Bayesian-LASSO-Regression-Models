# 0. Set seed

rm(list=ls()); set.seed(123)

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(coda)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(corrplot)))
suppressMessages(suppressWarnings(library(aricode)))
suppressMessages(suppressWarnings(library(GIGrvg)))

# 2. Import dataset

Data <- read_xlsx(path = "~/Trabajo de grado/Database.xlsx")

# 2.1 Select response variable and explanatory variables

Data <- Data %>%
  dplyr::select(CODE, # Country code
         GR6096, # Average growth rate of Gross Domestic Product (GDP) per capita between 1960 and 1996
         GDPCH60L, # GDP per capita in 1960 (logaritmic scale)
         LIFE060, # Life expectancy in 1960
         P60, # Primary school enrollment rate in 1960
         SAFRICA, # Dummy variable for Sub-Sahara Africa
         LAAM, # Dummy variable for Latin America
         EUROPE, # Dummy variable for Europe
         EAST, # Dummy variable for Asia
         SCOUT, # Dummy variable for Outward orientation
         DPOP6090, # Growth rate of population between 1960 and 1990
         H60, # Higher education enrollment rate in 1960
         YRSOPEN, # Number of years an economy has been open between 1950 and 1994
         REVCOUP, # Number of military coups and revolutions
         WARTORN, # Dummy variable for countries that have been involved in war any time between 1960 and 1990
         PRIGHTS, # Index of political rights
         CIV72, # Index of civil liberties
         ABSLATIT, # Absolute latitude
         AVELF, # Index of ethnolinguistic fractionalization
         PRIEXP70, # Fraction of primary exports in total exports in 1970
         RERD, # Real exchange rate distortions
         BRIT, # Dummy variable for former British colonies
         SPAIN, # Dummy variable for former Spanish colonies
         BUDDHA, # Percentage of the population that follows the Buddhist religion
         CATH00, # Percentage of the population that follows the Catholic religion
         CONFUC, # Percentage of the population that follows Confucian religion
         HINDU00, # Percentage of the population that follows the Hindu religion
         MUSLIM00, # Percentage of the population that follows the Muslim religion
         PROT00, # Percentage of the population that follows the Protestant religion
         MINING, # Percentage of the GDP in mining
         ECORG, # Index of degree in which economies favor capitalist forms of production
         OTHFRAC, # Percentage of the population speaking foreign language
         ENGFRAC # Percentage of the population speaking english
  )

# 2.2 Remove observations with some missing values

Data <- Data[complete.cases(Data), ]

# 3. Bayesian Ridge regression

y <- Data$GR6096 # Set the response variable

x <- Data %>%
  dplyr::select(-c(CODE, GR6096)) %>% # Set the matrix containing the explanatory variables
  scale(center = TRUE, scale = TRUE) %>% # Standardize the explanatory variables
  as.matrix()

x <- cbind(x, 1) # Create the intercept column

n <- length(y) # Sample size

p <- ncol(x) # Number of explanatory variables

# 3.1 Hyperparameter elicitation

e <- 3 # Shape parameter of inverse-gamma distribution

f <- e*sigma2_OLS # Scale parameter of inverse-gamma distribution

g <- 1 # Shape parameter of gamma distribution

h <- 1 # Rate parameter of gamma distribution

l <- 1 # Shape parameter of gamma distribution

m <- 1 # Shape parameter of gamma distribution

# Before executing the following code, Bayesian Nonparametric LASSO.r must be executed, since it displays
# the Gibbs sampling algorithm

# 3.2 Gibbs sampling algorithm implementation

M5 <- Gibbs_lassonp(y,
                    n_burn = 1000, # Set the number of burn-in samples
                    n_sams = 10000, # Set the number of effective samples
                    n_skip = 10, # Accounting for Markov chain autocorrelation will require systematic sampling
                    e, f, g, h, l, m)

# 4. Bayesian inference

# 4.1 Display the log-likelihood chain

plot(M5$LL, type = "p", pch = ".", cex = 1.1, col = "deeppink3", xlab = "Iteración", ylab = "Log-verosimilitud", main = "",
     ylim = c(300, 340))
abline(h = mean(M5$LL), lwd = 3, col = "deeppink3")

# 4.2 Inference on the number of clusters

K <- apply(M5$XI, 1, function(x) length(unique(x))) # Compute the number of clusters at each iteration

K_table <- as.data.frame(table(K)/length(K))

# Plot the posterior distribution of K

ggplot(K_table, aes(x = K, y = Freq)) +
  geom_segment(aes(x = K, xend = K, y = 0, yend = Freq),
               color = "cyan2", lwd = 1.5) +
  ylim(c(0, 1)) +
  ylab("Densidad") + xlab("Número de clústeres") +
  theme_bw(base_size = 16) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# 4.3 Inference for Bayesian nonparametric Ridge

hat_np <- function(model, K, p, ite){
  permu <- gtools::permutations(n = K, r = K) # Permutations
  
  # Objects where the samples of beta, and sigma2 will be stored
  BETA_corrected <- vector("list", p)
  SIGMA2_corrected <- NULL
  
  # Posterior distribution of K
  k <- apply(X = model$XI, MARGIN = 1, function(x) length(unique(x)))
  k_tab <- table(factor(x = k, levels = k, labels = k))
  
  # Objects where the posterior mean of beta, and sigma2 will be stored
  beta_pos <- matrix(0, nrow = K, ncol = p)
  sigma2_pos <- 0
  
  for (b in 1:ite) {
    if (length(unique(model$XI[b,])) == K) {
      beta_pos <- beta_pos +
        model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE] / max(k_tab)
      
      sigma2_pos <- sigma2_pos + model$SIGMA[[b]] / max(k_tab)
    }
  }
  
  # Average over the permuted spaces
  for (b in 1:ite) {
    if (length(table(model$XI[b,])) == K) {
      for (j in 1:p) {
        beta_current <- model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE][,j]
        # Reorder according to the permutations, and compute the distance of each sample to its posterior mean
        dist <- apply(X = permu, MARGIN = 1, 
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
      dist <- apply(X = permu, MARGIN = 1, 
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
  inf <- matrix(data = NA, nrow = p + 1, ncol = 5*K)
  id <- function(k)(5*(k - 1) + 1):(5*k)
  
  for (j in 1:p) {
    beta_mean <- colMeans(BETA_corrected[[j]])
    beta_median <- apply(BETA_corrected[[j]], MARGIN = 2, FUN = median)
    beta_sd <- apply(BETA_corrected[[j]], MARGIN = 2, FUN = sd)
    beta_ic <- apply(X = BETA_corrected[[j]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
    
    for (k in 1:K) {
      inf[j, id(k)] <- c(beta_mean[k], beta_median[k], beta_sd[k], beta_ic[1, k], beta_ic[2, k])
    }
  }
  sigma2_mean <- colMeans(SIGMA2_corrected)
  sigma2_median <- apply(SIGMA2_corrected, MARGIN = 2, FUN = median)
  sigma2_sd <- apply(SIGMA2_corrected, MARGIN = 2, FUN = sd)
  sigma2_ic <- apply(SIGMA2_corrected, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
  
  for (k in 1:K) {
    inf[p + 1, id(k)] <- c(sigma2_mean[k], sigma2_median[k], sigma2_sd, sigma2_ic[1, k], sigma2_ic[2, k])
  }
  return(list(inference = inf))
}

inference <- hat_np(M5, K = 1, p, ite = 10000)$inference

# Bayesian inference for lambda

LAMBDA_MEAN <- round(mean(M5$LAMBDA), 4) # Posterior mean

LAMBDA_MEDIAN <- round(median(M5$LAMBDA), 4) # Posterior median

LAMBDA_SD <- round(sd(M5$LAMBDA), 4) # Posterior standard deviation

CI_LAMBDA <- round(quantile(x = M5$LAMBDA, probs = c(0.025, 0.975)), 4) # 95% credible interval

# 4.3.4 Bayesian inference for alpha

ALPHA_MEAN <- round(mean(M5$ALPHA), 4) # Posterior mean

ALPHA_MEDIAN <- round(median(M5$ALPHA), 4) # Posterior median

ALPHA_SD <- round(sd(M5$ALPHA), 4) # Posterior standard deviation

CI_ALPHA <- round(quantile(x = M5$ALPHA, probs = c(0.025, 0.975)), 4) # 95% credible interval

# 4.4 Compute information criterion

# Deviance Information Criterion

# Posterior mean of beta, and sigma 2 corrected

beta_mean <- round(inference[-nrow(inference),1], 4)
sigma2_mean <- round(inference[nrow(inference),1], 5)

LL_HAT <- sum(dnorm(x = y, mean = x%*%beta_mean, sd = sqrt(sigma2_mean), log = TRUE))

LL_B <- M5$LL

pDIC <- 2*(LL_HAT - mean(LL_B))

DIC <- -2*LL_HAT + 2*pDIC

# Watanabe-Akaike Information Criterion

compute_WAIC_np <- function(model, y, x, n, ite){
  TMP <- matrix(data = NA, nrow = ite, ncol = n)
  TMP_2 <- matrix(data = NA, nrow = ite, ncol = n)
  
  for (i in 1:n) {
    for (b in 1:ite) {
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
  
  return(WAIC)
}

WAIC <- compute_WAIC_np(M5, y, x, n, ite = 10000)

# 4.5 Bayesian inference for density function

posterior_density_estimate <- function(model, y_seq, x) {
  B <- length(model$BETA)
  n <- ncol(model$XI)
  M <- length(y_seq)
  
  # Store density estimates
  FE <- matrix(NA, nrow = B, ncol = M)
  
  for (b in 1:B) {
    xi <- model$XI[b, ]
    beta <- model$BETA[[b]]
    sigma2 <- model$SIGMA[[b]]
    
    # Unique cluster labels
    cluster_labels <- sort(unique(xi))
    cluster_counts <- as.numeric(table(xi)[as.character(cluster_labels)])
    
    # Cluster probabilities
    weight_b <- cluster_counts / sum(cluster_counts)
    
    for (i in 1:M) {
      FE[b, i] <- sum(weight_b * dnorm(y_seq[i], mean = x%*%t(beta), sd = sqrt(sigma2)))
    }
  }
  
  f_hat <- colMeans(FE)  # Posterior mean density
  f_inf <- apply(FE, 2, quantile, probs = 0.025)  # 2.5%  credible interval
  f_sup <- apply(FE, 2, quantile, probs = 0.975)  # 97.5% credible interval
  
  return(list(f_hat = f_hat, f_inf = f_inf, f_sup = f_sup))
}

# Define a sequence of y values for density estimation
y_seq <- seq(min(y), max(y), length.out = 150)

# Define a grid of x values
x_seq <- colMeans(x)

# Compute posterior density estimate and credible intervals
density_estimate <- posterior_density_estimate(M5, y_seq, x_seq)

f_hat <- density_estimate$f_hat
f_inf <- density_estimate$f_inf
f_sup <- density_estimate$f_sup

# Plot the histogram
hist(x = y, freq = FALSE, xlim = c(-0.06, 0.08), ylim = c(0, 50),
     ylab = "Densidad", main = "",
     col = alpha("grey", 0.3), cex.lab = 1.5, cex.axis  = 1.5)
# Overlay the posterior density estimate as a blue line
polygon(c(y_seq, rev(y_seq)), c(f_inf, rev(f_sup)), col = alpha("deeppink3", 0.3), border = NA)
lines(y_seq, f_hat, lwd = 2, col = "deeppink3")
