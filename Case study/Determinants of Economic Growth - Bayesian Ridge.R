# 0. Set seed

rm(list=ls()); set.seed(123)

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(coda)))
suppressMessages(suppressWarnings(library(ggplot2)))

# 2. Import dataset

Data <- read_xlsx(path = "~/Trabajo de grado/Database.xlsx")

# 2.1 Select response variable and explanatory variables

Data <- Data %>%
  select(CODE, # Country code
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
  select(-c(CODE, GR6096)) %>% # Set the matrix containing the explanatory variables
  scale(center = TRUE, scale = TRUE) %>% # Standardize the explanatory variables
  as.matrix()

x <- cbind(x, 1) # Create the intercept column
  
n <- length(y) # Sample size

p <- ncol(x) # Number of explanatory variables

# 3.1 Hyperparameter elicitation

a <- 3 # Shape parameter of inverse-gamma distribution

b <- a*sigma2_OLS # Scale parameter of inverse-gamma distribution

c <- 1 # Shape parameter of gamma distribution

d <- 1 # Rate parameter of gamma distribution

# Before executing the following code, Bayesian Ridge.r must be executed, since it displays
# the Gibbs sampling algorithm

# 3.2 Gibbs sampling algorithm implementation

M2 <- Gibbs_ridge(y, x, a, b, c, d,
                  n_skip = 10, # Accounting for Markov chain autocorrelation will require systematic sampling,
                  n_sams = 10000, # Set the number of effective samples
                  n_burn = 1000) # Set the number of burn-in samples

# Compute the effective sample size for model parameters

TEM_beta <- coda::effectiveSize(M2$BETA); summary(TEM_beta) # beta

TEM_sigma2 <- coda::effectiveSize(M2$SIGMA); summary(TEM_sigma2) # sigma2

TEM_lambda <- coda::effectiveSize(M2$LAMBDA); summary(TEM_lambda) # lambda

# Compute the Monte Carlo standard error for model parameters

EEMC_beta <- apply(X = M2$BETA, MARGIN = 2, FUN = sd)/sqrt(TEM_beta); round(summary(EEMC_beta), 4) # beta

EEMC_sigma2 <- sd(M2$SIGMA)/sqrt(TEM_sigma2); round(summary(EEMC_sigma2), 4) # sigma2

EEMC_lambda <- sd(M2$LAMBDA)/sqrt(TEM_lambda); round(summary(EEMC_lambda), 3) # lambda

# 4. Bayesian inference

# 4.1 Bayesian inference for beta

colnames(M2$BETA) <- paste0("beta", 1:p)

BETA_MEAN <- round(apply(M2$BETA, MARGIN = 2, FUN = mean), 4) # Posterior mean

BETA_MEDIAN <- round(apply(M2$BETA, MARGIN = 2, FUN = median), 4) # Posterior median

BETA_SD <- round(apply(M2$BETA, MARGIN = 2, FUN = sd), 4) # Posterior standard deviation

CI_BETA <- round(apply(M2$BETA, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)), 4) # 95% credible interval

# 4.2 Bayesian inference for sigma2

SIGMA2_MEAN <- round(mean(M2$SIGMA), 5) # Posterior mean

SIGMA2_MEDIAN <- round(median(M2$SIGMA), 5) # Posterior median

SIGMA2_SD <- round(sd(M2$SIGMA), 5) # Posterior standard deviation

CI_SIGMA <- round(quantile(x = M2$SIGMA, probs = c(0.025, 0.975)), 5) # 95% credible interval

# 4.3 Bayesian inference for lambda

LAMBDA_MEAN <- round(mean(M2$LAMBDA), 4) # Posterior mean

LAMBDA_MEDIAN <- round(median(M2$LAMBDA), 4) # Posterior median

LAMBDA_SD <- round(sd(M2$LAMBDA), 4) # Posterior standard deviation

CI_LAMBDA <- round(quantile(x = M2$LAMBDA, probs = c(0.025, 0.975)), 4) # 95% credible interval

# 4.4 Compute information criterion

# Deviance Information Criterion

LL_HAT <- sum(dnorm(x = y, mean = x%*%BETA_MEAN, sd = sqrt(SIGMA2_MEAN), log = TRUE))

LL_B <- M2$LL

pDIC <- 2*(LL_HAT - mean(LL_B))

DIC <- -2*LL_HAT + 2*pDIC

# Watanabe-Akaike Information Criterion

LPPD <- 0

pWAIC <- 0

for (i in 1:n) {
  # LPPD
  TMP <- dnorm(x = y[i], mean = x[i,]%*%t(M2$BETA), sd = sqrt(M2$SIGMA))
  LPPD <- LPPD + log(mean(TMP))
  # pWAIC
  TMP_2 <- dnorm(x = y[i], mean = x[i,]%*%t(M2$BETA), sd = sqrt(M2$SIGMA), log = TRUE)
  pWAIC <- pWAIC + 2*(log(mean(TMP)) - mean(TMP_2))
}

WAIC <- -2*LPPD + 2*pWAIC

# 4.5 Display log-likelihood chain

plot(M2$LL, type = "p", pch = ".", cex = 1.1, col = "#00CD66", xlab = "Iteración", ylab = "Log-verosimilitud", main = "",
     ylim = c(300, 340))
abline(h = mean(M2$LL), lwd = 3, col = "#00CD66")

# 5. Monte Carlo samples from the posterior predictive distribution of test statistics

# Create test statistics function

test_stats <- function(x) {
  c(mean = mean(x),
    median = median(x),
    sd = sd(x),
    iqr = diff(quantile(x, c(0.25, 0.75))),
    min = min(x),
    max = max(x)
  )
}

ts_display <- c("Media", "Mediana", "Desviación estándar", "Rango intercuartílico",
                "Mínimo", "Máximo")

ts <- NULL # Object where the test statistics will be stored

# Simulated statistics

for (b in 1:length(M2$SIGMA)) {
  # Samples from the posterior distribution
  beta <- M2$BETA[b, ]
  sigma2 <- M2$SIGMA[b]
  
  # Posterior predictive datasets, each of size n
  y_tilde <- rnorm(n = n, mean = x%*%beta, sd = sqrt(sigma2))
  
  # Samples from the posterior predictive distribution of test statistics
  ts <- rbind(ts, test_stats(y_tilde)) # Compute test statistics
}

ts_hat <- test_stats(y) # Observed test statistics

# Comparison plots between simulated and observed test statistics

par(mfrow = c(3, 2), mar = c(3, 3, 2, 1), mgp = c(1.75, 0.75, 0))
for (j in 1:length(ts_hat)) {
  test_statistics  <- ts[, j]
  test_statistics_hat <- ts_hat[j]
  
  # Plot histogram
  hist(
    x = test_statistics, freq = FALSE, nclass = 30,
    col = "gray90", border = "gray90",
    xlab = ts_display[j], ylab = "Densidad",
    main = ts_display[j]
  )
  
  abline(
    v = quantile(test_statistics, c(0.025, 0.5, 0.975)),
    col = c(4, 2, 4), lty = c(4, 2, 4), lwd = c(2, 1, 2)
  )
  
  abline(v = test_statistics_hat, lwd = 2)
}

# Compute posterior predictive p-value

ppp <- NULL

for (j in 1:length(ts_hat)) {
  ppp[j] <- round(mean(ts[,j] < ts_hat[j]), 3)
}
