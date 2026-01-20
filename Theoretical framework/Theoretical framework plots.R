# 0. Set seed

set.seed(1)

# Figure 1----

# Figure 1 displays the estimated coefficient \beta under Ordinary Least Squares
# and G prior.

# 1. Load necessary libraries
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(coda)))
suppressMessages(suppressWarnings(library(ISLR2)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(lmtest)))
suppressMessages(suppressWarnings(library(sandwich)))

# 2. Import dataset from ISLR2 library

Hitters <- na.omit(Hitters) # Remove observations with some missing values

y <- Hitters$Salary # Set the response variable

x <- model.matrix(Salary ~ ., Hitters)[, -1] # Set the matrix containing the explanatory variables

x <- cbind(1, x) # Create the intercept column

n <- length(y) # Sample size

p <- ncol(x) # # Number of explanatory variables

# 3. Display OLS estimation

OLS <- lm(y ~ x - 1)

beta <- coeftest(OLS, vcov. = vcovHC, type = "HC3")[,1] # Heteroskedasticity-robust standard errors

CI <- round(coefci(OLS, vcov. = vcovHC, type = "HC3"), 5) # 95% confidence interval

# 4. Hyperparameter elicitation

beta_OLS <- solve(t(x)%*%x)%*%t(x)%*%y

residuals <- y - x%*%beta_OLS

sigma2_OLS <- sum(residuals^2)/(n - p)

nu_0 <- 1

sigma2_0 <- sigma2_OLS

g <- n # It reflects information equivalent to that provided by an observation

# Before executing the following code, G prior.r must be executed, since it displays
# the Monte Carlo algorithm

# 5. Monte Carlo algorithm implementation

M1 <- G_prior(shape, rate, mu, Sigma, 
              n_sams = 10000, # Set the number of effective samples
              n_skip = 1, # Accounting for Markov chain autocorrelation will require systematic sampling
              n_burn = 1000) # Set the number of burn-in samples

# 6. Bayesian inference for beta

BETA_HAT <- apply(M1$BETA, MARGIN = 2, FUN = mean) # Posterior mean

CI <- apply(M1$BETA, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)) # 95% credible interval

# Display Bayesian estimation

Data <- data.frame(Position = 1:p, Estimate = BETA_HAT, Lower = CI[1, ],
                   Upper = CI[2, ], Frecuentist = beta_OLS)

Plot_Bayesian <- ggplot(Data, aes(x = 1:p)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1) +
  geom_point(aes(y = Estimate, color = "Bayesiano"), shape = 16, size = 2) +
  geom_point(aes(y = Frecuentist, color = "Frecuentista"), shape = 15, size = 1) +
  labs(
    main = "Modelo de regresión normal Bayesiano",
    sub = "Previa g",
    x = "Coeficientes estimados",
    y = expression(beta[j]),
    color = ""
  ) +
  scale_x_continuous(breaks = 1:p) +
  scale_y_continuous(breaks = seq(-200, 400, by = 100)) + 
  scale_color_manual(
    values = c("Bayesiano" = "gold1", "Frecuentista" = "firebrick1"),
    labels = c("Bayesiano" = "Media posterior", "Frecuentista" = "OLS")
  ) +
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='bottom', legend.text = element_text(size = 14))

# Figure 2----

# Figure 2 displays the estimated coefficient \beta for the shrinkage parameter 
# \lambda in the Ridge regression model.

# 1. Load necessary library
suppressMessages(suppressWarnings(library(ISLR2)))
suppressMessages(suppressWarnings(library(glmnet)))

# 2. Import dataset Hitters
Hitters <- na.omit(Hitters)
Hitters <- Hitters[, c(1:7, 19)]

x <- model.matrix(Salary ~ ., Hitters)[, -1] # Qualitative variables are transformed into dummy variables
y <- Hitters$Salary # Response varibale

# 3. Set the grid for the shrinkage parameter
grid <- 10^seq(10, -2, length = 100)

# 4. Estimate the Ridge regression model
ridge <- glmnet(x, y, alpha = 0, lambda = grid) # The glmnet function standardises explanatory variables

# 5. Plot the estimated coefficient \beta
plot(ridge, xvar = "lambda", col = rainbow(8)[-4], ylab = "Coeficiente estimado",
     cex.lab = 1.25)
legend("topright", legend = c(expression(hat(beta)[1]), 
                              expression(hat(beta)[2]),
                              expression(hat(beta)[3]),
                              expression(hat(beta)[4]),
                              expression(hat(beta)[5]),
                              expression(hat(beta)[6]),
                              expression(hat(beta)[7])), 
       col = c(rainbow(8)[1], rainbow(8)[2], rainbow(8)[3],
               rainbow(8)[5], rainbow(8)[6], rainbow(8)[7],
               rainbow(8)[8]), 
       lty = 1, bty = "n", cex = 1)

# Figure 4----

# Figure 4 displays the estimated coefficient \beta for the shrinkage parameter 
# \lambda in the LASSO regression model.

# 1. Load necessary library
suppressMessages(suppressWarnings(library(ISLR2)))
suppressMessages(suppressWarnings(library(glmnet)))

# 2. Import dataset Hitters
Hitters <- na.omit(Hitters)
Hitters <- Hitters[, c(1:7, 19)]

x <- model.matrix(Salary ~ ., Hitters)[, -1] # Qualitative variables are transformed into dummy variables
y <- Hitters$Salary # Response varibale

# 3. Set the grid for the shrinkage parameter
grid <- 10^seq(10, -2, length = 100)

# 4. Estimate the Ridge regression model
lasso <- glmnet(x, y, alpha = 1, lambda = grid) # The glmnet function standardises explanatory variables

# 5. Plot the estimated coefficient \beta
plot(lasso, xvar = "lambda", col = rainbow(8)[-4], ylab = "Coeficiente estimado",
     cex.lab = 1.25)
legend("topright", legend = c(expression(hat(beta)[1]), 
                              expression(hat(beta)[2]),
                              expression(hat(beta)[3]),
                              expression(hat(beta)[4]),
                              expression(hat(beta)[5]),
                              expression(hat(beta)[6]),
                              expression(hat(beta)[7])), 
       col = c(rainbow(8)[1], rainbow(8)[2], rainbow(8)[3],
               rainbow(8)[5], rainbow(8)[6], rainbow(8)[7],
               rainbow(8)[8]), 
       lty = 1, bty = "n", cex = 1)

# Figure 6----

# Figure 6 corresponds to realizations from a $\textsf{DP}(\alpha, G_0$), 
# for $\alpha \in \{ 1, 10, 100 \}$ and $G_0 = \textsf{N}(0,1)$. 
# The solid black line corresponds to $G_0$, while the color lines represent the 
# realizations.

# 1. Load necessary library
suppressMessages(suppressWarnings(library(gtools)))

# 2. Set parameters
n <- 1000  # Number of points
alpha <- c(1, 10, 100)  # $\alpha$ values
G0 <- function(x) pnorm(x)  # Base measure: Standard normal distribution

# 3. Set graphical parameters
par(mfrow = c(1, 3), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

# 4. Simulation for different alpha values
for (i in alpha) {
  # Set up plot
  plot(NA, NA, xlim = c(-3, 3), ylim = c(0, 1), 
       xlab = "x", ylab = "G(x)", 
       main = bquote(i == .(i)))
  
  # Generate 10 realizations
  for (l in 1:10) {
    # Sample x from a uniform distribution x ~ U(-3, 3) and sort x values
    x <- sort(runif(n = n, min = -3, max = 3))
    
    # Compute concentration parameters
    a <- numeric(n + 1)
    a[1] <- i * G0(x[1]) # First component of the parameter vector from the Dirichlet distribution
    a[n + 1] <- i * (1 - G0(x[n])) # Last component of the parameter vector from the Dirichlet distribution
    
    for (j in 2:n) {
      a[j] <- i * (G0(x[j]) - G0(x[j - 1]))
    }
    
    # Sample from Dirichlet distribution
    u <- gtools::rdirichlet(n = 1, alpha = a)
    
    # Plot cumulative sum of sampled weights
    lines(x, cumsum(u)[-length(u)], type = "l", col = which(alpha == i))
  }
  
  # Add base measure curve
  curve(G0, from = -3, to = 3, n = 1000, lwd = 2, add = TRUE)
}


# Figure 9----

# Figure 9 corresponds to a realization from a $\textsf{DP}(\alpha, G_0$)
# for $\alpha = 25$ and $G_0 = \textsf{N}(0,1)$. The solid black line corresponds to $G_0$,
# while the blue line represents the realization.

# The spiked lines are located at 1000 draws from $G_0$ with heights given by the 
# stick-breaking weights.

# 1. Load necessary library
suppressMessages(suppressWarnings(library(stats)))

# 2. Set parameters
alpha <- 25
epsilon <- 1e-6
J <- round(log(epsilon) / log(alpha / (alpha + 1)))
G0 <- function(x) pnorm(x) # Base measure: Standard normal distribution

# 3. Construct atoms and weights using stick-breaking process
vartheta <- rnorm(n = J)
V <- rbeta(n = J, shape1 = 1, shape2 = alpha)

# 4. Compute stick-breaking weights
p <- numeric(J)

p[1] <- V[1] # p_{1} = V_{1}

for (j in 2:(J - 1)) {
  p[j] <- V[j] * prod(1 - V[1:(j - 1)]) # p_{k - 1} = (1 - \sum_{j = 1}^{k - 2} p_{j}) V_{k - 1}
}

p[J] <- 1 - sum(p[1:(J - 1)]) # p_{k} = 1 - \sum_{j = 1}^{J - 1} p_{j}

# 5. Set graphical parameters
par(mfrow = c(1, 2), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

# Panel 1: Plot atoms and weights
plot(x = vartheta, y = p, type = "h", xlim = c(-3, 3), 
     xlab = "x", ylab = expression(p), col ="deepskyblue", 
     main = "")

# Panel 2: Plot the cumulative distribution function
x <- seq(from = -3, to = 3, length.out = 1000) # Sequence of x values
G <- sapply(x, function(x) sum(p[vartheta <= x]))

plot(x, G, type = "l", xlim = c(-3, 3), col = "darkorchid2", 
     xlab = "x", ylab = "G(x)", main = "")
curve(G0, from = -3, to = 3, n = 1000, lwd = 2, add = TRUE) # Add the base measure G0 as reference

# Figure 9----

# Figure 9 corresponds to a realization from a $\textsf{GEM}(\alpha$)
# for $\alpha = \{ 0,1, 1, 10, 100 \}$ (Dilber, 2018).

# 1. Load necessary library
suppressMessages(suppressWarnings(library(ggplot2)))

# 2. Sample from GEM distribution

rGem <- function(alpha){
  p <- rbeta(n = 1, shape1 = 1, shape2 = alpha) # V_{1} ~ B(1, \alpha), p_{1} = V_{1}
  while(sum(p) < 1){ # p is a probability measure
    p <- c(p, (1 - sum(p))*rbeta(n = 1, shape1 = 1, shape2 = alpha)) # p_{k} = (1 - \sum_{i = 1}^{k - 1} p_{i}) V_{k}
  }
  return(p)
}

# 3. Normalize vector length

set_same_len <- function(vector, max_len){
  length(vector) <- max_len
  vector[is.na(vector)] <- 0 # Replace NAs with zero.
  
  return(vector)
}

# 4. Simulate multiple GEM distributions and compute the average of the samples

rGEM <- function(n_sim, alpha){
  p_list <- vector(mode = "list", n_sim)
  for (i in 1:n_sim) {
    p_list[[i]] <- rGem(alpha) # Simulate n_sim GEM distributions
  }
  max_len <- max(sapply(p_list, length)) # Maximum length vector
  simul <- sapply(p_list, function(x) set_same_len(x, max_len)) # Replace NAs with zero.
  return(rowMeans(simul)) # Compute the average of the samples
}

# 5. Set parameter alpha

alpha <- c(0.1, 1, 10, 100)

X <- lapply(alpha, function(x) cumsum(rGEM(n_sim = 10, x))) # GEM cumulative distribution function

label <- c()
prob <- c()

for (i in 1:length(X)) {
  label <- c(label, rep(alpha[i], length(X[[i]])))
  prob <- c(prob, X[[i]])
}

df <- data.frame(alpha = as.character(label), probability = prob)

# 6. Plot cumulative probabilities p

ggplot(df, aes(alpha, probability)) + 
  #theme(legend.position = "none") +
  geom_jitter(width = 0.025, aes(colour = alpha)) +
  xlab(expression(alpha)) +
  ylab("Probabilidad acumulada") +
  theme_classic()
 
