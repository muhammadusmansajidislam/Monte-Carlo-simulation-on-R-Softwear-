# Monte-Carlo-simulation-on-R-Softwear-
This Monte Carlo simulation study assesses the performance of OLS regression under multicollinearity and heteroscedasticity using 1,000 simulated datasets. Results show that OLS estimates remain approximately unbiased, but multicollinearity increases estimation variability and heteroscedasticity makes conventional standard errors unreliable, highlighting the need for robust inference methods.
# OLS under Multicollinearity and Heteroscedasticity

library(car)
library(lmtest)
library(sandwich)

set.seed(123)

# Settings
nsim <- 1000
n <- 100
beta <- c(2, 3, 1.5)

# Storage
B <- matrix(0, nsim, 3)

# Monte Carlo Simulation
for(i in 1:nsim){

  x1 <- rnorm(n)
  x2 <- 0.8*x1 + rnorm(n, 0, 0.2)      # Multicollinearity

  sigma <- 1 + 0.5*abs(x1)             # Heteroscedasticity
  e <- rnorm(n, 0, sigma)

  y <- beta[1] + beta[2]*x1 + beta[3]*x2 + e

  B[i, ] <- coef(lm(y ~ x1 + x2))
}

# Performance Measures
Mean <- colMeans(B)
Bias <- Mean - beta
MSE  <- colMeans((B - matrix(beta, nsim, 3, byrow = TRUE))^2)
RMSE <- sqrt(MSE)

cat("\nMean Estimates:\n"); print(Mean)
cat("\nBias:\n"); print(Bias)
cat("\nMSE:\n"); print(MSE)
cat("\nRMSE:\n"); print(RMSE)

# Diagnostic Sample
x1 <- rnorm(n)
x2 <- 0.8*x1 + rnorm(n,0,0.2)
sigma <- 1 + 0.5*abs(x1)
y <- beta[1] + beta[2]*x1 + beta[3]*x2 + rnorm(n,0,sigma)

model <- lm(y ~ x1 + x2)

# Diagnostics
summary(model)
vif(model)
dwtest(model)
bptest(model)

# Robust SE
coeftest(model, vcov = vcovHC(model, type = "HC3"))

# Diagnostic Plots
par(mfrow = c(2,2))
plot(model)
