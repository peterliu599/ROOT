library(devtools)
library(tidyverse)
#library(pbapply)

remove.packages("ROOT")
.rs.restartR()

devtools::document()
res <- devtools::check()
res$errors
res$warnings
res$notes

devtools::install()

library(ROOT)
devtools::build_manual()


### A simple generalizability example
### Using the diabetes_data.rda dataset in the data folder

ROOT.output <- ROOT(diabetes_data, generalizability_path = TRUE, seed=123)
summary(ROOT.output)
print(ROOT.output)
plot(ROOT.output)

char.output <- characterizing_underrep(diabetes_data,generalizability_path = TRUE, seed = 123)
summary(char.output)
print(char.output)
plot(char.output)

### A simple optimization example
### ROOT minimizes variance globally
set.seed(123)
n <- 1000

X1 <- sample(0:1, n, replace = TRUE)
X2 <- sample(0:1, n, replace = TRUE)

# XOR pattern: low variance on diagonal (X1==X2), high variance off-diagonal
Y <- ifelse(X1 == X2,
            rnorm(n, mean = 5, sd = 1),
            rnorm(n, mean = 5, sd = 10))

data <- data.frame(X1 = X1, X2 = X2, v = Y)

variance_objective <- function(D) {
  w <- D$w
  if (sum(w, na.rm = TRUE) < 2) return(Inf)

  y_kept <- D$v[w == 1]
  sd(y_kept)
}

root.result <- ROOT(
  data = data,
  global_objective_fn = variance_objective,
  generalizability_path = FALSE
)
summary(root.result)
print(root.result)
plot(root.result)

underrep.result <- characterizing_underrep(
  data = data,
  global_objective_fn = variance_objective,
  generalizability_path = FALSE
)

summary(underrep.result)
plot(underrep.result)


### A simple optimization example with default global objective function
library(ROOT)
set.seed(123)

n_assets <- 100

# Asset features
volatility  <- runif(n_assets, 0.05, 0.40)   # annualised vol
beta        <- runif(n_assets, 0.5,  1.8)     # market beta
sector      <- sample(c("Tech", "Finance", "Energy", "Health"), n_assets, replace = TRUE)

# Simulate returns correlated with beta and volatility
market      <- rnorm(1000, 0.0005, 0.01)
returns_mat <- sapply(seq_len(n_assets), function(i)
  beta[i] * market + rnorm(1000, 0, volatility[i] / sqrt(252))
)
vsq <- apply(returns_mat, 2, var)   # per-asset return variance

dat_portfolio <- data.frame(
  vsq      = vsq,
  vol      = volatility,
  beta     = beta,
  sector   = as.integer(factor(sector))   # encode as integer for splitting
)

portfolio_fit <- ROOT(
  data        = dat_portfolio,
  num_trees   = 20,
  top_k_trees = TRUE,
  k           = 10,
  seed        = 123
)

print(portfolio_fit)
plot(portfolio_fit)


### Another simple generalizability example
### Using the diabetes_data.rda dataset in the data folder
data(diabetes_data, package = "ROOT")

gen_fit <- characterizing_underrep(
  data                  = diabetes_data,
  generalizability_path = TRUE,
  num_trees             = 20,
  top_k_trees           = TRUE,
  k                     = 15,
  seed                  = 12
)

print(gen_fit)
plot(gen_fit)
