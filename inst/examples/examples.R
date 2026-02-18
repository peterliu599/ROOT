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
set.seed(123)

n  <- 200
X1 <- runif(n, -1, 1)
X2 <- runif(n, -1, 1)
# True optimal subgroup has XOR structure:
# top-right (X1>0, X2>0) and bottom-left (X1<0, X2<0)
true_w <- as.integer((X1 > 0 & X2 > 0) | (X1 < 0 & X2 < 0))
v      <- rnorm(n, mean = true_w, sd = 2)

dat_xor <- data.frame(v = v, X1 = X1, X2 = X2)

xor_fit <- ROOT(
  data      = dat_xor,
  num_trees = 20,
  top_k_trees = TRUE,
  k         = 10,
  seed      = 42
)

print(xor_fit)
plot(xor_fit)


### Another simple generalizability example
### Using the diabetes_data.rda dataset in the data folder
data(diabetes_data, package = "ROOT")

gen_fit <- characterizing_underrep(
  data                  = diabetes_data,
  generalizability_path = TRUE,
  num_trees             = 20,
  top_k_trees           = TRUE,
  k                     = 10,
  seed                  = 123
)

print(gen_fit)
plot(gen_fit)
