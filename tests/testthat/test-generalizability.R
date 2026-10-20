# Tests for the generalizability-path corrections:
#  (1) IPW scores carry the pi/(1-pi) prefactor, so mean(v) estimates the TATE;
#  (2) the default generalizability objective is the centered, kappa-scaled SE;
#  (3) kappa_w re-normalizes to pi_w = P(S = 1 | w = 1);
#  (5) general mode with only `v` centers at the weighted mean.

sim_gen <- function(seed = 1, n = 10000) {
  set.seed(seed)
  x1 <- rbinom(n, 1, 0.3)                       # under-sampled, effect modifier
  x2 <- rbinom(n, 1, 0.5)                       # not an effect modifier
  S  <- rbinom(n, 1, plogis(-1.2 - 1.0 * x1 + 0.3 * x2))
  Tr <- ifelse(S == 1, rbinom(n, 1, 0.5), NA)
  tau <- 1 + 2 * x1
  Y  <- ifelse(S == 1, 0.5 + x2 + Tr * tau + rnorm(n), NA)
  list(df = data.frame(x1 = x1, x2 = x2, S = S, Tr = Tr, Y = Y),
       tate = mean(tau[S == 0]))
}

test_that("mean of transport scores over trial units estimates the TATE", {
  cts <- getFromNamespace("compute_transport_scores", "ROOT")
  est <- vapply(1:20, function(s) {
    d <- sim_gen(s)
    sc <- suppressWarnings(cts(d$df, "Y", "Tr", "S"))
    c(mean(sc$v[d$df$S == 1]), d$tate)
  }, numeric(2))
  expect_lt(abs(mean(est[1, ]) - mean(est[2, ])), 0.1)
})

test_that("objective_generalizability is zero for a constant-score kept set", {
  obj <- getFromNamespace("objective_generalizability", "ROOT")
  D <- data.frame(v = c(2, 2, 2, 9), sel_odds = c(1, 2, 3, 4), w = c(1, 1, 1, 0))
  expect_equal(obj(D), 0)
})

test_that("objective_generalizability has kappa = 1 when w is all ones", {
  obj <- getFromNamespace("objective_generalizability", "ROOT")
  v <- c(1, 3, 5, 7); D <- data.frame(v = v, sel_odds = c(1, 2, 3, 4), w = 1)
  expect_equal(obj(D), sqrt(sum((v - mean(v))^2)) / length(v))
})

test_that("wtate_kappa equals ratio of refined to full selection odds", {
  kap <- getFromNamespace("wtate_kappa", "ROOT")
  r <- c(1, 1, 3, 3); A <- c(TRUE, TRUE, FALSE, FALSE)
  # refined odds n_A / sum_A r = 2/2 = 1; full odds 4/8 = 0.5
  expect_equal(kap(r, A), 2)
})

test_that("objective_default centers v at the weighted mean when vsq is absent", {
  obj <- getFromNamespace("objective_default", "ROOT")
  D <- data.frame(v = c(5, 5, 5, -100), w = c(1, 1, 1, 0))
  expect_equal(obj(D), 0)
})

test_that("ROOT() general mode with only `v` does not freeze centering", {
  set.seed(3)
  # Two covariates: MASS::lm.ridge() cannot fit a single-predictor model,
  # which would trigger ROOT's uniform-sampling fallback warning.
  df <- data.frame(x = rbinom(200, 1, 0.5), x2 = rbinom(200, 1, 0.5))
  df$v <- ifelse(df$x == 1, 10 + rnorm(200, 0, 0.01), rnorm(200, 0, 5))
  fit <- ROOT(data = df, num_trees = 5, top_k_trees = TRUE, k = 3, seed = 1)
  expect_false("vsq" %in% names(fit$D_forest))
})

test_that("ROOT() generalizability estimate is on the treatment-effect scale", {
  d <- sim_gen(7)
  fit <- suppressWarnings(ROOT(data = d$df, generalizability_path = TRUE,
                               num_trees = 5, top_k_trees = TRUE, k = 3, seed = 1))
  expect_lt(abs(fit$estimate$value_unweighted - d$tate), 0.5)
  expect_true(is.numeric(fit$estimate$kappa_w))
})
