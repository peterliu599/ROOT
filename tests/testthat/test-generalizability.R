# Tests for the generalizability-path corrections:
#  (1) IPW scores carry the pi/(1-pi) prefactor, so mean(v) estimates the TATE;
#  (2) the default generalizability objective is the centered, kappa-scaled SE;
#  (3) kappa_w re-normalizes to pi_w = P(S = 1 | w = 1);
#  (5) general mode with only `v` centers at the weighted mean.

# Note: with the default vote_threshold = 2/3, a Rashomon set of k = 3 (or any
# multiple of 3) trees can produce a mean vote of exactly 2/3, which ROOT
# rejects as undefined. Tests that use the default threshold therefore use
# k = 5, for which mean votes are multiples of 0.2 and can never equal 2/3.

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
  fit <- ROOT(data = df, num_trees = 5, top_k_trees = TRUE, k = 5, seed = 1)
  expect_false("vsq" %in% names(fit$D_forest))
})

test_that("ROOT() generalizability estimate is on the treatment-effect scale", {
  d <- sim_gen(7)
  fit <- suppressWarnings(ROOT(data = d$df, generalizability_path = TRUE,
                               num_trees = 5, top_k_trees = TRUE, k = 5, seed = 1))
  expect_lt(abs(fit$estimate$value_unweighted - d$tate), 0.5)
  expect_true(is.numeric(fit$estimate$kappa_w))
})

# ---------------------------------------------------------------------------
# Additional coverage for new branches
# ---------------------------------------------------------------------------

test_that("objective_generalizability errors when required columns are missing", {
  obj <- getFromNamespace("objective_generalizability", "ROOT")
  expect_error(obj(data.frame(w = 1, v = 1)), "sel_odds")
})

test_that("objective_generalizability returns Inf with fewer than two kept units", {
  obj <- getFromNamespace("objective_generalizability", "ROOT")
  D <- data.frame(v = c(1, 2, 3), sel_odds = c(1, 1, 1), w = c(1, 0, 0))
  expect_identical(obj(D), Inf)
})

test_that("objective_generalizability returns Inf when the value is not finite", {
  obj <- getFromNamespace("objective_generalizability", "ROOT")
  # sel_odds of 0 in the kept set gives kappa_w = Inf and Inf * 0 = NaN
  D <- data.frame(v = c(2, 2, 5), sel_odds = c(0, 0, 1), w = c(1, 1, 0))
  expect_identical(obj(D), Inf)
})

test_that("generalizability estimate handles an empty kept set", {
  d <- sim_gen(11)
  fit <- suppressWarnings(suppressMessages(ROOT(
    data = d$df, generalizability_path = TRUE,
    num_trees = 4, top_k_trees = TRUE, k = 3, seed = 2,
    vote_threshold = function(votes) integer(nrow(votes))   # all w_opt = 0
  )))
  expect_true(is.na(fit$estimate$value_weighted))
  expect_true(is.na(fit$estimate$kappa_w))
  expect_identical(fit$estimate$n_A, 0L)
  expect_true(is.finite(fit$estimate$value_unweighted))
})

test_that("verbose output notes a user-supplied objective in generalizability mode", {
  d <- sim_gen(12)
  my_obj <- function(D) {
    w <- D$w
    if (sum(w) < 2) return(Inf)
    stats::sd(D$v[w == 1]) / sqrt(sum(w))
  }
  msgs <- capture_messages(suppressWarnings(ROOT(
    data = d$df, generalizability_path = TRUE, global_objective_fn = my_obj,
    num_trees = 4, top_k_trees = TRUE, k = 3, seed = 3,
    vote_threshold = function(votes) rep(1L, nrow(votes)),  # keep everyone
    verbose = TRUE
  )))
  expect_true(any(grepl("STATE \\(unweighted\\)", msgs)))
  expect_true(any(grepl("You supplied a custom global_objective_fn", msgs)))
})

test_that("general mode without vsq, v, or objective uses vsq = 1 with a message", {
  set.seed(5)
  df <- data.frame(x1 = rbinom(60, 1, 0.5), x2 = rnorm(60))
  expect_message(
    suppressWarnings(ROOT(data = df, num_trees = 5, top_k_trees = TRUE, k = 5, seed = 1)),
    "using vsq = 1"
  )
})

test_that("no vsq fallback message when a custom objective is supplied", {
  set.seed(6)
  df <- data.frame(x1 = rbinom(60, 1, 0.5), x2 = rnorm(60))
  obj <- function(D) if (sum(D$w) == 0) Inf else 1 / sum(D$w)
  msgs <- capture_messages(suppressWarnings(ROOT(
    data = df, global_objective_fn = obj,
    num_trees = 5, top_k_trees = TRUE, k = 5, seed = 1
  )))
  expect_false(any(grepl("using vsq = 1", msgs)))
})

test_that("covariance-based portfolio objective runs in general mode", {
  set.seed(123)
  n <- 30
  beta <- runif(n, 0.5, 1.8); vol <- runif(n, 0.05, 0.4)
  mkt <- rnorm(200, 0, 0.01)
  R <- sapply(seq_len(n), function(i) beta[i] * mkt + rnorm(200, 0, vol[i] / sqrt(252)))
  Sigma_hat <- stats::cov(R)
  pv <- function(D) {
    w <- D$w
    if (sum(w) == 0) return(Inf)
    sum(Sigma_hat[w == 1, w == 1]) / sum(w)^2
  }
  dat <- data.frame(vsq = diag(Sigma_hat), beta = beta, vol = vol)
  fit <- suppressWarnings(suppressMessages(ROOT(
    data = dat, global_objective_fn = pv,
    num_trees = 5, top_k_trees = TRUE, k = 5, seed = 4
  )))
  expect_s3_class(fit, "ROOT")
  expect_true(all(fit$D_rash$w_opt %in% 0:1))
  expect_false("vsq" %in% c(names(fit$f$variable.importance)))  # vsq never a split feature
})
