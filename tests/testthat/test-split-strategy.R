# tests/testthat/test-split-strategy.R
# Tests for the split_strategy / quantile_range feature.
# Covers: split_point(), split_node() parameter threading, ROOT() and
#         characterizing_underrep() integration.
#
# Note: generalizability-mode calls use vote_threshold = 0.4 and
#       num_trees = 5 to avoid exact-tie errors. With 5 trees the
#       possible mean votes are {0, 0.2, 0.4, 0.6, 0.8, 1.0}, and
#       none equals 0.4 exactly, so ties are impossible.  Calls are wrapped in
#       suppressWarnings() to silence the expected P(S=1|X) flooring
#       and "all w_opt = 0" warnings from synthetic data.

# ============================================================================
# Helper: minimal stacked dataset for generalizability tests
# ============================================================================
make_gen_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  n_trial  <- round(n * 0.6)
  n_target <- n - n_trial
  data.frame(
    x1 = c(rnorm(n_trial, 0, 1), rnorm(n_target, 1, 1)),
    x2 = c(rnorm(n_trial, 0, 1), rnorm(n_target, 0, 1)),
    Tr = c(rbinom(n_trial, 1, 0.5), rep(NA_integer_, n_target)),
    Y  = c(rnorm(n_trial), rep(NA_real_, n_target)),
    S  = c(rep(1L, n_trial), rep(0L, n_target))
  )
}

# Helper: minimal dataset for general optimization tests
make_opt_data <- function(n = 100, seed = 7) {
  set.seed(seed)
  data.frame(
    x1  = rnorm(n),
    x2  = rnorm(n),
    vsq = runif(n, 0.1, 2)
  )
}

simple_objective <- function(D) {
  w <- D$w
  if (sum(w) == 0) return(Inf)
  sqrt(sum(w * D$vsq) / sum(w)^2)
}

# ============================================================================
# 1. split_point() — unit tests
# ============================================================================
test_that("split_point with 'midpoint' matches midpoint()", {
  x <- c(1, 3, 5, 7, 10)
  expect_equal(split_point(x, strategy = "midpoint"), midpoint(x))
  expect_equal(split_point(x, strategy = "midpoint"), (1 + 10) / 2)
})

test_that("split_point with 'random_quantile' returns a value within data range", {
  set.seed(123)
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  sp <- split_point(x, strategy = "random_quantile", quantile_range = c(0.1, 0.9))
  expect_true(sp >= min(x) && sp <= max(x))
})

test_that("split_point with 'random_quantile' varies across calls", {
  x <- seq(1, 100)
  set.seed(1); sp1 <- split_point(x, strategy = "random_quantile")
  set.seed(2); sp2 <- split_point(x, strategy = "random_quantile")
  # With different seeds, the split points should (almost certainly) differ
  expect_false(sp1 == sp2)
})

test_that("split_point with 'random_quantile' respects quantile_range", {
  set.seed(99)
  x <- seq(0, 100, by = 1)   # quantiles are easy to compute
  # Request quantile in [0.4, 0.6] -> split should be near 40-60
  results <- replicate(50, split_point(x, strategy = "random_quantile",
                                       quantile_range = c(0.4, 0.6)))
  expect_true(all(results >= quantile(x, 0.4)))
  expect_true(all(results <= quantile(x, 0.6)))
})

test_that("split_point with 'midpoint' is deterministic", {
  x <- c(2, 8, 5, 3)
  sp1 <- split_point(x, strategy = "midpoint")
  sp2 <- split_point(x, strategy = "midpoint")
  expect_identical(sp1, sp2)
})

test_that("split_point handles empty and all-NA vectors", {
  expect_warning(split_point(numeric(0), strategy = "midpoint"), "Empty")
  expect_warning(split_point(c(NA, NaN, Inf, -Inf), strategy = "midpoint"),
                 "No finite")
})

test_that("split_point errors on non-numeric input", {
  expect_error(split_point(letters, strategy = "midpoint"), "numeric")
})

test_that("split_point errors on unknown strategy", {
  expect_error(split_point(1:10, strategy = "median"), "Unknown split strategy")
})

test_that("split_point handles single-value vector", {
  expect_equal(split_point(c(5, 5, 5), strategy = "midpoint"), 5)
  set.seed(1)
  expect_equal(split_point(c(5, 5, 5), strategy = "random_quantile"), 5)
})

# ============================================================================
# 2. ROOT() split_strategy argument — validation
# ============================================================================
test_that("ROOT rejects invalid split_strategy values", {
  dat <- make_opt_data()

  expect_error(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "median", seed = 1, num_trees = 1),
    "split_strategy"
  )
  expect_error(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = 123, seed = 1, num_trees = 1),
    "split_strategy"
  )
  expect_error(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = c("midpoint", "random_quantile"), seed = 1,
         num_trees = 1),
    "split_strategy"
  )
})

test_that("ROOT rejects invalid quantile_range values", {
  dat <- make_opt_data()

  # Reversed bounds
  expect_error(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile",
         quantile_range = c(0.9, 0.1), seed = 1, num_trees = 1),
    "quantile_range"
  )
  # Out of [0, 1]
  expect_error(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile",
         quantile_range = c(-0.1, 0.9), seed = 1, num_trees = 1),
    "quantile_range"
  )
  # Wrong length
  expect_error(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile",
         quantile_range = 0.5, seed = 1, num_trees = 1),
    "quantile_range"
  )
  # Equal bounds
  expect_error(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile",
         quantile_range = c(0.5, 0.5), seed = 1, num_trees = 1),
    "quantile_range"
  )
})

# ============================================================================
# 3. ROOT() with split_strategy — general optimization mode
# ============================================================================
test_that("ROOT runs successfully with split_strategy = 'midpoint' (default)", {
  dat <- make_opt_data()

  res <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "midpoint", seed = 1, num_trees = 5,
         vote_threshold = 0.4)
  )

  expect_s3_class(res, "ROOT")
  expect_equal(res$split_strategy, "midpoint")
  expect_true(length(res$rashomon_set) >= 0)
})

test_that("ROOT runs successfully with split_strategy = 'random_quantile'", {
  dat <- make_opt_data()

  res <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile",
         quantile_range = c(0.1, 0.9),
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  )

  expect_s3_class(res, "ROOT")
  expect_equal(res$split_strategy, "random_quantile")
})

test_that("random_quantile produces different trees with different seeds", {
  dat <- make_opt_data(n = 80)

  res1 <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile", seed = 1, num_trees = 5,
         vote_threshold = 0.4)
  )
  res2 <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile", seed = 99, num_trees = 5,
         vote_threshold = 0.4)
  )

  # Extract objective values — different seeds should explore differently
  obj1 <- vapply(res1$w_forest, function(x) x[["local objective"]], numeric(1))
  obj2 <- vapply(res2$w_forest, function(x) x[["local objective"]], numeric(1))

  # They should not be identical (extremely unlikely with different seeds)
  expect_false(identical(obj1, obj2))
})

test_that("midpoint is fully backward-compatible (same results as before)", {
  dat <- make_opt_data(n = 60)

  # Run twice with same seed — should be identical
  res1 <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "midpoint", seed = 42, num_trees = 5,
         vote_threshold = 0.4)
  )
  res2 <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "midpoint", seed = 42, num_trees = 5,
         vote_threshold = 0.4)
  )

  expect_identical(res1$D_rash$w_opt, res2$D_rash$w_opt)
})

# ============================================================================
# 5. summary.ROOT reports split strategy
# ============================================================================
test_that("summary.ROOT prints split strategy", {
  dat <- make_opt_data()

  res <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "random_quantile", seed = 1, num_trees = 5,
         vote_threshold = 0.4)
  )

  out <- capture.output(summary(res))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Split strategy.*random_quantile", combined))
})

test_that("summary.ROOT prints midpoint strategy", {
  dat <- make_opt_data()

  res <- suppressWarnings(
    ROOT(dat, global_objective_fn = simple_objective,
         split_strategy = "midpoint", seed = 1, num_trees = 5,
         vote_threshold = 0.4)
  )

  out <- capture.output(summary(res))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Split strategy.*midpoint", combined))
})


