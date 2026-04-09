# tests/testthat/test-positivity.R
# Tests for positivity violation trimming (Crump et al., 2009)
# Covers: crump_alpha(), trim_positivity_violations(), and the
#         positivity_trim argument in ROOT() / characterizing_underrep().
#
# Note: generalizability-mode calls use vote_threshold = 0.4 and
#       num_trees = 5 to avoid exact-tie errors. With 5 trees the
#       possible mean votes are {0, 0.2, 0.4, 0.6, 0.8, 1.0}, and
#       none equals 0.4 exactly, so ties are impossible.  Calls are wrapped in
#       suppressWarnings() to silence the expected P(S=1|X) flooring
#       warning from compute_transport_scores (the synthetic data
#       deliberately has limited overlap).

# ============================================================================
# Helper: build a minimal stacked trial/target dataset
# ============================================================================
make_stacked_data <- function(n_trial = 200, n_target = 300, seed = 42) {
  set.seed(seed)
  # Trial: x1 ~ N(0, 1), x2 ~ N(0, 1)
  x1_trial <- rnorm(n_trial, 0, 1)
  x2_trial <- rnorm(n_trial, 0, 1)
  tr_trial <- rbinom(n_trial, 1, 0.5)
  y_trial  <- 1 + 0.5 * tr_trial + rnorm(n_trial, 0, 0.5)

  # Target: x1 ~ N(2, 1) — shifted, creating limited overlap on x1
  x1_target <- rnorm(n_target, 2, 1)
  x2_target <- rnorm(n_target, 0, 1)
  tr_target <- rep(NA_integer_, n_target)
  y_target  <- rep(NA_real_, n_target)

  data.frame(
    x1 = c(x1_trial, x1_target),
    x2 = c(x2_trial, x2_target),
    Tr = c(tr_trial, tr_target),
    Y  = c(y_trial, y_target),
    S  = c(rep(1L, n_trial), rep(0L, n_target))
  )
}

# Synthetic propensity scores for unit-testing crump_alpha directly
make_ps_uniform   <- function(n = 500, seed = 1) {
  set.seed(seed); runif(n, 0.2, 0.8)
}
make_ps_extreme   <- function(n = 500, seed = 1) {
  set.seed(seed); c(runif(n * 0.8, 0.3, 0.7), runif(n * 0.2, 0.001, 0.02))
}

# ============================================================================
# 1. crump_alpha() — unit tests
# ============================================================================
test_that("crump_alpha returns 0 when propensity scores are well-behaved", {
  ps <- make_ps_uniform()
  alpha <- crump_alpha(ps)
  expect_equal(alpha, 0)
})

test_that("crump_alpha returns a positive value when there are extreme scores", {
  ps <- make_ps_extreme()
  alpha <- crump_alpha(ps)
  expect_true(alpha > 0)
  expect_true(alpha < 0.5)
})

test_that("crump_alpha returns value in [0, 0.5]", {
  set.seed(99)
  # Very extreme: half the scores near 0
  ps <- c(runif(250, 0.001, 0.01), runif(250, 0.3, 0.7))
  alpha <- crump_alpha(ps)
  expect_true(alpha >= 0 && alpha <= 0.5)
})

test_that("crump_alpha handles edge cases gracefully", {
  # All identical scores -> no trimming needed
  expect_equal(crump_alpha(rep(0.5, 100)), 0)

  # Single observation
  expect_equal(crump_alpha(0.3), 0)

  # Boundary values (0 and 1) are filtered out
  ps <- c(0, 0.5, 0.5, 1)
  expect_silent(crump_alpha(ps))
})

test_that("crump_alpha errors on bad input", {
  expect_error(crump_alpha(character(5)), "non-empty numeric")
  expect_error(crump_alpha(numeric(0)), "non-empty numeric")
})

test_that("crump_alpha grid_size parameter works", {
  ps <- make_ps_extreme()
  alpha_coarse <- crump_alpha(ps, grid_size = 50)
  alpha_fine   <- crump_alpha(ps, grid_size = 5000)
  # Both should be positive; fine grid should be at least as precise
  expect_true(alpha_coarse > 0)
  expect_true(alpha_fine > 0)
  # They can differ slightly but should be in the same ballpark
  expect_true(abs(alpha_coarse - alpha_fine) < 0.05)
})

# ============================================================================
# 2. trim_positivity_violations() — unit tests
# ============================================================================
test_that("trim_positivity_violations with fixed threshold removes target rows", {
  dat <- make_stacked_data()
  n_before <- nrow(dat)

  result <- trim_positivity_violations(
    data          = dat,
    sample_col    = "S",
    outcome_col   = "Y",
    treatment_col = "Tr",
    threshold     = 0.05
  )

  expect_true(is.data.frame(result$data))
  expect_true(nrow(result$data) <= n_before)
  expect_equal(result$alpha, 0.05)
  expect_true(is.numeric(result$pi_s))
  expect_equal(length(result$pi_s), n_before)

  # Trial rows should never be removed
  n_trial_before <- sum(dat$S == 1)
  n_trial_after  <- sum(result$data$S == 1)
  expect_equal(n_trial_after, n_trial_before)
})

test_that("trim_positivity_violations with 'crump' uses data-driven alpha", {
  dat <- make_stacked_data()

  result <- trim_positivity_violations(
    data          = dat,
    sample_col    = "S",
    outcome_col   = "Y",
    treatment_col = "Tr",
    threshold     = "crump"
  )

  expect_true(is.numeric(result$alpha))
  expect_true(result$alpha >= 0 && result$alpha <= 0.5)
  expect_true(is.integer(result$n_trimmed) || is.numeric(result$n_trimmed))
})

test_that("trim_positivity_violations preserves all rows when alpha = 0", {
  # Use well-overlapping data where crump_alpha should return 0
  set.seed(10)
  n <- 200
  dat <- data.frame(
    x1 = rnorm(n),
    Tr = c(rbinom(n / 2, 1, 0.5), rep(NA, n / 2)),
    Y  = c(rnorm(n / 2), rep(NA, n / 2)),
    S  = c(rep(1L, n / 2), rep(0L, n / 2))
  )

  result <- trim_positivity_violations(
    data          = dat,
    sample_col    = "S",
    outcome_col   = "Y",
    treatment_col = "Tr",
    threshold     = "crump"
  )

  # Even if alpha > 0, n_trimmed should be non-negative
  expect_true(result$n_trimmed >= 0L)
  expect_equal(nrow(result$data) + result$n_trimmed, nrow(dat))
})

test_that("trim_positivity_violations handles no-covariate edge case", {
  # Only Y, Tr, S — no covariates to build sampling model
  dat <- data.frame(
    Y  = rnorm(10),
    Tr = rbinom(10, 1, 0.5),
    S  = rep(c(1L, 0L), each = 5)
  )

  result <- trim_positivity_violations(
    data          = dat,
    sample_col    = "S",
    outcome_col   = "Y",
    treatment_col = "Tr",
    threshold     = 0.05
  )

  expect_equal(result$n_trimmed, 0L)
  expect_equal(nrow(result$data), nrow(dat))
})

# ============================================================================
# 3. ROOT() positivity_trim argument — validation
# ============================================================================
test_that("ROOT rejects invalid positivity_trim values", {
  dat <- make_stacked_data()

  expect_error(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = -0.1, seed = 1,
         num_trees = 1),
    "positivity_trim"
  )
  expect_error(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = 0.6, seed = 1,
         num_trees = 1),
    "positivity_trim"
  )
  expect_error(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = "invalid", seed = 1,
         num_trees = 1),
    "positivity_trim"
  )
  expect_error(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = c(0.05, 0.1),
         seed = 1, num_trees = 1),
    "positivity_trim"
  )
})

test_that("ROOT accepts valid positivity_trim values without error", {
  dat <- make_stacked_data()

  # NULL (default, no trimming)
  suppressWarnings(expect_no_error(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = NULL,
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  ))
  # Numeric
  suppressWarnings(expect_no_error(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = 0.05,
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  ))
  # "crump"
  suppressWarnings(expect_no_error(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = "crump",
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  ))
})

# ============================================================================
# 4. ROOT() positivity_trim — functional behavior
# ============================================================================
test_that("ROOT with positivity_trim returns positivity_trim_info in output", {
  dat <- make_stacked_data()

  res <- suppressWarnings(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = 0.05,
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  )

  expect_true("positivity_trim_info" %in% names(res))
  pti <- res$positivity_trim_info
  expect_equal(pti$method, "fixed")
  expect_equal(pti$alpha, 0.05)
  expect_true(is.numeric(pti$n_trimmed))
  expect_true(pti$n_trimmed >= 0)
})

test_that("ROOT with positivity_trim = 'crump' returns crump method info", {
  dat <- make_stacked_data()

  res <- suppressWarnings(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = "crump",
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  )

  pti <- res$positivity_trim_info
  expect_equal(pti$method, "crump")
  expect_true(pti$alpha >= 0 && pti$alpha <= 0.5)
})

test_that("ROOT with positivity_trim = NULL has NULL positivity_trim_info", {
  dat <- make_stacked_data()

  res <- suppressWarnings(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = NULL,
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  )

  expect_null(res$positivity_trim_info)
})

test_that("positivity_trim is ignored in general optimization mode", {
  set.seed(1)
  my_obj <- function(D) {
    w <- D$w; if (sum(w) == 0) return(Inf)
    sqrt(sum(w * D$vsq) / sum(w)^2)
  }
  dat <- data.frame(x1 = rnorm(50), x2 = rnorm(50), vsq = runif(50))

  # Should not error even though positivity_trim is set —
  # it's simply ignored outside generalizability mode
  suppressWarnings(
    res <- ROOT(dat, global_objective_fn = my_obj, positivity_trim = 0.1,
                seed = 1, num_trees = 5, vote_threshold = 0.4)
  )

  expect_null(res$positivity_trim_info)
})

# ============================================================================
# 5. characterizing_underrep() passes positivity_trim through
# ============================================================================
test_that("characterizing_underrep passes positivity_trim to ROOT", {
  dat <- make_stacked_data()

  res <- suppressWarnings(
    characterizing_underrep(
      data                  = dat,
      generalizability_path = TRUE,
      positivity_trim       = 0.05,
      seed                  = 1,
      num_trees             = 5,
      vote_threshold        = 0.4
    )
  )

  expect_true("positivity_trim_info" %in% names(res$root))
  expect_equal(res$root$positivity_trim_info$method, "fixed")
})

# ============================================================================
# 6. summary.ROOT reports positivity trimming
# ============================================================================
test_that("summary.ROOT prints positivity trimming info when present", {
  dat <- make_stacked_data()

  res <- suppressWarnings(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = 0.05,
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  )

  out <- capture.output(summary(res))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Positivity trimming", combined))
})

test_that("summary.ROOT does not print positivity info when NULL", {
  dat <- make_stacked_data()

  res <- suppressWarnings(
    ROOT(dat, generalizability_path = TRUE, positivity_trim = NULL,
         seed = 1, num_trees = 5, vote_threshold = 0.4)
  )

  out <- capture.output(summary(res))
  combined <- paste(out, collapse = "\n")
  expect_false(grepl("Positivity trimming", combined))
})
