# Tests for the extended vote_threshold parameter in ROOT() and
# characterizing_underrep().
#
# Coverage:
#   - numeric threshold (default and boundary values)
#   - string shortcuts: "majority", "mean", "median"
#   - custom aggregation function
#   - validation errors for every invalid input type
#   - custom function runtime errors (bad output shape, non-binary, internal error)
#   - vote_threshold stored on the returned ROOT object
#   - forwarded correctly through characterizing_underrep()

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

mk_vt_data <- function(n = 160, seed = 77) {
  set.seed(seed)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  S  <- rbinom(n, 1, plogis(0.3 * X1 - 0.2 * X2))
  Tr <- rbinom(n, 1, 0.5)
  Y  <- 1 + 0.8 * Tr + 0.5 * X1 + rnorm(n)
  data.frame(Y = Y, Tr = Tr, S = S, X1 = X1, X2 = X2)
}

# Suppress "all w_opt = 0/1" warnings and single-class messages that are
# unrelated to vote_threshold behaviour.
muffle_vt <- function(expr) {
  suppressMessages(
    withCallingHandlers(expr,
                        warning = function(w) {
                          if (grepl("w_opt = 0|w_opt = 1|No trees selected", conditionMessage(w)))
                            invokeRestart("muffleWarning")
                        }
    )
  )
}

# Run ROOT with small fixed settings so tests stay fast.
root_vt <- function(data, vote_threshold, seed = 1) {
  muffle_vt(ROOT(
    data                  = data,
    generalizability_path = TRUE,
    seed                  = seed,
    num_trees             = 6,
    top_k_trees           = TRUE,
    k                     = 4,
    vote_threshold        = vote_threshold
  ))
}

# ===========================================================================
# 1. Numeric threshold
# ===========================================================================

test_that("vote_threshold numeric 2/3: w_opt is binary and stored on object", {
  d   <- mk_vt_data()
  out <- root_vt(d, 2 / 3)

  expect_s3_class(out, "ROOT")
  expect_true(all(out$D_rash$w_opt %in% 0:1))
  expect_equal(out$vote_threshold, 2 / 3)
})

test_that("vote_threshold = 1: only unanimous rows get w_opt = 1", {
  d   <- mk_vt_data()
  out <- root_vt(d, 1.0)

  n_set         <- length(out$rashomon_set)
  unanimous     <- sum(out$D_rash$vote_count == n_set)
  expect_equal(sum(out$D_rash$w_opt == 1L), unanimous)
})

test_that("vote_threshold = 0.01: keeps more observations than 0.99 (permissive vs strict)", {
  d        <- mk_vt_data()
  out_low  <- root_vt(d, 0.01)
  out_high <- root_vt(d, 0.99)

  # A near-zero threshold is the most permissive numeric option:
  # any row with at least one vote=1 across the rashomon set passes it.
  # So it must keep at least as many observations as a near-one threshold.
  expect_true(mean(out_low$D_rash$w_opt) >= mean(out_high$D_rash$w_opt))
})

# ===========================================================================
# 2. String shortcuts
# ===========================================================================

test_that("vote_threshold 'majority': w_opt binary and string stored", {
  d   <- mk_vt_data()
  out <- root_vt(d, "majority")

  expect_s3_class(out, "ROOT")
  expect_true(all(out$D_rash$w_opt %in% 0:1))
  expect_equal(out$vote_threshold, "majority")
})

test_that("vote_threshold 'mean': w_opt binary and string stored", {
  d   <- mk_vt_data()
  out <- root_vt(d, "mean")

  expect_s3_class(out, "ROOT")
  expect_true(all(out$D_rash$w_opt %in% 0:1))
  expect_equal(out$vote_threshold, "mean")
})

test_that("vote_threshold 'median': w_opt binary and string stored", {
  d   <- mk_vt_data()
  out <- root_vt(d, "median")

  expect_s3_class(out, "ROOT")
  expect_true(all(out$D_rash$w_opt %in% 0:1))
  expect_equal(out$vote_threshold, "median")
})

test_that("'majority' and numeric 2/3 give identical w_opt (same seed)", {
  d       <- mk_vt_data()
  out_num <- root_vt(d, 2 / 3,      seed = 20)
  out_str <- root_vt(d, "majority",  seed = 20)

  expect_equal(out_num$D_rash$w_opt, out_str$D_rash$w_opt)
})

test_that("'mean' and numeric 0.5 give identical w_opt (same seed)", {
  d       <- mk_vt_data()
  out_num <- root_vt(d, 0.5,    seed = 21)
  out_str <- root_vt(d, "mean", seed = 21)

  expect_equal(out_num$D_rash$w_opt, out_str$D_rash$w_opt)
})

# ===========================================================================
# 3. Custom function
# ===========================================================================

test_that("vote_threshold custom fn: receives correct n x k binary matrix", {
  d        <- mk_vt_data()
  captured <- NULL

  out <- root_vt(d, function(votes) {
    captured <<- votes
    as.integer(rowMeans(votes) > 0.5)
  })

  expect_true(is.matrix(captured))
  expect_equal(nrow(captured), nrow(out$D_rash))
  expect_equal(ncol(captured), length(out$rashomon_set))
  expect_true(all(captured %in% 0:1))
})

test_that("vote_threshold custom fn always-keep: all w_opt = 1", {
  d   <- mk_vt_data()
  out <- muffle_vt(root_vt(d, function(v) rep(1L, nrow(v))))

  expect_true(all(out$D_rash$w_opt == 1L))
})

test_that("vote_threshold custom fn always-drop: all w_opt = 0", {
  d   <- mk_vt_data()
  out <- muffle_vt(root_vt(d, function(v) rep(0L, nrow(v))))

  expect_true(all(out$D_rash$w_opt == 0L))
})

test_that("vote_threshold custom fn is stored on returned object", {
  d  <- mk_vt_data()
  fn <- function(v) as.integer(rowMeans(v) > 0.6)
  out <- root_vt(d, fn)

  expect_true(is.function(out$vote_threshold))
})

# ===========================================================================
# 4. Validation errors
# ===========================================================================

test_that("vote_threshold = 0 is rejected", {
  d <- mk_vt_data()
  expect_error(
    ROOT(d, generalizability_path = TRUE, vote_threshold = 0),
    "single value in \\(0, 1\\]"
  )
})

test_that("vote_threshold = 1.1 is rejected", {
  d <- mk_vt_data()
  expect_error(
    ROOT(d, generalizability_path = TRUE, vote_threshold = 1.1),
    "single value in \\(0, 1\\]"
  )
})

test_that("vote_threshold as length-2 numeric vector is rejected", {
  d <- mk_vt_data()
  expect_error(
    ROOT(d, generalizability_path = TRUE, vote_threshold = c(0.5, 0.6)),
    "single value in \\(0, 1\\]"
  )
})

test_that("vote_threshold as unrecognised string is rejected", {
  d <- mk_vt_data()
  expect_error(
    ROOT(d, generalizability_path = TRUE, vote_threshold = "mode"),
    '"majority".*"mean".*"median"'
  )
})

test_that("vote_threshold as empty string is rejected", {
  d <- mk_vt_data()
  expect_error(
    ROOT(d, generalizability_path = TRUE, vote_threshold = ""),
    '"majority".*"mean".*"median"'
  )
})

test_that("vote_threshold as list is rejected", {
  d <- mk_vt_data()
  expect_error(
    ROOT(d, generalizability_path = TRUE, vote_threshold = list(0.5)),
    "numeric threshold"
  )
})

test_that("vote_threshold as logical is rejected", {
  d <- mk_vt_data()
  expect_error(
    ROOT(d, generalizability_path = TRUE, vote_threshold = TRUE),
    "numeric threshold"
  )
})

# ===========================================================================
# 5. Custom function runtime errors
# ===========================================================================

test_that("vote_threshold custom fn returning wrong length errors clearly", {
  d <- mk_vt_data()
  expect_error(
    root_vt(d, function(v) rep(1L, nrow(v) + 5L)),
    "integer vector of 0s and 1s"
  )
})

test_that("vote_threshold custom fn returning non-binary values errors clearly", {
  d <- mk_vt_data()
  expect_error(
    root_vt(d, function(v) rep(0.5, nrow(v))),
    "integer vector of 0s and 1s"
  )
})

test_that("vote_threshold custom fn that throws internally is caught cleanly", {
  d <- mk_vt_data()
  expect_error(
    root_vt(d, function(v) stop("boom inside custom fn")),
    "vote_threshold\\(\\) raised an error"
  )
})

# ===========================================================================
# 6. Forwarding through characterizing_underrep()
# ===========================================================================

test_that("characterizing_underrep forwards numeric vote_threshold", {
  d   <- mk_vt_data()
  out <- muffle_vt(characterizing_underrep(
    d, generalizability_path = TRUE,
    seed = 50, num_trees = 4, top_k_trees = TRUE, k = 3,
    vote_threshold = 0.5
  ))

  expect_s3_class(out, "characterizing_underrep")
  expect_equal(out$root$vote_threshold, 0.5)
  expect_true(all(out$root$D_rash$w_opt %in% 0:1))
})

test_that("characterizing_underrep forwards string vote_threshold", {
  d   <- mk_vt_data()
  out <- muffle_vt(characterizing_underrep(
    d, generalizability_path = TRUE,
    seed = 51, num_trees = 4, top_k_trees = TRUE, k = 3,
    vote_threshold = "median"
  ))

  expect_s3_class(out, "characterizing_underrep")
  expect_equal(out$root$vote_threshold, "median")
})

test_that("characterizing_underrep forwards custom function vote_threshold", {
  d  <- mk_vt_data()
  fn <- function(v) as.integer(rowMeans(v) > 0.5)
  out <- muffle_vt(characterizing_underrep(
    d, generalizability_path = TRUE,
    seed = 52, num_trees = 4, top_k_trees = TRUE, k = 3,
    vote_threshold = fn
  ))

  expect_s3_class(out, "characterizing_underrep")
  expect_true(is.function(out$root$vote_threshold))
  expect_true(all(out$root$D_rash$w_opt %in% 0:1))
})
