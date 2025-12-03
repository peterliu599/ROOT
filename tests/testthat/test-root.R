test_that("ROOT runs in two-sample mode and returns structured outputs", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("rpart")
  skip_if_not_installed("withr")
  skip_if_not_installed("gbm")   # used if feature_est = "GBM" in user code paths
  skip_if_not_installed("mlbench")

  sim <- get_data(n = 600, seed = 20)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  set.seed(42)
  out <- ROOT(
    data        = df,
    outcome     = "Yobs",
    treatment   = "Tr",
    sample      = "S",
    seed        = 99,
    num_trees   = 6,
    top_k_trees = TRUE,
    k           = 3,
    feature_est = "Ridge",
    verbose     = TRUE
  )

  expect_s3_class(out, "ROOT")
  expect_true(all(c("D_rash","D_forest","w_forest","rashomon_set","f","testing_data","estimate") %in% names(out)))
  expect_equal(length(out$w_forest), 6L)
  expect_true(length(out$rashomon_set) >= 1L)
  expect_true(all(out$D_rash$w_opt %in% 0:1))

  # Estimates present with labels and SDs (keys updated from se_* to sd_*)
  est <- out$estimate
  expect_true(all(c("estimand_unweighted","value_unweighted","se_unweighted",
                    "estimand_weighted","value_weighted","se_weighted",
                    "n_analysis","sum_w") %in% names(est)))
  expect_true(is.numeric(est$value_unweighted))
  expect_true(is.numeric(est$value_weighted) || is.na(est$value_weighted))

  # summary.ROOT prints without error
  #expect_visible(summary(out))
})

test_that("ROOT single-sample mode works when sample=NULL", {
  skip_if_not_installed("mlbench")
  sim <- get_data(n = 300, seed = 21)
  # Single-sample: keep only S==1 and drop sample
  dfS <- subset(sim$data, S == 1)
  dfS$Yobs <- dfS$Yobs

  out <- ROOT(
    data        = dfS[, c(grep("^X", names(dfS), value = TRUE), "Tr", "Yobs")],
    outcome     = "Yobs",
    treatment   = "Tr",
    sample      = NULL,
    num_trees   = 4,
    vote_threshold = 0.6
  )
  expect_s3_class(out, "ROOT")
  expect_true(out$single_sample_mode)
  expect_equal(out$estimate$estimand_unweighted, "SATE")
  expect_equal(out$estimate$estimand_weighted,   "WATE")
})

test_that("ROOT input validation catches errors", {
  sim <- get_data(n = 100, seed = 5)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  expect_error(ROOT(df, "Yobs", "Tr", sample = "S_missing"), "not found")
  expect_error(ROOT(df, "Yobs", "Tr", sample = "S", num_trees = 0), "must be positive")
  expect_error(ROOT(df, "Yobs", "Tr", sample = "S", vote_threshold = 1.1), "in \\(0, 1\\]")
})

mk_data_two_sample <- function(n=200, p=5, seed=11) {
  set.seed(seed)
  X <- as.data.frame(matrix(rnorm(n*p), n, p)); names(X) <- paste0("X", seq_len(p))
  Tr <- rbinom(n, 1, 0.5)
  S  <- rbinom(n, 1, 0.5)
  Y  <- rnorm(n)
  data.frame(X, Tr=Tr, S=S, Yobs=Y)
}

test_that("ROOT() validates arguments and 0/1 coercion", {
  d <- mk_data_two_sample()

  expect_error(ROOT(as.matrix(d), "Yobs", "Tr", "S"), "`data` must be a data frame")
  expect_error(ROOT(d, "Nope", "Tr", "S"), "`outcome` must be a single column name")
  expect_error(ROOT(d, "Yobs", "ZZ", "S"), "`treatment` must be a single column name")
  expect_error(ROOT(d, "Yobs", "Tr", "SS"), "`sample` column not found")

  expect_error(ROOT(d, "Yobs", "Tr", "S", leaf_proba = -0.1), "between 0 and 1")
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees = 0), "must be positive")
  expect_error(ROOT(d, "Yobs", "Tr", "S", vote_threshold = 0), "in \\(0, 1\\]")
  expect_error(ROOT(d, "Yobs", "Tr", "S", explore_proba = 2), "between 0 and 1")
  expect_error(ROOT(d, "Yobs", "Tr", "S", feature_est = 123), "must be \"Ridge\", \"GBM\", or a function")
  expect_error(ROOT(d, "Yobs", "Tr", "S", top_k_trees = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_error(ROOT(d, "Yobs", "Tr", "S", k = 0), "positive integer")
  expect_error(ROOT(d, "Yobs", "Tr", "S", cutoff = "nope"), "must be \"baseline\" or numeric")
  expect_error(ROOT(d, "Yobs", "Tr", "S", verbose = 1), "TRUE or FALSE")
  expect_error(ROOT(d, "Yobs", "Tr", "S", global_objective_fn = 1), "must be a function")

  # OK: treatment/sample coercion of strings
  d2 <- d
  d2$Tr <- sample(c("treated","control","yes","no","True","False"), nrow(d2), TRUE)
  d2$S  <- sample(c("1","0","yes","no","true","false"), nrow(d2), TRUE)
  # Won't error, but will stop later if estimate_dml needs numeric; just ensure coercion branch runs:
  expect_error(
    ROOT(d2, "Yobs", "Tr", "S", num_trees=1),
    regexp = NA
  )
})

test_that("ROOT single-sample mode via sample=NULL and via constant S", {
  d <- mk_data_two_sample()
  dS1 <- d; dS1$S <- 1L

  skip_if_not_installed("mlbench")

  # (a) sample=NULL path (SATE/WATE)
  set.seed(7)
  r1 <- ROOT(dS1[, c(grep("^X", names(dS1), value=TRUE), "Tr", "Yobs")],
             outcome="Yobs", treatment="Tr", sample=NULL, num_trees=2)
  expect_true(r1$single_sample_mode)
  expect_identical(r1$estimate$estimand_unweighted, "SATE")
  expect_identical(r1$estimate$estimand_weighted, "WATE")

  # (b) sample column present but constant -> also single_sample_mode TRUE
  set.seed(8)
  r2 <- ROOT(dS1, outcome="Yobs", treatment="Tr", sample="S", num_trees=2)
  expect_true(r2$single_sample_mode)
})

test_that("ROOT feature_est = Ridge / GBM / custom; bad custom errors", {
  d <- mk_data_two_sample(n=220, p=6)

  skip_if_not_installed("MASS")
  skip_if_not_installed("gbm")
  skip_if_not_installed("mlbench")

  # Ridge
  set.seed(1)
  rR <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, feature_est="Ridge")
  expect_s3_class(rR, "ROOT")

  # GBM
  set.seed(2)
  rG <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, feature_est="GBM")
  expect_s3_class(rG, "ROOT")

  # custom ok: returns named nonnegative importances for all X columns
  ok_imp <- function(X, y, ...) { setNames(abs(colSums(X^2)) + 1e-6, colnames(X)) }
  set.seed(3)
  rC <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, feature_est=ok_imp)
  expect_s3_class(rC, "ROOT")

  # custom bad: unnamed -> error from .norm_feat_prob
  bad1 <- function(X, y, ...) { as.numeric(abs(colSums(X))) }
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees=1, feature_est=bad1),
               "named numeric vector")

  # custom bad: missing a name -> error
  bad2 <- function(X, y, ...) {
    nm <- colnames(X); nm[length(nm)] <- "WRONG"
    setNames(abs(colSums(X)), nm)
  }
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees=1, feature_est=bad2),
               "Importance missing for some X_df columns")

  # custom bad: negative -> error
  bad3 <- function(X, y, ...) { setNames(c(-1, rep(1, ncol(X)-1)), colnames(X)) }
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees=1, feature_est=bad3),
               "must be non-negative")
})

test_that("Rashomon selection: top-k, baseline cutoff, and empty set warning", {
  skip_if_not_installed("mlbench")
  d <- mk_data_two_sample(n=260, p=6)

  set.seed(11)
  r_topk <- ROOT(d, "Yobs", "Tr", "S", num_trees=5, top_k_trees=TRUE, k=3)
  expect_length(r_topk$rashomon_set, 3)

  set.seed(12)
  r_base <- ROOT(d, "Yobs", "Tr", "S", num_trees=5, top_k_trees=FALSE, cutoff="baseline")
  expect_true(length(r_base$rashomon_set) >= 1)

  set.seed(13)
  # absurdly small cutoff forces empty set → warning; w_opt should exist but be all 0
  expect_warning({
    r_empty <- ROOT(d, "Yobs", "Tr", "S", num_trees=3, top_k_trees=FALSE, cutoff = -1e9)
    expect_length(r_empty$rashomon_set, 0)
    expect_true("w_opt" %in% names(r_empty$D_rash))
    expect_true(all(r_empty$D_rash$w_opt == 0))
  }, "No trees selected")
})

test_that("ROOT estimate fields populated; binary vs non-binary weighted SE behavior", {
  skip_if_not_installed("mlbench")
  d <- mk_data_two_sample(n=280, p=6)

  set.seed(21)
  r <- ROOT(d, "Yobs", "Tr", "S", num_trees=4)
  e <- r$estimate
  expect_true(all(c("estimand_unweighted","value_unweighted","se_unweighted",
                    "estimand_weighted","value_weighted","se_weighted",
                    "se_weighted_note","n_analysis","sum_w") %in% names(e)))
  # Binary case: either valid SE or NA when n_A=1
  expect_true(is.numeric(e$se_weighted) || is.na(e$se_weighted))

  # Force non-binary to hit omit-SE path in summary/print
  r2 <- r
  r2$D_rash$w_opt <- as.numeric(r2$D_rash$w_opt) + runif(nrow(r2$D_rash), 0, 0.1)
  # Call summary/print to ensure they don't error
  expect_silent(capture.output(summary(r2)))
  expect_silent(capture.output(print(r2)))
})

test_that("ROOT seed yields reproducible w_forest & estimates; verbose prints", {
  skip_if_not_installed("mlbench")
  d <- mk_data_two_sample(n=200, p=5)

  r1 <- ROOT(d, "Yobs", "Tr", "S", num_trees=3, seed=123, verbose=TRUE)
  r2 <- ROOT(d, "Yobs", "Tr", "S", num_trees=3, seed=123, verbose=TRUE)

  expect_equal(lapply(r1$w_forest, `[[`, "local objective"),
               lapply(r2$w_forest, `[[`, "local objective"))
  expect_equal(r1$estimate$value_unweighted, r2$estimate$value_unweighted, tolerance=1e-12)
  expect_equal(r1$estimate$value_weighted,  r2$estimate$value_weighted,  tolerance=1e-12)

  # verbose -> messages (do not assert exact strings, just that something prints)
  expect_output({
    r3 <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, seed=124, verbose=TRUE)
    print(r3)  # ensure print path runs too
  })
})
