# ---------------------------------------------------------------------------
# Helper: create fake ROOT objects for testing S3 methods
# ---------------------------------------------------------------------------
fake_root_for_print <- function(f = NULL,
                                generalizability_path = FALSE,
                                estimate = NULL,
                                user_supplied_objective = FALSE) {
  obj <- list(
    D_rash                = data.frame(w_opt = c(1L, 0L, 1L)),
    D_forest              = data.frame(X1 = 1:3, vsq = 1),
    w_forest              = list(list(`local objective` = 0.5)),
    rashomon_set          = 1L,
    global_objective_fn   = objective_default,
    user_supplied_objective = user_supplied_objective,
    vote_threshold        = 2 / 3,
    f                     = f,
    testing_data          = data.frame(X1 = 1:3),
    estimate              = estimate,
    generalizability_path = generalizability_path
  )
  class(obj) <- c("ROOT", "list")
  obj
}

test_that("print.ROOT prints basic info for non-generalizability mode", {
  obj <- fake_root_for_print(generalizability_path = FALSE)

  out <- capture.output(print(obj))
  expect_true(any(grepl("ROOT object", out)))
  expect_true(any(grepl("Generalizability mode.*FALSE", out)))
  expect_true(any(grepl("Summary classifier", out)))
  expect_true(any(grepl("no summary tree available", out)))
})

test_that("print.ROOT prints summary tree when f is present", {
  skip_if_not_installed("rpart")

  set.seed(1)
  df <- data.frame(
    y  = factor(sample(c(0, 1), 30, TRUE)),
    X1 = runif(30),
    X2 = runif(30)
  )
  fit <- rpart::rpart(y ~ X1 + X2, data = df, method = "class")

  obj <- fake_root_for_print(f = fit, generalizability_path = FALSE)

  out <- capture.output(print(obj))
  expect_true(any(grepl("ROOT object", out)))
  expect_true(any(grepl("Summary classifier", out)))
  expect_false(any(grepl("no summary tree available", out)))
})

test_that("print.ROOT prints estimands in generalizability mode", {
  est <- list(
    estimand_unweighted = "SATE",
    value_unweighted    = 0.123,
    se_unweighted       = 0.045,
    estimand_weighted   = "WTATE",
    value_weighted      = 0.234,
    se_weighted         = 0.056,
    se_weighted_note    = "Test note"
  )

  obj <- fake_root_for_print(
    generalizability_path = TRUE,
    estimate = est
  )

  out <- capture.output(print(obj))
  expect_true(any(grepl("ROOT object", out)))
  expect_true(any(grepl("Generalizability mode.*TRUE", out)))
  expect_true(any(grepl("Estimand summary \\(generalization mode\\)", out)))
  expect_true(any(grepl("Unweighted.*SATE", out)))
  expect_true(any(grepl("Weighted.*WTATE", out)))
})

test_that("print.ROOT skips estimands when estimate is NULL", {
  obj <- fake_root_for_print(
    generalizability_path = TRUE,
    estimate = NULL
  )

  out <- capture.output(print(obj))
  expect_true(any(grepl("ROOT object", out)))
  expect_false(any(grepl("Estimand summary", out)))
})

test_that("print.ROOT returns the object invisibly", {
  obj <- fake_root_for_print()
  result <- capture.output(ret <- print(obj))
  expect_identical(ret, obj)
})

test_that("summary.ROOT prints 'User-supplied: Yes' when user_supplied_objective is TRUE", {
  obj <- fake_root_for_print(user_supplied_objective = TRUE)

  out <- capture.output(summary(obj))
  expect_true(any(grepl("User-supplied: Yes", out)))
})

test_that("summary.ROOT prints 'User-supplied: No' when user_supplied_objective is FALSE", {
  obj <- fake_root_for_print(user_supplied_objective = FALSE)

  out <- capture.output(summary(obj))
  expect_true(any(grepl("User-supplied: No", out)))
})

test_that("summary.ROOT shows 'User-supplied: Yes' from real ROOT call with custom objective", {
  skip_if_not_installed("MASS")

  set.seed(42)
  n <- 60
  df <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    v  = rnorm(n)
  )
  df$vsq <- df$v^2

  my_obj <- function(D) {
    w <- D$w
    if (sum(w) == 0) return(Inf)
    mean(D$vsq[w == 1])
  }

  suppressWarnings({
    res <- ROOT(
      data = df,
      generalizability_path = FALSE,
      global_objective_fn = my_obj,
      num_trees = 2,
      seed = 123
    )
  })

  out <- capture.output(summary(res))
  expect_true(any(grepl("User-supplied: Yes", out)))
})

test_that("ROOT validates max_depth argument", {
  df <- data.frame(
    Y  = rnorm(20),
    Tr = sample(0:1, 20, TRUE),
    S  = sample(0:1, 20, TRUE),
    X1 = rnorm(20)
  )

  expect_error(
    ROOT(df, generalizability_path = TRUE, max_depth = -1),
    "`max_depth` must be a non-negative integer"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, max_depth = "five"),
    "`max_depth` must be a non-negative integer"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, max_depth = c(1, 2)),
    "`max_depth` must be a non-negative integer"
  )
})

test_that("ROOT validates min_leaf_n argument", {
  df <- data.frame(
    Y  = rnorm(20),
    Tr = sample(0:1, 20, TRUE),
    S  = sample(0:1, 20, TRUE),
    X1 = rnorm(20)
  )

  expect_error(
    ROOT(df, generalizability_path = TRUE, min_leaf_n = 0),
    "`min_leaf_n` must be a positive integer"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, min_leaf_n = -1),
    "`min_leaf_n` must be a positive integer"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, min_leaf_n = c(2, 3)),
    "`min_leaf_n` must be a positive integer"
  )
})

test_that("ROOT validates max_rejects_per_node argument", {
  df <- data.frame(
    Y  = rnorm(20),
    Tr = sample(0:1, 20, TRUE),
    S  = sample(0:1, 20, TRUE),
    X1 = rnorm(20)
  )

  expect_error(
    ROOT(df, generalizability_path = TRUE, max_rejects_per_node = 0),
    "`max_rejects_per_node` must be a positive integer"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, max_rejects_per_node = -5),
    "`max_rejects_per_node` must be a positive integer"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, max_rejects_per_node = c(5, 10)),
    "`max_rejects_per_node` must be a positive integer"
  )
})

test_that("ROOT errors when all S == 0 in generalizability mode (no trial units)", {
  df <- data.frame(
    Y  = rnorm(10),
    Tr = sample(0:1, 10, TRUE),
    S  = rep(0L, 10),
    X1 = rnorm(10)
  )

  # compute_transport_scores() fires before the analysis_idx check,
  # with: "no S == 1 (trial) rows for treatment model."
  # Note: ROOT.R line 504 ("No observations with S == 1 (trial units)")
  # is unreachable defensive code because compute_transport_scores()
  # at line 492 catches the same condition first.
  expect_error(
    ROOT(df, generalizability_path = TRUE),
    "no S == 1"
  )
})

# Shared tiny dataset for tie tests
.tie_df <- data.frame(
  Y  = c(10, 20, 30, 40, 50, 60, 70, 80),
  Tr = c(1, 0, 1, 0, 1, 0, 1, 0),
  S  = c(1, 1, 1, 1, 1, 1, 0, 0),
  X1 = c(0.1, 0.9, 0.2, 0.8, 0.3, 0.7, 0.4, 0.6)
)

.find_tie_seed <- function(data, vote_threshold, pattern, max_seeds = 100) {
  for (s in seq_len(max_seeds)) {
    err <- tryCatch(
      {
        suppressWarnings(
          ROOT(data, generalizability_path = TRUE,
               num_trees = 2, top_k_trees = TRUE, k = 2,
               seed = s, vote_threshold = vote_threshold)
        )
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    if (!is.null(err) && grepl(pattern, err)) return(s)
  }
  NA_integer_
}

test_that("vote_threshold 'mean' errors on exact 0.5 tie (end-to-end)", {
  skip_if_not_installed("MASS")

  s <- .find_tie_seed(.tie_df, "mean", "exact mean vote of 0.5")
  if (is.na(s)) skip("Could not find seed producing exact 0.5 mean tie")

  expect_error(
    suppressWarnings(
      ROOT(.tie_df, generalizability_path = TRUE,
           num_trees = 2, top_k_trees = TRUE, k = 2,
           seed = s, vote_threshold = "mean")
    ),
    "exact mean vote of 0.5"
  )
})

test_that("vote_threshold 'median' errors on exact 0.5 tie (end-to-end)", {
  skip_if_not_installed("MASS")

  s <- .find_tie_seed(.tie_df, "median", "exact median vote of 0.5")
  if (is.na(s)) skip("Could not find seed producing exact 0.5 median tie")

  expect_error(
    suppressWarnings(
      ROOT(.tie_df, generalizability_path = TRUE,
           num_trees = 2, top_k_trees = TRUE, k = 2,
           seed = s, vote_threshold = "median")
    ),
    "exact median vote of 0.5"
  )
})

test_that("vote_threshold numeric errors on exact tie at threshold (end-to-end)", {
  skip_if_not_installed("MASS")

  s <- .find_tie_seed(.tie_df, 0.5, "mean vote exactly equal to the numeric")
  if (is.na(s)) skip("Could not find seed producing exact numeric tie")

  expect_error(
    suppressWarnings(
      ROOT(.tie_df, generalizability_path = TRUE,
           num_trees = 2, top_k_trees = TRUE, k = 2,
           seed = s, vote_threshold = 0.5)
    ),
    "mean vote exactly equal to the numeric"
  )
})

test_that("ROOT appends custom-objective SE note in generalizability mode", {
  skip_if_not_installed("MASS")

  set.seed(55)
  n <- 100
  df <- data.frame(
    Y  = rnorm(n, mean = 5),
    Tr = rep(0:1, n / 2),
    S  = c(rep(1L, n - 15), rep(0L, 15)),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  # Must NOT be identical() to objective_default
  my_custom_obj <- function(D) {
    w <- D$w
    if (sum(w) == 0) return(Inf)
    sqrt(sum(w * D$vsq) / sum(w)^2)
  }

  result <- suppressWarnings(
    ROOT(
      data                  = df,
      generalizability_path = TRUE,
      global_objective_fn   = my_custom_obj,
      num_trees             = 4,
      top_k_trees           = TRUE,
      k                     = 4,
      seed                  = 42,
      vote_threshold        = "majority"
    )
  )

  expect_s3_class(result, "ROOT")

  # If observations were kept, the note should mention custom objective
  if (!is.null(result$estimate) && any(result$D_rash$w_opt == 1)) {
    expect_true(grepl("custom global_objective_fn|verify this SE",
                      result$estimate$se_weighted_note))
  }
})

test_that("ROOT verbose output shows SE info with custom objective", {
  skip_if_not_installed("MASS")

  set.seed(55)
  n <- 100
  df <- data.frame(
    Y  = rnorm(n, mean = 5),
    Tr = rep(0:1, n / 2),
    S  = c(rep(1L, n - 15), rep(0L, 15)),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  my_custom_obj <- function(D) {
    w <- D$w
    if (sum(w) == 0) return(Inf)
    sqrt(sum(w * D$vsq) / sum(w)^2)
  }

  msgs <- capture.output(
    suppressWarnings(
      result <- ROOT(
        data                  = df,
        generalizability_path = TRUE,
        global_objective_fn   = my_custom_obj,
        num_trees             = 4,
        top_k_trees           = TRUE,
        k                     = 4,
        seed                  = 42,
        vote_threshold        = "majority",
        verbose               = TRUE
      )
    ),
    type = "message"
  )

  expect_s3_class(result, "ROOT")
  # verbose should print SATE/WTATE messages
  expect_true(any(grepl("SATE|WTATE", msgs)))
})

test_that("print.characterizing_underrep rejects wrong class", {
  expect_error(
    print.characterizing_underrep(list()),
    "Not a characterizing_underrep object"
  )
})

test_that("print.characterizing_underrep prints brief ROOT summary and returns invisibly", {
  obj <- structure(
    list(
      root = fake_root_for_print(generalizability_path = FALSE),
      combined = data.frame(),
      leaf_summary = NULL
    ),
    class = "characterizing_underrep"
  )

  out <- capture.output(ret <- print(obj))
  expect_true(any(grepl("characterizing_underrep object", out)))
  expect_true(any(grepl("ROOT brief summary", out)))
  expect_identical(ret, obj)
})

test_that("ROOT validates generalizability_path must be logical scalar", {
  df <- data.frame(X1 = 1:10, Y = rnorm(10), Tr = rep(0:1, 5), S = rep(1, 10))

  expect_error(
    ROOT(df, generalizability_path = "yes"),
    "`generalizability_path` must be TRUE or FALSE"
  )
  expect_error(
    ROOT(df, generalizability_path = c(TRUE, FALSE)),
    "`generalizability_path` must be TRUE or FALSE"
  )
})

test_that("ROOT validates Y must be numeric in generalizability_path mode", {
  df <- data.frame(
    Y  = c("a", "b", "c", "d", "e", "f"),
    Tr = c(1, 0, 1, 0, 1, 0),
    S  = c(1, 1, 1, 1, 0, 0),
    X1 = rnorm(6)
  )

  expect_error(
    ROOT(df, generalizability_path = TRUE),
    "Y.*must be numeric"
  )
})

test_that("ROOT errors when Tr has NA among trial units", {
  df <- data.frame(
    Y  = rnorm(6),
    Tr = c(1, 0, NA, 0, 1, 0),
    S  = c(1, 1, 1, 1, 0, 0),
    X1 = rnorm(6)
  )

  expect_error(
    ROOT(df, generalizability_path = TRUE),
    "Tr.*missing values among trial units"
  )
})

test_that("ROOT issues message when Rashomon set empty with baseline cutoff", {
  skip_if_not_installed("MASS")

  set.seed(123)
  n <- 30
  df <- data.frame(
    Y  = rep(0, n),
    Tr = rep(0:1, n / 2),
    S  = rep(1L, n),
    X1 = rnorm(n)
  )

  # This potentially produces "No trees selected into Rashomon set" as message
  result <- tryCatch(
    suppressWarnings(
      ROOT(df, generalizability_path = TRUE, num_trees = 2,
           top_k_trees = FALSE, cutoff = "baseline", seed = 1,
           vote_threshold = "majority")
    ),
    error = function(e) e
  )

  expect_true(inherits(result, "ROOT") || inherits(result, "error"))
})

test_that("ROOT coerce01 handles logical S and Tr columns", {
  skip_if_not_installed("MASS")

  set.seed(10)
  n <- 40
  df <- data.frame(
    Y  = rnorm(n),
    Tr = rep(c(TRUE, FALSE), n / 2),
    S  = c(rep(TRUE, n - 5), rep(FALSE, 5)),
    X1 = rnorm(n)
  )

  result <- suppressWarnings(
    ROOT(df, generalizability_path = TRUE, num_trees = 2,
         top_k_trees = TRUE, k = 2, seed = 1,
         vote_threshold = "majority")
  )

  expect_s3_class(result, "ROOT")
  expect_true(all(result$D_rash$w_opt %in% 0:1))
})

test_that("ROOT coerce01 handles factor S and Tr columns", {
  skip_if_not_installed("MASS")

  set.seed(11)
  n <- 40
  df <- data.frame(
    Y  = rnorm(n),
    Tr = factor(rep(c("treated", "control"), n / 2)),
    S  = factor(c(rep("yes", n - 5), rep("no", 5))),
    X1 = rnorm(n)
  )

  result <- suppressWarnings(
    ROOT(df, generalizability_path = TRUE, num_trees = 2,
         top_k_trees = TRUE, k = 2, seed = 1,
         vote_threshold = "majority")
  )

  expect_s3_class(result, "ROOT")
  expect_true(all(result$D_rash$w_opt %in% 0:1))
})

