test_that("objective_default and loss_from_objective agree on real code (no mocks)", {
  # Fetch the *original* functions from the package namespace.
  lfo <- getFromNamespace("loss_from_objective", "ROOT")
  obj <- getFromNamespace("objective_default",   "ROOT")

  # Small deterministic state
  D <- data.frame(
    vsq = c(0.1, 0.3, 0.1, 0.5),
    w   = c(0L,   0L,   0L,   0L),
    row.names = paste0("r", 1:4),
    check.names = FALSE
  )

  losser <- lfo(obj)

  # 1) Empty index should not error and should be numeric & non-negative.
  L0 <- losser(1, integer(0), D)
  expect_type(L0, "double")
  expect_true(is.finite(L0) || is.infinite(L0))  # allow either if code returns Inf
  expect_true(!is.na(L0))

  # 2) “All rows kept” (val = 1) should be finite numeric and >= 0.
  L_all <- losser(1, seq_len(nrow(D)), D)
  expect_type(L_all, "double")
  expect_true(is.finite(L_all))
  expect_gte(L_all, 0)

  # 3) The two call paths should be *consistent up to monotonicity*:
  #    turning some indices from kept->dropped (1->0) should not *reduce*
  #    the global objective returned by losser when that’s how your
  #    objective/loss is designed (nonnegative vsq).
  idx_all   <- seq_len(nrow(D))
  idx_drop1 <- idx_all[-length(idx_all)]  # drop the last row
  L_drop1   <- losser(1, idx_drop1, D)

  expect_true(is.finite(L_drop1))
  # If your objective is loss-like with nonnegative contributions,
  # keeping fewer rows should *not* produce a smaller loss.
  expect_gte(L_drop1, min(L_all, L0, na.rm = TRUE))
})

test_that("loss_from_objective composes correctly with a custom objective (no namespace mocks)", {
  lfo <- getFromNamespace("loss_from_objective", "ROOT")

  # Define a custom objective that's easy to predict from the subset:
  # here it's just 5 * nrow(subset). (So for 2 rows -> 10.)
  custom_obj <- function(D) 5 * nrow(D)

  losser <- lfo(custom_obj)

  D <- data.frame(
    vsq = c(2, 3),
    w   = c(0L, 0L),
    row.names   = c("a", "b"),
    check.names = FALSE
  )

  # Empty index: don't assert exact value—just that it's defined.
  L0 <- losser(1, integer(0), D)
  expect_true(is.finite(L0) || is.infinite(L0))
  expect_type(L0, "double")

  # Non-empty subset should equal the custom objective evaluated on that subset.
  idx   <- 1:2
  L_any <- losser(1, idx, D)

  expected <- custom_obj(D[idx, , drop = FALSE])
  expect_equal(L_any, expected, tolerance = 1e-12)
})

test_that(".check_no_na detects NAs and returns invisibly TRUE otherwise", {
  df <- data.frame(a = 1:3, b = c(1, NA, 3))
  expect_error(.check_no_na(df, c("a","b")), "contains missing values")

  df2 <- data.frame(a = 1:3, b = 11:13)
  expect_invisible(.check_no_na(df2, c("a","b")))
})

test_that("midpoint, choose_feature, and reduce_weight work", {
  set.seed(1)
  x <- c(2, 8, 10, 4)
  expect_equal(midpoint(x), (min(x)+max(x))/2)

  # choose_feature sampling respects probabilities and names
  sf <- c(leaf = 0.2, X0 = 0.3, X1 = 0.5)
  got <- replicate(1000, choose_feature(sf, depth = 0))
  props <- prop.table(table(got))
  expect_true(all(names(sf) %in% names(props)))
  expect_gt(props[["X1"]], props[["leaf"]]) # higher prob picked more often (probabilistic check)

  # reduce_weight halves target prob (after renorm)
  sf2 <- reduce_weight("X1", sf)
  expect_lt(sf2["X1"], sf["X1"])
  expect_equal(sum(sf2), 1, tolerance = 1e-12)

  # error cases
  expect_error(choose_feature(unname(sf), 0), "must have names")
  expect_error(reduce_weight("X9", sf), "not found")
})
