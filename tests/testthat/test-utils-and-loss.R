test_that("objective_default and loss behave and are monotone in weights", {
  D <- data.frame(vsq = c(4, 1, 9, 16), w = c(1, 1, 1, 1))
  base <- objective_default(D)

  # If we zero-out a large-vsq row, objective should not increase
  D2 <- D; D2$w[which.max(D2$vsq)] <- 0
  expect_lte(objective_default(D2), base)

  # loss() matches objective_default() when indices cover all rows
  L0 <- loss(0, indices = seq_len(nrow(D)), D = D)  # all excluded → Inf
  expect_true(is.infinite(L0))
  L1 <- loss(1, indices = seq_len(nrow(D)), D = D)  # all included
  expect_equal(L1, objective_default(D), tolerance = 1e-12)

  # Local change reduces objective if we drop a high vsq row
  i_big <- which.max(D$vsq)
  L_drop_big <- loss(0, indices = i_big, D = D)
  expect_lte(L_drop_big, base)
})

test_that("check_no_na detects NAs and returns invisibly TRUE otherwise", {
  df <- data.frame(a = 1:3, b = c(1, NA, 3))
  expect_error(check_no_na(df, c("a","b")), "contains missing values")

  df2 <- data.frame(a = 1:3, b = 11:13)
  expect_invisible(check_no_na(df2, c("a","b")))
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
