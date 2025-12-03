dat_rct <- dat_2s[dat_2s$S == 1, ]
dat_tgt <- dat_2s[dat_2s$S == 0, setdiff(names(dat_2s), c("Yobs", "Tr", "S"))]
covs <- paste0("X", 0:9)

test_that("characterizing_underrep wrapper works", {
  suppressWarnings({
    res <- characterizing_underrep(
      DataRCT = dat_rct[1:50,],
      covariateColName_RCT = covs,
      trtColName_RCT = "Tr",
      outcomeColName_RCT = "Yobs",
      DataTarget = dat_tgt[1:50,],
      covariateColName_TargetData = covs,
      num_trees = 2, seed = 123
      # removed max_depth argument which was causing error
    )
  })

  expect_s3_class(res, "characterizing_underrep")
  expect_s3_class(res$root, "ROOT")
})

test_that("S3 methods for characterizing_underrep work", {
  suppressWarnings({
    res <- characterizing_underrep(
      DataRCT = dat_rct[1:50,], covs, "Tr", "Yobs",
      DataTarget = dat_tgt[1:50,], covs,
      num_trees = 2, seed = 123
    )
  })

  expect_output(print(res), "characterizing_underrep object")
  expect_output(summary(res), "ROOT summary")

  if (!is.null(res$root$f)) {
    expect_silent(grDevices::pdf(file = NULL))
    plot(res)
    grDevices::dev.off()
  }
})

test_that("S3 methods for ROOT work", {
  # Use dat_2s to ensure stable fit
  suppressWarnings({
    res_2s <- ROOT(dat_2s, "Yobs", "Tr", "S", num_trees = 3, seed = 123)
  })

  # Internals check
  expect_false(.root_is_single_sample(res_2s))

  # Methods check
  expect_output(print(res_2s), "TATE.+unweighted")
  expect_output(summary(res_2s), "Rashomon size")
})

test_that("characterize_tree() validates length and binarity of w", {
  skip_if_not_installed("rpart")

  X <- data.frame(X1 = rnorm(30), X2 = runif(30))

  # 1) Length mismatch should error with the length message
  expect_error(
    characterize_tree(X, w = 0:2),
    "Length of `w` must equal the number of rows in `X`"
  )

  # 2) Same length but 3 classes should trip the 'exactly two classes' check
  w3 <- c(0,1,2)[(seq_len(nrow(X)) - 1) %% 3 + 1]   # length 30, classes {0,1,2}
  expect_error(
    characterize_tree(X, w = w3),
    "exactly two classes"
  )

  # 3) Happy path
  set.seed(1)
  w <- sample(c(0,1), size = nrow(X), replace = TRUE)
  fit <- characterize_tree(X, w = w, max_depth = 2)
  expect_s3_class(fit, "rpart")
})

test_that("characterize_tree fits a shallow rpart classifier", {
  skip_if_not_installed("rpart")
  # Simple separable data
  X <- data.frame(x1 = c(0,0,1,1), x2 = c(0,1,0,1))
  w <- c(0,0,1,1)
  fit <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(fit, "rpart")
  expect_true(all(c("frame","where") %in% names(fit)))
})

test_that("characterize_tree errors when w is not binary", {
  X <- data.frame(x1 = rnorm(10))
  w <- rep(2, 10)  # single level
  expect_error(characterize_tree(X, w), "exactly two classes")
})
