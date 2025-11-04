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
