test_that("train()/estimate() enforce NA rules", {
  d <- data.frame(Y = rnorm(6), Tr = c(1,0,1,0,NA,NA), S = c(1,1,1,1,0,0),
                  X1 = rnorm(6), X2 = runif(6))
  # NA in S
  d2 <- d; d2$S[1] <- NA
  expect_error(train(d2, "Y", "Tr", "S"), "contains NA")
  # NA covariates
  d3 <- d; d3$X2[2] <- NA
  expect_error(train(d3, "Y", "Tr", "S"), "Covariates contain NA")
  # NA in Y among S==1
  d4 <- d; d4$Y[1] <- NA
  expect_error(train(d4, "Y", "Tr", "S"), "contains NA among S==1")
})

test_that("characterize_tree() fails on non-binary w and length mismatch", {
  X <- data.frame(X1 = rnorm(5))
  expect_error(characterize_tree(X, w = c(0,1,2,0,1)), "exactly two classes")
  expect_error(characterize_tree(X, w = c(0,1,1)), "Length of `w`")
})
