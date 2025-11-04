test_that("gen_XY produces sensible shapes and values", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  out <- gen_XY(n = 200, seed = 123)
  expect_s3_class(out$X, "data.frame")
  expect_s3_class(out$Y, "data.frame")
  expect_equal(nrow(out$X), 200L)
  expect_equal(nrow(out$Y), 200L)
  expect_true(all(c("Y0","Y1") %in% names(out$Y)))
})

test_that("gen_S lowers inclusion inside the (0.5,1)×(0.5,1) box", {
  XY <- gen_XY(n = 400, seed = 7)
  X  <- XY$X
  S  <- gen_S(X, seed = 7)
  in_box  <- (X$X0 > 0.5 & X$X0 < 1) & (X$X1 > 0.5 & X$X1 < 1)
  expect_true(mean(S$S[in_box] == 1) < mean(S$S[!in_box] == 1))
})

test_that("gen_T mixes randomized and observational treatment", {
  XY <- gen_XY(n = 300, seed = 11)
  S  <- gen_S(XY$X, seed = 11)
  T  <- gen_T(XY$X, S, seed = 11)
  expect_s3_class(T$Tr, "data.frame")
  expect_equal(nrow(T$Tr), nrow(XY$X))
  expect_equal(length(T$pi), nrow(XY$X))
  # probabilities in (0,1)
  expect_true(all(T$pi > 0 & T$pi < 1))
})

test_that("get_data returns aligned frames and observed outcomes", {
  sim <- get_data(n = 250, seed = 99)
  expect_equal(nrow(sim$data), 250L)
  expect_equal(nrow(sim$Y),    250L)
  expect_true(all(c("S","Tr","Yobs") %in% names(sim$data)))
  # Check Yobs relation with Y0/Y1
  idx <- 1:10
  with(sim, {
    Y0 <- Y$Y0[idx]; Y1 <- Y$Y1[idx]; Tr <- data$Tr[idx]
    expect_equal(data$Yobs[idx], Tr*Y1 + (1-Tr)*Y0, tolerance = 1e-8)
  })
})
