test_that("train() fits nuisance models and estimate() returns finite v,a,b", {
  # Build a small two-sample dataset
  sim <- get_data(n = 400, seed = 123)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)  # already inside sim$data

  # Some S==1 rows need non-NA Tr & Yobs — get_data provides that structure
  trn <- train(df, outcome = "Yobs", treatment = "Tr", sample = "S")
  expect_true(is.numeric(trn$pi) && trn$pi > 0 && trn$pi < 1)
  expect_s3_class(trn$pi_m, "glm")
  expect_s3_class(trn$e_m,  "glm")

  est <- estimate(df, outcome = "Yobs", treatment = "Tr", sample = "S",
                  pi = trn$pi, pi_m = trn$pi_m, e_m = trn$e_m)
  expect_length(est$v, nrow(df))
  expect_length(est$a, nrow(df))
  expect_length(est$b, nrow(df))
  expect_true(all(is.finite(est$v[ df$S == 1 ])))
})

test_that("estimate_dml returns aligned df_v and data2 for S==1 rows", {
  sim <- get_data(n = 500, seed = 321)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  out <- estimate_dml(df, outcome = "Yobs", treatment = "Tr", sample = "S", crossfit = 5)
  expect_true(all(c("te","a","b","te_sq","a_sq","primary_index") %in% names(out$df_v)))
  expect_equal(nrow(out$data2), nrow(out$df_v))  # aligned
})

test_that("estimate_dml_single works in single-sample mode", {
  # Emulate single-sample by filtering S==1 and dropping S
  sim <- get_data(n = 300, seed = 777)
  dfS <- subset(sim$data, S == 1)
  dfS$Yobs <- dfS$Yobs  # explicit

  out <- estimate_dml_single(data = dfS[, c(grep("^X", names(dfS), value = TRUE), "Tr", "Yobs")],
                             outcome = "Yobs", treatment = "Tr", crossfit = 4)
  expect_true(all(c("te","a","b","te_sq","a_sq","primary_index") %in% names(out$df_v)))
  expect_true(all(out$df_v$b == 1))
})

test_that("DML errors are thrown for degenerate or NA cases", {
  sim <- get_data(n = 200, seed = 55)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  # All S==1 or S==0 → train() should error
  df_bad <- df; df_bad$S <- 1L
  expect_error(train(df_bad, "Yobs", "Tr", "S"), "no variation")

  # NA in covariates should error
  df_na <- df; nm <- grep("^X", names(df_na), value = TRUE)[1]; df_na[1, nm] <- NA
  expect_error(train(df_na, "Yobs", "Tr", "S"), "Covariates contain NA")

  # estimate() invalid pi
  trn <- train(df, "Yobs", "Tr", "S")
  expect_error(estimate(df, "Yobs", "Tr", "S", pi = 1, pi_m = trn$pi_m, e_m = trn$e_m),
               "must be between 0 and 1")
})

test_that("stratified_kfold works correctly", {
  S <- c(rep(0, 10), rep(1, 10))
  folds <- stratified_kfold(S, K = 5)
  expect_length(folds, 5)
  # Check all indices are present exactly once
  expect_equal(sort(unlist(folds)), 1:20)

  # Test K > N adjustment
  expect_warning(folds_big <- stratified_kfold(S[1:5], K = 10), "Requested K")
  expect_length(folds_big, 5)

  # Test error inputs
  expect_error(stratified_kfold(S, K = 0), "positive integer")
  expect_error(stratified_kfold(data.frame(S=S)), "vector or factor")
})

test_that("train function validates inputs and fits models", {
  # Happy path
  train_idx <- sample(nrow(dat_2s), 150)
  train_data <- dat_2s[train_idx, ]
  res <- train(train_data, outcome="Yobs", treatment="Tr", sample="S")

  expect_type(res$pi, "double")
  expect_s3_class(res$pi_m, "glm")
  expect_s3_class(res$e_m, "glm")

  # Input validation errors
  expect_error(train(as.matrix(train_data), "Yobs", "Tr", "S"), "must be a data frame")
  expect_error(train(train_data, "WrongY", "Tr", "S"), "Column 'WrongY' not found")

  # NA handling errors
  dat_na <- train_data
  dat_na$S[1] <- NA
  expect_error(train(dat_na, "Yobs", "Tr", "S"), "`S` contains NA")

  dat_na <- train_data
  dat_na$X0[1] <- NA
  expect_error(train(dat_na, "Yobs", "Tr", "S"), "Covariates contain NA")

  # Lack of variation errors
  dat_no_var <- train_data
  dat_no_var$S <- 1
  expect_error(train(dat_no_var, "Yobs", "Tr", "S"), "has no variation")
})

test_that("estimate function calculates pseudo-outcomes correctly", {
  # Setup training models
  train_idx <- 1:100
  test_idx <- 101:200
  trn <- train(dat_2s[train_idx,], outcome="Yobs", treatment="Tr", sample="S")
  test_data <- dat_2s[test_idx,]

  # Happy path
  est <- estimate(test_data, outcome="Yobs", treatment="Tr", sample="S",
                  pi = trn$pi, pi_m = trn$pi_m, e_m = trn$e_m)

  expect_equal(length(est$v), nrow(test_data))
  expect_equal(length(est$a), nrow(test_data))
  expect_equal(length(est$b), nrow(test_data))
  expect_true(all(is.finite(est$v)))

  # Test probability clamping (hard to trigger naturally, verify logic exists)
  # If we force an extreme prediction, it shouldn't result in Inf weights due to clamping in code
  # Mock a model that predicts near 0
  mock_em <- trn$e_m
  mock_em$coefficients[] <- -100 # forces prob near 0
  est_clamp <- estimate(test_data, outcome="Yobs", treatment="Tr", sample="S",
                        pi = trn$pi, pi_m = trn$pi_m, e_m = mock_em)
  expect_true(all(is.finite(est_clamp$a)))

  # Input validation errors
  expect_error(estimate(test_data, "Yobs", "Tr", "S", pi = 1.5, trn$pi_m, trn$e_m), "Invalid `pi`")
  expect_error(estimate(test_data, "Yobs", "Tr", "S", pi = 0.5, "not_a_model", trn$e_m), "must be model objects")
})

test_that("estimate_dml performs cross-fitting correctly", {
  # Happy path (two-sample)
  res <- estimate_dml(dat_2s, outcome="Yobs", treatment="Tr", sample="S", crossfit = 3)
  expect_s3_class(res$df_v, "data.frame")
  expect_true(all(c("te", "a", "b") %in% names(res$df_v)))
  # Should only contain S=1 rows
  expect_equal(nrow(res$df_v), sum(dat_2s$S == 1))

  # Validation errors
  expect_error(estimate_dml(dat_2s, "Yobs", "Tr", "S", crossfit = 1), "integer >= 2")

  dat_no_var <- dat_2s
  dat_no_var$S <- 1
  expect_error(estimate_dml(dat_no_var, "Yobs", "Tr", "S"), "has no variation")
})

test_that("Single sample DML functions work", {
  # train_single
  trn_s <- train_single(dat_1s[1:100,], "Yobs", "Tr")
  expect_s3_class(trn_s$e_m, "glm")

  # estimate_single
  est_s <- estimate_single(dat_1s[101:200,], "Yobs", "Tr", trn_s$e_m)
  expect_equal(est_s$b, rep(1, 100)) # b should be 1 in single sample

  # estimate_dml_single
  res_s <- estimate_dml_single(dat_1s, "Yobs", "Tr", crossfit = 3)
  expect_equal(nrow(res_s$df_v), nrow(dat_1s))
  expect_true(all(res_s$df_v$b == 1))
})
