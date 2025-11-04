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
