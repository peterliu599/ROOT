# helper in tests/test-helper.R (or at top of the test file)
fake_root <- function(D_forest,
                      D_rash = data.frame(),
                      rashomon_set = integer(0),
                      estimate = NULL,
                      f = NULL,
                      w_forest = NULL,
                      global_objective_fn = objective_default,
                      single_sample_mode = NULL) {
  obj <- list(
    D_forest = D_forest,
    D_rash = D_rash,
    rashomon_set = rashomon_set,
    testing_data = D_forest,           # used by diagnostics
    estimate = estimate,
    f = f,
    w_forest = w_forest,
    global_objective_fn = global_objective_fn,
    single_sample_mode = single_sample_mode
  )
  class(obj) <- "ROOT"
  obj
}

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

test_that(".root_is_single_sample respects flag and lX fallback", {
  obj_flag_true  <- list(single_sample_mode = TRUE)
  obj_flag_false <- list(single_sample_mode = FALSE)
  expect_true(ROOT:::.root_is_single_sample(obj_flag_true))
  expect_false(ROOT:::.root_is_single_sample(obj_flag_false))

  # Fallback: all lX NA -> TRUE; any non-NA -> FALSE
  Dna  <- data.frame(lX = c(NA_real_, NA_real_))
  Dmix <- data.frame(lX = c(NA_real_,  0.5))
  obj_na  <- list(D_forest = Dna)
  obj_mix <- list(D_forest = Dmix)
  expect_true(ROOT:::.root_is_single_sample(obj_na))
  expect_false(ROOT:::.root_is_single_sample(obj_mix))
})

test_that(".root_covariate_names excludes v,vsq,S,lX and w_tree_*", {
  Df <- data.frame(
    X1 = 1, X2 = 2, v = 3, vsq = 4, S = 1, lX = NA_real_,
    w_tree_1 = 0L, w_tree_2 = 1L
  )
  obj <- list(D_forest = Df)
  covs <- ROOT:::.root_covariate_names(obj)
  expect_equal(sort(covs), c("X1", "X2"))
})

test_that(".root_baseline_loss uses sqrt(sum(vsq)/n^2) with na.rm", {
  Df <- data.frame(vsq = c(1, 4, NA_real_), X1 = 1:3)
  obj <- list(D_forest = Df)
  got <- ROOT:::.root_baseline_loss(obj)
  expect_equal(got, sqrt(sum(c(1,4), na.rm = TRUE) / (nrow(Df)^2)))
})

test_that(".root_selected_objectives handles empty and names outputs", {
  # Empty forest
  obj_empty <- list(w_forest = list(), rashomon_set = integer(0))
  expect_length(ROOT:::.root_selected_objectives(obj_empty), 0)

  # Some trees missing the field -> NA, but rashomon_set empty -> numeric(0)
  obj_none <- list(
    w_forest = list(list(), list(`local objective` = 0.25)),
    rashomon_set = integer(0)
  )
  expect_length(ROOT:::.root_selected_objectives(obj_none), 0)

  # Select 2nd only; ensure name = "w_tree_2"
  obj_sel <- list(
    w_forest = list(list(`local objective` = 0.5),
                    list(`local objective` = 0.25),
                    list(`local objective` = 0.75)),
    rashomon_set = 2L
  )
  out <- ROOT:::.root_selected_objectives(obj_sel)
  expect_equal(unname(out), 0.25)
  expect_equal(names(out), "w_tree_2")
})

test_that("summary.ROOT prints precomputed estimates (with/without SE)", {
  Df <- data.frame(v = rnorm(3), vsq = 1, S = c(1L,1L,1L), lX = NA_real_,
                   X1 = 1:3, w_tree_1 = 1L)

  obj_with_se <- fake_root(
    D_forest = Df,
    D_rash   = data.frame(w_opt = c(1L,1L,0L)),
    rashomon_set = 1L,
    estimate = list(
      estimand_unweighted = "ATE in RCT",
      value_unweighted    = 0.1,
      se_unweighted       = 0.2,
      estimand_weighted   = "WATE",
      value_weighted      = 0.3,
      se_weighted         = 0.4,
      se_weighted_note    = "ok"
    )
  )
  expect_output(summary(obj_with_se), "ATE in RCT \\(unweighted\\).*SE = 0.200000")
  expect_output(summary(obj_with_se), "WATE \\(weighted\\).*SE = 0.400000")
  expect_output(summary(obj_with_se), "Note: ok")

  obj_no_se <- obj_with_se
  obj_no_se$estimate$se_weighted <- NA_real_
  expect_output(summary(obj_no_se), "WATE \\(weighted\\)   = 0.300000")
})

test_that("summary.ROOT fallback paths: binary vs nonbinary and no-kept", {
  Df <- data.frame(
    v = c(1, 2, 3, 4), vsq = 1, S = c(1L,1L,1L,0L), lX = 0,
    X1 = 0, w_tree_1 = 1L
  )

  # (A) Binary weights with some kept
  objA <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(1L,0L,1L,NA)))
  expect_output(summary(objA), "WTATE \\(weighted\\)")

  # (B) Binary weights but no kept
  objB <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(0L,0L,0L,NA)))
  expect_output(summary(objB), "NA \\(no kept observations\\)")
  expect_output(summary(objB), "SE omitted because no kept observations")

  # (C) Non-binary weights
  objC <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(0.2, 0.8, 0.6, NA)))
  expect_output(summary(objC), "WTATE \\(weighted\\)   = ")
  expect_output(summary(objC), "SE omitted: non-binary w_opt detected")

  # (D) Custom objective -> advisory note
  custom_obj <- function(D) mean(D$vsq)
  objD <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(1L,0L,1L,NA)),
                    global_objective_fn = custom_obj)
  expect_output(summary(objD), "custom global_objective_fn")
})

test_that("summary.ROOT diagnostics block prints baseline/selected loss and kept %", {
  Df <- data.frame(v = rnorm(5), vsq = c(1,2,3,4,5), S = c(1,1,1,1,0), lX = 0,
                   X1 = 0, w_tree_1 = 1L, w_tree_2 = 0L)
  obj <- fake_root(
    D_forest = Df,
    D_rash   = data.frame(w_opt = c(1,0,1,0,NA)),
    rashomon_set = c(1L, 2L),
    w_forest = list(
      list(`local objective` = 0.33),
      list(`local objective` = 0.10)
    )
  )
  out <- capture.output(summary(obj))
  expect_true(any(grepl("^  Baseline loss:\\s+\\d", out)))
  expect_true(any(grepl("^  Selected loss:\\s+min/median = ", out)))
  expect_true(any(grepl("^  Kept \\(w_opt=1\\):", out)))
})

test_that("print.ROOT covers weighted/unweighted branches and kept%", {
  Df <- data.frame(v = c(5,6,7), vsq = 1, S = c(1L,1L,0L), lX = 0,
                   X1 = 0, w_tree_1 = 1L)
  obj <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(1L, 0L, NA)),
                   rashomon_set = 1L)
  expect_output(print(obj), "TATE \\(unweighted\\)|ATE in RCT \\(unweighted\\)")
  expect_output(print(obj), "Kept \\(w_opt=1\\):")
})

test_that("plot.ROOT handles x$f NULL and plots when present", {
  # NULL branch still returns invisibly with a message
  obj_null <- fake_root(D_forest = data.frame(), f = NULL)
  expect_invisible(plot(obj_null))  # message is fine

  skip_if_not_installed("rpart")
  skip_if_not_installed("rpart.plot")

  # Use a classification tree so default extra=109 is valid
  set.seed(1)
  df  <- data.frame(y = factor(sample(c(0, 1), 60, TRUE)),
                    X1 = runif(60), X2 = runif(60))
  fit <- rpart::rpart(y ~ X1 + X2, data = df, method = "class")
  obj <- fake_root(D_forest = df, f = fit)

  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off(), add = TRUE)

  # We only care that it doesn't error; suppress incidental layout warnings.
  expect_error(suppressWarnings(plot(obj)), NA)
})
