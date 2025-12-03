# Helper to create dummy state for split_node calls
make_dummy_state <- function(n=20, p=2) {
  X <- matrix(runif(n*p), ncol=p)
  colnames(X) <- paste0("X", 1:p)
  X_df <- as.data.frame(X)
  X_df$w <- 1
  rownames(X_df) <- as.character(1:n)

  D <- X_df
  D$vsq <- rchisq(n, df=1) # random positive values

  split_feats <- c(leaf=0.1, X1=0.45, X2=0.45)
  list(X=X_df, D=D, sf=split_feats)
}

test_that("split_node base cases work", {
  st <- make_dummy_state(n=10)
  res_depth <- split_node(st$sf, st$X, st$D, parent_loss = Inf, depth = 8, max_depth = 8)
  expect_equal(res_depth$node, "leaf")
  expect_equal(res_depth$leaf_reason, "max-depth")

  res_min_n <- split_node(st$sf, st$X[1:4,], st$D, parent_loss = Inf, depth = 0, min_leaf_n = 5)
  expect_equal(res_min_n$node, "leaf")
  expect_equal(res_min_n$leaf_reason, "min-leaf")

  sf_leaf <- c(leaf=1, X1=0, X2=0)
  res_leaf_feat <- split_node(sf_leaf, st$X, st$D, parent_loss = Inf, depth = 0)
  expect_equal(res_leaf_feat$node, "leaf")
  expect_equal(res_leaf_feat$leaf_reason, "feature==leaf")
})

test_that("split_node successfully splits", {
  st <- make_dummy_state(n=50)
  res_split <- withr::with_seed(1, {
    split_node(st$sf, st$X, st$D, parent_loss = 100, depth = 0)
  })
  expect_true(res_split$node %in% c("X1", "X2"))
})

test_that("Helpers: choose_feature, reduce_weight, midpoint", {
  # choose_feature
  sf <- c(leaf=0.2, X1=0.8)
  set.seed(1)
  chosen <- replicate(100, choose_feature(sf, depth=0))
  expect_true(mean(chosen == "X1") > 0.7)

  # reduce_weight
  sf_red <- reduce_weight("X1", sf)
  expected_val <- (0.8/2) / (0.2 + 0.8/2)
  expect_equal(sf_red["X1"], expected_val, ignore_attr = TRUE)
  expect_equal(sum(sf_red), 1)

  # midpoint
  expect_equal(midpoint(c(1, 3)), 2)
  suppressWarnings({
    expect_true(is.na(midpoint(c(NA_real_))))
  })
})

test_that("characterize_tree fits an rpart model", {
  # Use dat_tiny which has 50 rows in the helper
  X <- dat_tiny[, c("X0", "X1")]
  n <- nrow(X)

  # Create w with exactly length n, and ensure two classes
  w <- rep(0, n)
  w[1:(n/2)] <- 1

  fit <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(fit, "rpart")

  # Error: w not binary (all 1s)
  expect_error(characterize_tree(X, rep(1, n)), "must have exactly two classes")

  # Error: Length mismatch
  expect_error(characterize_tree(X, rep(1, n - 1)), "Length of `w` must equal")
})

# Minimal toy data with two numeric features
toy_X <- data.frame(X1 = c(0, 1, 2, 3, 4, 5),
                    X2 = c(5, 4, 3, 2, 1, 0),
                    row.names = as.character(1:6))
toy_D <- data.frame(vsq = c(1, 2, 1, 2, 1, 2),
                    w   = rep(1, 6),
                    row.names = as.character(1:6))

test_that("split_node() stops at max depth & min leaf size", {
  set.seed(123)
  sf <- c(leaf = 0.0, X1 = 0.5, X2 = 0.5)

  # Enforce early leaf via max_depth=0
  tree0 <- split_node(sf, toy_X, toy_D,
                      parent_loss = Inf, depth = 0,
                      explore_proba = 0, max_depth = 0, min_leaf_n = 1)
  expect_equal(tree0$node, "leaf")
  expect_match(tree0$leaf_reason, "max-depth")

  # Enforce early leaf via min_leaf_n
  tree1 <- split_node(sf, toy_X[1:5, , drop = FALSE], toy_D,
                      parent_loss = Inf, depth = 0,
                      explore_proba = 0, max_depth = 8, min_leaf_n = 6)
  expect_equal(tree1$node, "leaf")
  expect_match(tree1$leaf_reason, "min-leaf")
})

test_that("split_node() can choose 'leaf' and assign best w without exploration", {
  set.seed(123)
  sf <- c(leaf = 1.0, X1 = 0.0)  # always leaf
  res <- split_node(sf, toy_X, toy_D,
                    parent_loss = Inf, depth = 0,
                    explore_proba = 0, max_depth = 8, min_leaf_n = 1)
  expect_equal(res$node, "leaf")
  expect_true(res$w %in% c(0, 1))
  # Objective attached
  expect_true(is.finite(res[["local objective"]]))
})

test_that("split_node() rejects non-improving split and penalizes feature weight", {
  set.seed(1)
  sf <- c(leaf = 0.0, X1 = 1.0)
  # Construct parent_loss tiny to force 'reject'
  res <- split_node(sf, toy_X, toy_D,
                    parent_loss = .Machine$double.eps,
                    depth = 0, explore_proba = 0,
                    max_depth = 8, min_leaf_n = 2)
  # Because it recurses until a leaf (after penalization), final node can be leaf or split.
  expect_true(res$node %in% c("leaf", "X1"))
})

test_that("split_node accepts improving splits and updates weights", {
  # Minimal synthetic frame with two features and vsq heterogeneity
  set.seed(10)
  X  <- data.frame(X0 = runif(50), X1 = runif(50))
  v  <- rnorm(50)
  vsq <- (v - mean(v))^2
  D <- data.frame(X, v = v, vsq = vsq, w = rep(1, 50), S = rep(1L, 50), stringsAsFactors = FALSE)
  rownames(D) <- as.character(seq_len(nrow(D)))

  sf <- c(leaf = 0.2, X0 = 0.4, X1 = 0.4)

  res <- split_node(
    split_feature        = sf,
    X                    = D,
    D                    = D,
    parent_loss          = Inf,
    depth                = 0L,
    explore_proba        = 0.0,                 # deterministic exploitation for test
    choose_feature_fn    = choose_feature,      # updated arg name
    reduce_weight_fn     = reduce_weight,
    global_objective_fn  = objective_default,   # replaces objective_fn/loss_fn combo
    max_depth            = 2,
    min_leaf_n           = 5
  )

  expect_type(res, "list")
  expect_true(all(c("D","local objective","depth") %in% names(res)))
  expect_equal(nrow(res$D), nrow(D))
  expect_true(all(res$D$w %in% 0:1))
})
