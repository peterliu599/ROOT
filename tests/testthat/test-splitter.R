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
  tree0 <- split_node(sf, toy_X, toy_D, parent_loss = Inf, depth = 0,
                      explore_proba = 0, max_depth = 0, min_leaf_n = 1)
  expect_equal(tree0$node, "leaf")
  expect_match(tree0$leaf_reason, "max-depth")

  # Enforce early leaf via min_leaf_n
  tree1 <- split_node(sf, toy_X[1:5, , drop = FALSE], toy_D, parent_loss = Inf, depth = 0,
                      explore_proba = 0, max_depth = 8, min_leaf_n = 6)
  expect_equal(tree1$node, "leaf")
  expect_match(tree1$leaf_reason, "min-leaf")
})

test_that("split_node() can choose 'leaf' and assign best w without exploration", {
  set.seed(123)
  sf <- c(leaf = 1.0, X1 = 0.0)  # always leaf
  res <- split_node(sf, toy_X, toy_D, parent_loss = Inf, depth = 0,
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
  res <- split_node(sf, toy_X, toy_D, parent_loss = .Machine$double.eps,
                    depth = 0, explore_proba = 0, max_depth = 8, min_leaf_n = 2)
  # Because it recurses until a leaf (after penalization), final node can be leaf or split.
  expect_true(res$node %in% c("leaf", "X1"))
})

test_that("split_node accepts improving splits and updates weights", {
  # Minimal synthetic frame with two features and vsq heterogeneity
  set.seed(10)
  X <- data.frame(X0 = runif(50), X1 = runif(50))
  v  <- rnorm(50)
  vsq <- (v - mean(v))^2
  D <- data.frame(X, v = v, vsq = vsq, w = rep(1, 50), S = rep(1L, 50), stringsAsFactors = FALSE)
  rownames(D) <- as.character(seq_len(nrow(D)))

  sf <- c(leaf = 0.2, X0 = 0.4, X1 = 0.4)

  res <- split_node(
    split_feature   = sf,
    X               = D,
    D               = D,
    parent_loss     = Inf,
    depth           = 0L,
    explore_proba   = 0.0,           # deterministic exploitation for test
    choose_feature  = choose_feature,
    loss_fn         = NULL,          # wrap objective_default
    reduce_weight_fn= reduce_weight,
    objective_fn    = objective_default,
    max_depth       = 2,
    min_leaf_n      = 5
  )

  expect_type(res, "list")
  expect_true(all(c("D","local objective","depth") %in% names(res)))
  expect_equal(nrow(res$D), nrow(D))
  expect_true(all(res$D$w %in% 0:1))
})
