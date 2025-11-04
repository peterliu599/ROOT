#' Recursive split builder for weighted tree (internal function)
#'
#' Recursively builds a weighted decision tree to optimize a global objective,
#' using an exploration/exploitation trade-off. Internal; used by ROOT().
#'
#' @param split_feature Named numeric vector of feature selection probabilities (must include "leaf").
#' @param X Data frame of current observations (includes candidate split feature columns; may include a working copy of weights `w`).
#' @param D Data frame representing the global state (must include columns `w` and `vsq`; row names align to observations).
#' @param parent_loss Numeric, the loss value of the parent node (used to decide if a split improves the objective).
#' @param depth Integer, current tree depth.
#' @param explore_proba Numeric, the probability (between 0 and 1) of flipping the exploit choice at a leaf.
#' @param choose_feature_fn Function to choose next feature (default `choose_feature`).
#' @param reduce_weight_fn Function to penalize last-tried feature on rejected split (default `reduce_weight`).
#' @param global_objective_fn Function `function(D) -> numeric` scoring the **entire** state.
#' @param max_depth Integer max depth (stop and make leaf at this depth).
#' @param min_leaf_n Integer min rows to attempt a split; else make leaf.
#' @param log_fn Function for logging; default no-op.
#' @param max_rejects_per_node Safety budget of rejected splits before forcing a leaf.
#' @return A list representing the (sub)tree; includes updated `D` and `local objective`.
#' @importFrom stats rbinom
split_node <- function(split_feature, X, D, parent_loss, depth,
                       explore_proba = 0.05,
                       choose_feature_fn = choose_feature,   # (renamed arg for clarity)
                       reduce_weight_fn = reduce_weight,
                       global_objective_fn = objective_default,
                       max_depth = 8,
                       min_leaf_n = 5,
                       log_fn = function(...) {},
                       max_rejects_per_node = 1000) {
  # Derive micro-evaluator: loss(val, indices, D) from global objective
  loss_fn <- loss_from_objective(global_objective_fn)

  # Helpers
  .log <- function(fmt, ...) log_fn(sprintf(fmt, ...))
  nearly_leq <- function(a,b,tol_abs=1e-12,tol_rel=1e-12) a <= b + tol_abs + tol_rel*max(1,abs(b))

  make_leaf <- function(X_sub, D_sub, depth_cur, reason, fj = NA_character_, cj = NA_real_) {
    losses <- c(loss_fn(0, rownames(X_sub), D_sub), loss_fn(1, rownames(X_sub), D_sub))
    w_exploit <- which.min(losses) - 1
    w_explore <- stats::rbinom(1, 1, 0.5)
    explore_flip <- stats::rbinom(1, 1, explore_proba)
    final_w <- if (explore_flip == 1) w_explore else w_exploit

    .log("[depth=%d] LEAF (%s): feature=%s, cut=%s, n=%d, losses={%.4f, %.4f}, w=%d",
         depth_cur, reason, as.character(fj),
         ifelse(is.na(cj), "NA", sprintf("%.4f", cj)),
         nrow(X_sub), losses[1], losses[2], final_w)

    if (nrow(X_sub) > 0) {
      D_sub[rownames(X_sub), "w"] <- final_w
      X_sub$w <- final_w
    }
    list(
      node = "leaf", w = final_w,
      `local objective` = min(losses),
      depth = depth_cur, D = D_sub,
      leaf_reason = reason, feature = fj, cut = cj
    )
  }

  # Stopping
  if (depth >= max_depth) return(make_leaf(X, D, depth, reason = "max-depth"))
  if (nrow(X) <= min_leaf_n) return(make_leaf(X, D, depth, reason = "min-leaf"))

  fj <- choose_feature_fn(split_feature, depth)
  if (identical(fj, "leaf")) return(make_leaf(X, D, depth, reason = "feature==leaf", fj = fj))

  cj <- midpoint(X[[fj]])
  X_left  <- X[X[[fj]] <= cj, , drop = FALSE]
  X_right <- X[X[[fj]] >  cj, , drop = FALSE]
  if (nrow(X_left) == 0 || nrow(X_right) == 0) {
    return(make_leaf(X, D, depth, reason = "empty-child", fj = fj, cj = cj))
  }

  # Evaluate best decision on children leaves
  loss_left  <- c(loss_fn(0, rownames(X_left),  D), loss_fn(1, rownames(X_left),  D))
  loss_right <- c(loss_fn(0, rownames(X_right), D), loss_fn(1, rownames(X_right), D))
  min_left   <- min(loss_left)
  min_right  <- min(loss_right)
  new_loss   <- (nrow(X_left) * min_left + nrow(X_right) * min_right) / nrow(X)

  if (nearly_leq(new_loss, parent_loss)) {
    .log("[depth=%d] SPLIT: feature=%s, cut=%.4f | n_left=%d, n_right=%d | new_loss=%.6f, parent_loss=%.6f",
         depth, fj, cj, nrow(X_left), nrow(X_right), new_loss, parent_loss)
    .log("    left losses={%.6f, %.6f}, right losses={%.6f, %.6f}",
         loss_left[1], loss_left[2], loss_right[1], loss_right[2])

    w_left  <- which.min(loss_left)  - 1
    w_right <- which.min(loss_right) - 1
    D[rownames(X_left),  "w"] <- w_left;  X_left$w  <- w_left
    D[rownames(X_right), "w"] <- w_right; X_right$w <- w_right

    # Recurse (random order to break symmetry)
    if (stats::rbinom(1, 1, 0.5) == 1) {
      left_res  <- split_node(split_feature, X_left,  D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- left_res$D
      right_res <- split_node(split_feature, X_right, D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- right_res$D
    } else {
      right_res <- split_node(split_feature, X_right, D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- right_res$D
      left_res  <- split_node(split_feature, X_left,  D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- left_res$D
    }

    local_obj <- global_objective_fn(D); if (!is.finite(local_obj) || is.nan(local_obj)) local_obj <- Inf
    return(list(
      node = fj, split = cj,
      left_tree  = left_res, right_tree = right_res,
      `local objective` = local_obj, depth = depth, D = D
    ))
  }

  # Rejection loop with budget (stack-safe)
  .log("[depth=%d] REJECTED: feature=%s, cut=%.4f | new_loss=%.16f > parent_loss=%.16f",
       depth, fj, cj, new_loss, parent_loss)

  rej_sf <- split_feature
  for (attempt in seq_len(max_rejects_per_node)) {
    # Penalize the last tried feature
    rej_sf <- reduce_weight_fn(fj, rej_sf)

    # If "leaf" is now effectively certain, bail out as leaf
    if ("leaf" %in% names(rej_sf) && rej_sf["leaf"] >= 1 - 1e-8) {
      return(make_leaf(X, D, depth, reason = "forced-leaf-after-rejects", fj = fj, cj = cj))
    }

    # Choose a new feature and recompute children
    fj <- choose_feature_fn(rej_sf, depth)
    if (identical(fj, "leaf")) {
      return(make_leaf(X, D, depth, reason = "feature==leaf-after-rejects", fj = fj))
    }

    cj <- midpoint(X[[fj]])
    X_left  <- X[X[[fj]] <= cj, , drop = FALSE]
    X_right <- X[X[[fj]] >  cj, , drop = FALSE]
    if (nrow(X_left) == 0 || nrow(X_right) == 0) {
      return(make_leaf(X, D, depth, reason = "empty-child-after-rejects", fj = fj, cj = cj))
    }

    loss_left  <- c(loss_fn(0, rownames(X_left),  D), loss_fn(1, rownames(X_left),  D))
    loss_right <- c(loss_fn(0, rownames(X_right), D), loss_fn(1, rownames(X_right), D))
    min_left   <- min(loss_left); min_right <- min(loss_right)
    new_loss   <- (nrow(X_left) * min_left + nrow(X_right) * min_right) / nrow(X)

    if (nearly_leq(new_loss, parent_loss)) {
      .log("[depth=%d] SPLIT(after %d rejects): feature=%s, cut=%.4f | n_left=%d, n_right=%d | new_loss=%.6f, parent_loss=%.6f",
           depth, attempt, fj, cj, nrow(X_left), nrow(X_right), new_loss, parent_loss)

      w_left  <- which.min(loss_left)  - 1
      w_right <- which.min(loss_right) - 1
      D[rownames(X_left),  "w"] <- w_left;  X_left$w  <- w_left
      D[rownames(X_right), "w"] <- w_right; X_right$w <- w_right

      if (stats::rbinom(1, 1, 0.5) == 1) {
        left_res  <- split_node(rej_sf, X_left,  D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- left_res$D
        right_res <- split_node(rej_sf, X_right, D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- right_res$D
      } else {
        right_res <- split_node(rej_sf, X_right, D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- right_res$D
        left_res  <- split_node(rej_sf, X_left,  D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- left_res$D
      }

      local_obj <- global_objective_fn(D); if (!is.finite(local_obj) || is.nan(local_obj)) local_obj <- Inf
      return(list(
        node = fj, split = cj,
        left_tree  = left_res, right_tree = right_res,
        `local objective` = local_obj, depth = depth, D = D
      ))
    }

    .log("[depth=%d] REJECTED (%d/%d): feature=%s, cut=%.4f | new_loss=%.16f > parent_loss=%.16f",
         depth, attempt, max_rejects_per_node, fj, cj, new_loss, parent_loss)
  }

  # Budget exhausted -> become a leaf
  make_leaf(X, D, depth, reason = "max-rejects", fj = fj, cj = cj)
}

#' Randomly choose a split feature based on provided probabilities
#'
#' Given a probability distribution over features (and possibly a "leaf" option), selects one feature at random according to those probabilities.
#'
#' @param split_feature A named numeric vector of feature selection probabilities. Names should correspond to feature IDs (and may include a special "leaf" entry).
#' @param depth Current tree depth (an integer, used for parity with Python implementation but not affecting probabilities in this implementation).
#' @return A single feature name (or "leaf") chosen randomly according to the provided probability weights.
#' @note The factor \eqn{2^{(0 * depth / 4)}} present in the code is effectively 1 (no effect on the first element's weight) and is included only for parity with an equivalent Python implementation. All probabilities are normalized to sum to 1 before sampling.
choose_feature <- function(split_feature, depth) {
  # Input validation
  if (!is.numeric(split_feature) || length(split_feature) == 0) {
    stop("`split_feature` must be a non-empty numeric vector of probabilities.", call. = FALSE)
  }
  if (is.null(names(split_feature))) {
    stop("`split_feature` must have names for each feature (use 'leaf' for the special leaf option if applicable).", call. = FALSE)
  }
  if (!is.numeric(depth) || length(depth) != 1) {
    stop("`depth` must be a numeric scalar (current tree depth).", call. = FALSE)
  }

  # New strict probability checks
  if (anyNA(split_feature)) {
    stop("`split_feature` contains NA values.", call. = FALSE)
  }
  if (any(split_feature < 0)) {
    bad <- names(split_feature)[split_feature < 0]
    stop(sprintf(
      "All probabilities must be >= 0. Negative entries found in: %s",
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }

  split_prob <- as.numeric(split_feature)
  # Adjust the first element probability by the (inactive) depth factor (kept for parity)
  if (length(split_prob) > 0) {
    split_prob[1] <- split_prob[1] * (2^(0 * depth / 4))
  }
  # Normalize to a probability distribution
  total <- sum(split_prob)
  if (!is.finite(total)) {
    stop("Invalid probabilities: sum is non-finite.", call. = FALSE)
  }
  split_prob <- split_prob / total

  # Perform random sampling
  choices <- names(split_feature)
  chosen <- sample(choices, size = 1, replace = TRUE, prob = split_prob)
  return(chosen)
}

#' Reduce a feature's selection weight by half and renormalize
#'
#' Lowers the probability weight of a given feature by 50%, and then re-normalizes the entire probability vector.
#'
#' @param fj A feature name (character string) present in the names of `split_feature`.
#' @param split_feature A named numeric vector of probabilities for features (as used in splitting).
#' @return A numeric vector of the same length as `split_feature`, giving the updated probabilities that sum to 1.
#' @details This is typically used when a particular feature split was rejected; the feature's probability is halved to reduce its chance of being chosen again immediately, encouraging exploration of other features. If `fj` is "leaf", its weight is also halved similarly.
reduce_weight <- function(fj, split_feature) {
  # Input validation
  if (!is.character(fj) || length(fj) != 1) {
    stop("`fj` must be a single feature name (string).", call. = FALSE)
  }
  if (!is.numeric(split_feature) || is.null(names(split_feature))) {
    stop("`split_feature` must be a named numeric vector of probabilities.", call. = FALSE)
  }
  if (!(fj %in% names(split_feature))) {
    stop("Feature '", fj, "' not found in names(split_feature).", call. = FALSE)
  }
  # New strict probability checks
  if (anyNA(split_feature)) {
    stop("`split_feature` contains NA values.", call. = FALSE)
  }
  if (any(split_feature < 0)) {
    bad <- names(split_feature)[split_feature < 0]
    stop(sprintf(
      "All probabilities must be >= 0. Negative entries found in: %s",
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }

  total_before <- sum(split_feature)
  if (!is.finite(total_before)) {
    stop("Invalid probabilities: sum is non-finite.", call. = FALSE)
  }

  # Reduce the specified feature's weight by half
  split_feature[fj] <- split_feature[fj] / 2
  # Renormalize to sum to 1
  split_feature <- split_feature / sum(split_feature)
  return(split_feature)
}

#' Compute the midpoint of a numeric vector
#'
#' Calculates the midpoint defined as \eqn{( \max(X) + \min(X) ) / 2}, ignoring any NA values.
#'
#' @param X A numeric vector.
#' @return A numeric scalar giving the midpoint of the finite values in `X`. If `X` is empty or has no finite values, `NA` is returned.
midpoint <- function(X) {
  # Input validation
  if (!is.numeric(X)) {
    stop("`X` must be a numeric vector.", call. = FALSE)
  }
  if (length(X) == 0) {
    warning("Empty vector provided to midpoint(). Returning NA.", call. = FALSE)
    return(NA_real_)
  }
  check_no_na(X, colnames(X))
  # Remove NAs for calculation
  rng <- range(X, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) {
    # All values are NA or infinite
    warning("No finite values in `X`. Returning NA.", call. = FALSE)
    return(NA_real_)
  }
  return((rng[1] + rng[2]) / 2)
}

#' Fit a shallow decision tree to characterize learned weights `w`
#'
#' Trains a classification tree on the covariates `X` to predict the binary membership `w`.
#' This provides an interpretable summary of how the weighted subgroup can be distinguished by `X`.
#'
#' @param X A data frame of covariates (features).
#' @param w A binary vector (0/1 or a factor with two levels) indicating class membership for each observation (e.g., whether an observation is in the selected subgroup).
#' @param max_depth Integer, the maximum tree depth (default 3).
#' @return An `rpart` object representing the fitted decision tree.
#' @details The tree is grown using the Gini index (classification) and is not pruned (complexity parameter `cp = 0`), relying solely on `max_depth` to control complexity. This mirrors the default behavior of scikit-learn's `DecisionTreeClassifier(max_depth=...)`.
#'   If `w` is not already a factor, it will be converted internally. The tree's rules can be interpreted to understand which covariates (and what splits) best separate the two classes defined by `w`.
characterize_tree <- function(X, w, max_depth = 3) {
  # Input validation
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame of covariates.", call. = FALSE)
  }
  check_no_na(X, colnames(X))
  if (length(w) != nrow(X)) {
    stop("Length of `w` must equal the number of rows in `X`.", call. = FALSE)
  }
  # Coerce w to factor and check it has two levels
  w_factor <- as.factor(w)
  if (nlevels(w_factor) != 2) {
    stop("`w` must have exactly two classes (binary).", call. = FALSE)
  }

  # Prepare data for rpart
  df <- data.frame(w = w_factor, X)
  # Fit classification tree
  fit <- rpart::rpart(
    w ~ .,
    data = df,
    method = "class",
    parms = list(split = "gini"),
    control = rpart::rpart.control(maxdepth = max_depth, cp = 0.0, minsplit = 2, minbucket = 1),
    model = TRUE
  )
  return(fit)
}
