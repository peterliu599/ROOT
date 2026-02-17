#' Compute transport influence scores for generalization mode
#'
#' Internal helper used in the \code{generalization} path to construct
#' glm-based, IPW-style scores for transporting trial effects to a target
#' population without double machine learning.
#'
#' The function treats \code{data} as a stacked dataset with a sample
#' indicator \eqn{S} (\code{sample}) taking value 1 in the randomized trial
#' and 0 in the target sample. It proceeds in three steps:
#' \enumerate{
#'   \item Fit a sampling model \eqn{P(S = 1 | X)} using logistic regression
#'         on all rows of \code{data}.
#'   \item Within the trial subset \eqn{S = 1}, fit a treatment model
#'         \eqn{P(T = 1 | X, S = 1)} using logistic regression.
#'   \item For each row, form the density ratio
#'         \eqn{r(X) = P(S = 0 | X) / P(S = 1 | X)}
#'         and compute a Horvitz-Thompson-style transported score
#'         \deqn{v(X, T, Y) = r(X) \left[ \frac{T Y}{e(X)} - \frac{(1 - T) Y}{1 - e(X)} \right],}
#'         where \eqn{e(X) = P(T = 1 | X, S = 1)}.
#' }
#'
#' The resulting score vector \code{v} and its squared version \code{vsq}
#' can be used as pseudo-outcomes for tree-based search over transported
#' treatment effects.
#'
#' @param data A \code{data.frame} containing the outcome, treatment indicator,
#'   sample indicator, and covariates. All columns except \code{outcome},
#'   \code{treatment}, and \code{sample} are treated as covariates \eqn{X}
#'   and must be suitable for use in \code{stats::glm()} with
#'   \code{family = binomial()}.
#' @param outcome A length-1 character string giving the name of the outcome
#'   column in \code{data}. Must be numeric (e.g., a continuous outcome or
#'   a 0/1 indicator).
#' @param treatment A length-1 character string giving the name of the
#'   treatment indicator column in \code{data}. Must be coded 0/1.
#' @param sample A length-1 character string giving the name of the sample
#'   indicator column in \code{data}, with 1 for trial rows and 0 for target
#'   rows.
#' @return A list with two numeric vectors of length \code{nrow(data)}:
#' \describe{
#'   \item{\code{v}}{The transported influence-style score.}
#'   \item{\code{vsq}}{The element-wise square of \code{v}, i.e., \code{v^2}.}
#' }
#' @keywords internal
compute_transport_scores <- function(data, outcome, treatment, sample) {
  covars <- setdiff(colnames(data), c(outcome, treatment, sample))
  if (length(covars) == 0L) {
    stop("compute_transport_scores(): need at least one covariate.", call. = FALSE)
  }

  # P(S = 1 | X)
  formula_s <- stats::as.formula(paste(sample, "~", paste(covars, collapse = "+")))
  model_s   <- stats::glm(formula_s, data = data, family = stats::binomial())
  pi_s      <- stats::predict(model_s, type = "response")

  # P(Tr = 1 | X, S = 1)
  trial_data <- data[data[[sample]] == 1, , drop = FALSE]
  if (nrow(trial_data) == 0L) {
    stop("compute_transport_scores(): no S == 1 (trial) rows for treatment model.",
         call. = FALSE)
  }
  formula_t <- stats::as.formula(paste(treatment, "~", paste(covars, collapse = "+")))
  model_t   <- stats::glm(formula_t, data = trial_data, family = stats::binomial())
  e_x       <- stats::predict(model_t, newdata = data, type = "response")

  Y     <- data[[outcome]]
  T_ind <- data[[treatment]]

  # The influence score v uses 1/l(X) = P(S=0|X) / P(S=1|X) directly, where
  # l(X) = P(S=1|X) / P(S=0|X) is the selection ratio defined in Parikh et al.
  # (2025). 1/l(X) is unbounded when pi_s -> 0, i.e. when some observations are
  # almost never selected into the trial, producing numerically extreme or
  # infinite influence scores v. We therefore floor pi_s at 0.01 before dividing.
  .floor <- 0.01
  n_floored <- sum(pi_s < .floor, na.rm = TRUE)
  if (n_floored > 0L) {
    warning(
      sprintf(
        paste0(
          "%d observation(s) had an estimated trial-inclusion probability ",
          "P(S=1|X) below the numerical floor of 0.01. ",
          "P(S=1|X) has been floored at 0.01 for those observations to prevent ",
          "1/l(X) = P(S=0|X)/P(S=1|X) from becoming unbounded ",
          "(values below 0.01 would produce 1/l(X) > 99), where l(X) = P(S=1|X)/P(S=0|X) ",
          "is the selection ratio defined in Parikh et al. (2025). ",
          "This typically indicates near-perfect separation in the sampling model ",
          "(i.e. certain covariate regions are almost entirely absent from the trial). ",
          "Consider checking the overlap between your trial and target populations."
        ),
        n_floored
      ),
      call. = FALSE
    )
    pi_s <- pmax(pi_s, .floor)
  }

  r_x <- (1 - pi_s) / pi_s

  # Horvitz-Thompson-style transported score
  v   <- r_x * ((T_ind * Y / e_x) - ((1 - T_ind) * Y / (1 - e_x)))
  vsq <- v^2

  list(v = v, vsq = vsq)
}

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
#'
#' @return A list representing the (sub)tree; includes updated `D` and `local objective`.
#' @importFrom stats rbinom
#' @keywords internal
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

  # Evaluate the global objective for all four joint (c_left, c_right) combinations,
  # per Algorithm 1 of Parikh et al. (2025): new_loss <- min_{(c_left,c_right)} loss(w', Dn).
  # The global objective is non-additive across the partition, so each combination
  # must be evaluated jointly rather than independently per child.
  combos       <- list(c(0L, 0L), c(0L, 1L), c(1L, 0L), c(1L, 1L))
  joint_losses <- vapply(combos, function(co) {
    D_tmp <- D
    D_tmp[rownames(X_left),  "w"] <- co[1]
    D_tmp[rownames(X_right), "w"] <- co[2]
    global_objective_fn(D_tmp)
  }, numeric(1))
  best_idx <- which.min(joint_losses)
  new_loss <- joint_losses[best_idx]
  w_left   <- combos[[best_idx]][1]
  w_right  <- combos[[best_idx]][2]

  if (nearly_leq(new_loss, parent_loss)) {
    .log("[depth=%d] SPLIT: feature=%s, cut=%.4f | n_left=%d, n_right=%d | new_loss=%.6f, parent_loss=%.6f",
         depth, fj, cj, nrow(X_left), nrow(X_right), new_loss, parent_loss)
    .log("    joint losses={(0,0)=%.6f, (0,1)=%.6f, (1,0)=%.6f, (1,1)=%.6f}, best=(%d,%d)",
         joint_losses[1], joint_losses[2], joint_losses[3], joint_losses[4],
         w_left, w_right)

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

    combos       <- list(c(0L, 0L), c(0L, 1L), c(1L, 0L), c(1L, 1L))
    joint_losses <- vapply(combos, function(co) {
      D_tmp <- D
      D_tmp[rownames(X_left),  "w"] <- co[1]
      D_tmp[rownames(X_right), "w"] <- co[2]
      global_objective_fn(D_tmp)
    }, numeric(1))
    best_idx <- which.min(joint_losses)
    new_loss <- joint_losses[best_idx]
    w_left   <- combos[[best_idx]][1]
    w_right  <- combos[[best_idx]][2]

    if (nearly_leq(new_loss, parent_loss)) {
      .log("[depth=%d] SPLIT(after %d rejects): feature=%s, cut=%.4f | n_left=%d, n_right=%d | new_loss=%.6f, parent_loss=%.6f",
           depth, attempt, fj, cj, nrow(X_left), nrow(X_right), new_loss, parent_loss)
      .log("    joint losses={(0,0)=%.6f, (0,1)=%.6f, (1,0)=%.6f, (1,1)=%.6f}, best=(%d,%d)",
           joint_losses[1], joint_losses[2], joint_losses[3], joint_losses[4],
           w_left, w_right)

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
#' @keywords internal
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
#' @keywords internal
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
#' Calculates midpoint = (min + max)/2 using finite values only.
#' If there are no finite values, returns NA with a single warning.
#'
#' @param X A numeric vector.
#' @return Numeric scalar midpoint, or NA_real_ if no finite values.
#' @keywords internal
midpoint <- function(X) {
  if (!is.numeric(X)) {
    stop("`X` must be a numeric vector.", call. = FALSE)
  }
  if (length(X) == 0L) {
    warning("Empty vector provided to midpoint(). Returning NA.", call. = FALSE)
    return(NA_real_)
  }

  x_fin <- X[is.finite(X)]
  if (length(x_fin) == 0L) {
    warning("No finite values in `X`. Returning NA.", call. = FALSE)
    return(NA_real_)
  }

  rng <- range(x_fin)
  (rng[1] + rng[2]) / 2
}


#' Fit a shallow decision tree to characterize learned weights `w`
#'
#' @param X A data frame of covariates (features).
#' @param w A binary vector (0/1, TRUE/FALSE, or factor/character encoding 0/1)
#'   with exactly two classes present.
#' @param max_depth Integer, the maximum tree depth (default 3).
#' @return An `rpart` object representing the fitted decision tree.
#' @keywords internal
characterize_tree <- function(X, w, max_depth = 3) {

  # ---- Validate X ----
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame of covariates.", call. = FALSE)
  }
  check_no_na(X, colnames(X))

  # ---- Validate length(w) ----
  if (length(w) != nrow(X)) {
    stop("Length of `w` must equal the number of rows in `X`.", call. = FALSE)
  }

  # ---- Coerce w -> integer 0/1 with strict checks ----
  w0 <- w
  if (is.logical(w0)) w0 <- as.integer(w0)

  if (is.factor(w0)) w0 <- as.character(w0)
  if (is.character(w0)) {
    w_chr <- trimws(w0)
    ok <- w_chr %in% c("0", "1") | is.na(w_chr)
    if (!all(ok)) {
      stop("`w` must be binary with values 0/1 and contain exactly two classes.", call. = FALSE)
    }
    w0 <- suppressWarnings(as.integer(w_chr))
  }

  if (is.numeric(w0)) {
    bad <- !is.na(w0) & !(w0 %in% c(0, 1))
    if (any(bad)) {
      stop("`w` must be binary with values 0/1 and contain exactly two classes.", call. = FALSE)
    }
    w0 <- as.integer(w0)
  }

  if (!is.integer(w0)) {
    stop("`w` must be binary with values 0/1 and contain exactly two classes.", call. = FALSE)
  }

  # Do not allow NA in w for tree fitting
  if (anyNA(w0)) {
    stop("`w` must be binary with values 0/1 and contain exactly two classes.", call. = FALSE)
  }

  # ---- Must have exactly two classes present ----
  u <- sort(unique(w0))
  if (length(u) != 2L || !all(u %in% c(0L, 1L))) {
    stop("`w` must be binary with values 0/1 and contain exactly two classes.", call. = FALSE)
  }

  # ---- Fit rpart classifier ----
  w_factor <- factor(w0, levels = c(0, 1))
  df <- data.frame(w = w_factor, X)

  fit <- rpart::rpart(
    w ~ .,
    data = df,
    method = "class",
    parms = list(split = "gini"),
    control = rpart::rpart.control(
      maxdepth = max_depth, cp = 0.0, minsplit = 2, minbucket = 1
    ),
    model = TRUE
  )
  return(fit)
}
