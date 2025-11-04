#' Ensemble of weighted trees (loss/objective-agnostic) and Rashomon selection
#'
#' Builds multiple weighted trees, then identifies a "Rashomon set" of
#' top-performing trees and aggregates their weight assignments by majority vote.
#'
#' @param data A data frame containing the dataset (must include outcome, treatment, sample indicator).
#' @param outcome Name of the outcome column in `data`.
#' @param treatment Name of the treatment indicator column (0/1) in `data`.
#' @param sample Name of the sample indicator column (0/1) in `data`. Use NULL for single-sample SATE mode.
#' @param leaf_proba Probability mass for the "leaf" option in each tree (default 0.25).
#' @param seed Integer seed for reproducibility (default NULL).
#' @param num_trees Number of trees to grow in the forest (default 10).
#' @param vote_threshold Majority vote threshold in (0.5, 1] for final weight=1 (default 2/3).
#' @param explore_proba Probability of exploration at leaves in each tree (default 0.05).
#' @param feature_est "Ridge", "GBM", or a function(X, y, ...) returning a named, non-negative
#'   vector of importances; normalized to probabilities (default "Ridge").
#' @param feature_est_args Named list of extra args for a user-supplied `feature_est` function.
#' @param top_k_trees If TRUE, select top-k trees by objective; else use `cutoff` (default FALSE).
#' @param k Number of top trees if `top_k_trees = TRUE` (default 10).
#' @param cutoff If `top_k_trees = FALSE`, numeric cutoff or "baseline" (default "baseline").
#'   With "baseline", the cutoff is computed by evaluating `objective_fn` on the state with all `w=1`.
#' @param verbose If TRUE, prints 2 lines with (unweighted and weighted) estimate + SE. Default FALSE.
#' @param plot_tree If TRUE, plots the characterized tree (default TRUE). Guarded by `interactive()`.
#' @param plot_tree_args Named list forwarded to `rpart.plot::rpart.plot()`.
#' @param objective_fn Function \code{function(D) -> numeric} that scores the **entire state** (minimize).
#'   Default `objective_default()` reproduces prior behavior (SE proxy of PATE/TATE using `vsq` and `w`).
#' @param loss_fn Optional micro-evaluator \code{function(val, indices, D) -> numeric} that returns the
#'   objective after hypothetically setting `w=val` on `indices`. If `NULL`, ROOT wraps `objective_fn`
#'   via `loss_from_objective(objective_fn)`.
#'
#' @return S3 object of class "ROOT" with components:
#'   D_rash, D_forest, w_forest, rashomon_set, f, testing_data, tree_plot, estimate
#' @export
ROOT <- function(data,
                 outcome,
                 treatment,
                 sample,
                 leaf_proba = 0.25,
                 seed = NULL,
                 num_trees = 10,
                 vote_threshold = 2 / 3,
                 explore_proba = 0.05,
                 feature_est = "Ridge",
                 feature_est_args = list(),
                 top_k_trees = FALSE,
                 k = 10,
                 cutoff = "baseline",
                 verbose = FALSE,
                 objective_fn = objective_default,
                 loss_fn = NULL,
                 plot_tree = TRUE,
                 plot_tree_args = list(
                   type = 2, extra = 109, under = TRUE, faclen = 0, tweak = 1.1,
                   fallen.leaves = TRUE, box.palette = c("pink", "palegreen3"),
                   shadow.col = c("gray"), branch.lty = 3,
                   main = "Final Characterized Tree from Rashomon Set"
                 )
) {

  # Helpers
  coerce01 <- function(x, allow_na = TRUE) {
    if (is.numeric(x) || is.logical(x)) return(as.integer(x))
    if (is.factor(x)) x <- as.character(x)
    x_chr <- trimws(tolower(as.character(x)))
    map1 <- c("1","yes","treated","t","true")
    map0 <- c("0","no","control","c","false")
    out <- rep(NA_integer_, length(x_chr))
    out[x_chr %in% map1] <- 1L; out[x_chr %in% map0] <- 0L
    if (!allow_na && any(is.na(out))) stop("Non 0/1 values found.")
    out
  }
  .norm_feat_prob <- function(imp, X_df) {
    if (!is.numeric(imp) || is.null(names(imp))) {
      stop("Custom feature importance must be a *named* numeric vector.", call. = FALSE)
    }
    imp <- imp[colnames(X_df)]
    if (any(is.na(imp))) stop("Importance missing for some X_df columns.", call. = FALSE)
    if (any(imp < 0)) stop("Importances must be non-negative.", call. = FALSE)
    s <- sum(imp)
    if (!is.finite(s) || s <= 0) {
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- imp / s
    }
    feat_prob
  }

  if (!is.function(objective_fn)) stop("`objective_fn` must be a function(D)->numeric.")
  if (is.null(loss_fn)) loss_fn <- loss_from_objective(objective_fn)

  covariate_cols <- setdiff(names(data), c(outcome, treatment, sample))
  if (length(covariate_cols) == 0L) {
    stop("ROOT(): no covariate columns found (need at least one feature besides outcome/treatment/sample).",
         call. = FALSE)
  }

  # Basic tests
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.character(outcome)   || length(outcome)   != 1L || !(outcome   %in% names(data)))
    stop("`outcome` must be a single column name present in `data`.", call. = FALSE)
  if (!is.character(treatment) || length(treatment) != 1L || !(treatment %in% names(data)))
    stop("`treatment` must be a single column name present in `data`.", call. = FALSE)
  if (!is.null(sample)) {
    if (!is.character(sample) || length(sample) != 1L)
      stop("`sample` must be NULL or a single column name string.", call. = FALSE)
    if (!(sample %in% names(data)))
      stop("`sample` column not found; pass `sample = NULL` to run single-sample mode.", call. = FALSE)
  }
  if (!is.numeric(leaf_proba) || leaf_proba < 0 || leaf_proba > 1)
    stop("`leaf_proba` must be between 0 and 1.", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L))
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  if (!is.numeric(num_trees) || num_trees < 1) stop("`num_trees` must be positive.", call. = FALSE)
  num_trees <- as.integer(num_trees)
  if (!is.numeric(vote_threshold) || vote_threshold <= 0 || vote_threshold > 1)
    stop("`vote_threshold` must be in (0, 1].", call. = FALSE)
  if (!is.numeric(explore_proba) || explore_proba < 0 || explore_proba > 1)
    stop("`explore_proba` must be between 0 and 1.", call. = FALSE)
  if (!(is.character(feature_est) || is.function(feature_est)))
    stop("`feature_est` must be \"Ridge\", \"GBM\", or a function.", call. = FALSE)
  if (!is.logical(top_k_trees) || length(top_k_trees) != 1L)
    stop("`top_k_trees` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(k) || k < 1) stop("`k` must be a positive integer.", call. = FALSE)
  k <- as.integer(k)
  if (!(is.numeric(cutoff) || (is.character(cutoff) && cutoff == "baseline")))
    stop("`cutoff` must be \"baseline\" or numeric.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L)
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  if (!is.function(loss_fn)) stop("`loss_fn` must be a function (val, indices, D).", call. = FALSE)

  # Coerce treatment/sample to 0/1
  data[[treatment]] <- coerce01(data[[treatment]], allow_na = TRUE)
  if (!is.null(sample)) data[[sample]] <- coerce01(data[[sample]], allow_na = FALSE)
  Tr_check <- data[[treatment]]
  if (!all(Tr_check[!is.na(Tr_check)] %in% c(0L,1L)))
    stop("`treatment` column must be 0/1 (after coercion).", call. = FALSE)
  if (!is.null(sample)) {
    S_check <- data[[sample]]
    if (!all(S_check %in% c(0L,1L)))
      stop("`sample` column must be 0/1.", call. = FALSE)
  }

  # Single vs two-sample mode
  single_sample_mode <- is.null(sample) ||
    (all(data[[sample]] %in% c(0,1)) && (all(data[[sample]] == 1) || all(data[[sample]] == 0)))
  if (!single_sample_mode && (all(data[[sample]] == 0) || all(data[[sample]] == 1)))
    stop("Sample indicator `", sample, "` has no variation (all 0 or all 1).", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  # Estimate orthogonal scores (v)
  if (single_sample_mode) {
    dml <- estimate_dml_single(data, outcome = outcome, treatment = treatment, crossfit = 5)
    df_v <- dml$df_v
    testing_data <- dml$data2
    S_vec <- rep(1L, nrow(testing_data))      # all target in single-sample path
    lX_col <- NA_real_                        # no overlap factor
  } else {
    dml <- estimate_dml(data, outcome = outcome, treatment = treatment, sample = sample, crossfit = 5)
    df_v <- dml$df_v
    testing_data <- dml$data2
    S_vec <- testing_data[[sample]]
    lX_col <- 1 / df_v$b                      # overlap factor = 1/b
  }

  Y_vec <- testing_data[[outcome]]
  Tr_vec <- testing_data[[treatment]]
  drop_cols <- c(outcome, treatment, if (!is.null(sample) && (sample %in% names(testing_data))) sample)
  X_df <- testing_data[, setdiff(names(testing_data), drop_cols), drop = FALSE]

  n <- nrow(testing_data)
  v_vals <- df_v$te
  vsq_vals <- df_v$te_sq

  # Feature probabilities for splitting
  no_cov <- (ncol(X_df) == 0)
  features <- c("leaf", colnames(X_df))

  compute_split_feature <- function() {
    if (ncol(X_df) < 1L) {
      stop("compute_split_feature(): no features available for weighting.", call. = FALSE)
    }
    if (is.function(feature_est)) {
      imp <- if (is.null(seed)) {
        do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
      } else {
        withr::with_seed(seed + 777L,
                         do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
        )
      }
      feat_prob <- .norm_feat_prob(imp, X_df)
    } else if (tolower(feature_est) == "ridge") {
      if (ncol(X_df) == 0 || all(abs(vsq_vals) < .Machine$double.eps) || !is.finite(sum(abs(vsq_vals)))) {
        feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
        if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
      } else {
        ridge_fit <- MASS::lm.ridge(vsq_vals ~ ., data = data.frame(X_df), lambda = 1)
        beta <- as.numeric(ridge_fit$coef); names(beta) <- colnames(X_df)
        if (all(beta == 0) || !is.finite(sum(abs(beta)))) {
          feat_prob <- rep(1 / ncol(X_df), ncol(X_df)); names(feat_prob) <- colnames(X_df)
        } else {
          feat_prob <- abs(beta) / sum(abs(beta))
        }
      }
    } else { # GBM
      df_gbm <- data.frame(X_df, vsq = vsq_vals)
      fit_gbm <- function() {
        gbm::gbm(vsq ~ ., data = df_gbm, distribution = "gaussian",
                 n.trees = 100, interaction.depth = 3, shrinkage = 0.1,
                 bag.fraction = 0.5, train.fraction = 1.0, n.minobsinnode = 2, verbose = FALSE)
      }
      gbm_fit <- if (is.null(seed)) fit_gbm() else withr::with_seed(seed + 1000L, fit_gbm())
      rel_inf <- gbm::relative.influence(gbm_fit, n.trees = 100, sort. = FALSE)
      names(rel_inf) <- colnames(X_df)
      if (sum(rel_inf) <= 0 || !is.finite(sum(rel_inf))) {
        feat_prob <- rep(1 / ncol(X_df), ncol(X_df)); names(feat_prob) <- colnames(X_df)
      } else {
        feat_prob <- rel_inf / sum(rel_inf)
      }
    }
    proba <- c(leaf_proba, feat_prob); proba <- proba / sum(proba)
    stats::setNames(proba, features)
  }

  split_feature <- if (no_cov) stats::setNames(1, "leaf") else compute_split_feature()

  # Build forest
  D_forest <- data.frame(X_df, v = v_vals, vsq = vsq_vals, S = S_vec, stringsAsFactors = FALSE)
  D_forest$lX <- lX_col

  w_forest <- vector("list", num_trees)
  for (t_idx in seq_len(num_trees)) {
    D_tree <- data.frame(X_df, v = v_vals, vsq = vsq_vals, w = rep(1, n), S = S_vec, stringsAsFactors = FALSE)
    rownames(D_tree) <- as.character(seq_len(n)); rownames(X_df) <- as.character(seq_len(n))
    grow_one <- function() {
      split_node(
        split_feature = split_feature,
        X = D_tree, D = D_tree, parent_loss = Inf, depth = 0L,
        explore_proba = explore_proba,
        choose_feature_fn = choose_feature,
        loss_fn = loss_fn,
        reduce_weight_fn = reduce_weight,
        objective_fn = objective_fn
      )
    }
    tree_res <- if (is.null(seed)) grow_one() else withr::with_seed(seed + t_idx, grow_one())
    D_updated <- tree_res$D
    D_forest[[paste0("w_tree_", t_idx)]] <- D_updated$w
    w_forest[[t_idx]] <- tree_res
  }

  # Rashomon set selection
  obj_values <- vapply(w_forest, function(res) res[["local objective"]], numeric(1))
  tree_indices <- seq_along(w_forest)
  if (top_k_trees) {
    if (k > num_trees) { warning("k > num_trees; using k = num_trees.", call. = FALSE); k <- num_trees }
    ord <- order(obj_values, na.last = NA)
    rashomon_set <- utils::head(ord, k)
  } else {
    cutoff_val <- if (identical(cutoff, "baseline")) {
      D_base <- D_forest
      D_base$w <- 1
      objective_fn(D_base)
    } else as.numeric(cutoff)
    if (!is.finite(cutoff_val)) cutoff_val <- Inf
    rashomon_set <- tree_indices[obj_values < cutoff_val]
  }
  if (length(rashomon_set) == 0) warning("No trees selected into Rashomon set.", call. = FALSE)
  not_in_set <- setdiff(tree_indices, rashomon_set)

  # Keep selected trees and votes
  D_rash <- D_forest
  if (length(not_in_set) > 0) {
    drop_cols <- paste0("w_tree_", not_in_set)
    keep_cols <- setdiff(names(D_rash), drop_cols)
    D_rash <- D_rash[, keep_cols, drop = FALSE]
  }
  weight_cols <- grep("^w_tree_", names(D_rash), value = TRUE)
  D_weights <- if (length(weight_cols) > 0) D_rash[, weight_cols, drop = FALSE] else data.frame()
  if (ncol(D_weights) > 0) {
    row_means <- rowMeans(D_weights)
    D_rash$w_opt <- as.integer(row_means > vote_threshold)
    D_rash$vote_count <- rowSums(D_weights)
  } else {
    D_rash$w_opt <- integer(nrow(D_rash))
    D_rash$vote_count <- integer(nrow(D_rash))
  }

  # Final summarized tree
  final_classifier <- if (no_cov) NULL else characterize_tree(X_df, as.factor(D_rash$w_opt))

  # Optional plot (unguarded)
  tree_plot <- NULL
  if (isTRUE(plot_tree) && !is.null(final_classifier)) {
    if (!requireNamespace("rpart.plot", quietly = TRUE)) {
      warning("Package 'rpart.plot' is not installed; skipping tree plot.", call. = FALSE)
    } else {
      # draw and capture base graphics
      do.call(rpart.plot::prp, c(list(final_classifier), plot_tree_args))
      tree_plot <- grDevices::recordPlot()

      # optionally print immediately in console/knit:
      grDevices::replayPlot(tree_plot)
    }
  }

  # Estimands: unweighted vs weighted (with SEs only)
  in_S  <- if (single_sample_mode) rep(TRUE, n) else (S_vec == 1L)
  v_sel <- v_vals[in_S]
  w_sel <- if ("w_opt" %in% names(D_rash)) D_rash$w_opt[in_S] else rep(1L, length(v_sel))

  ## Unweighted (ATE in RCT or TATE)
  est_label_unw <- if (single_sample_mode) "SATE" else "TATE"
  mu_unw <- mean(v_sel, na.rm = TRUE)
  n_eff_unw <- sum(!is.na(v_sel))
  # SE = sqrt( var / n ), where var() uses Bessel correction (n-1)
  se_unw <- if (n_eff_unw > 1) sqrt(stats::var(v_sel, na.rm = TRUE) / n_eff_unw) else NA_real_

  ## Weighted (WATE in RCT or WTATE) — binary weights assumed
  est_label_w <- if (single_sample_mode) "WATE" else "WTATE"
  den_w <- sum(w_sel, na.rm = TRUE)  # effective n when w ∈ {0,1}
  if (isTRUE(den_w > 0)) {
    mu_w <- sum(w_sel * v_sel, na.rm = TRUE) / den_w
    # SE of weighted mean (binary weights): sqrt( sum w*(v - mu)^2 / den_w^2 )
    se_w <- sqrt( sum(w_sel * (v_sel - mu_w)^2, na.rm = TRUE) / (den_w^2) )
  } else {
    mu_w <- NA_real_; se_w <- NA_real_
  }

  if (verbose) {
    message(sprintf("%s (unweighted) = %.6f, SE = %.6f", est_label_unw, mu_unw, se_unw))
    message(sprintf("%s (weighted)   = %.6f, SE = %.6f", est_label_w,   mu_w,   se_w))
  }

  # Assemble result (SEs only; SDs removed)
  results <- list(
    D_rash = D_rash,
    D_forest = D_forest,
    w_forest = w_forest,
    rashomon_set = rashomon_set,
    f = final_classifier,
    testing_data = testing_data,
    tree_plot = tree_plot,
    estimate = list(
      estimand_unweighted = est_label_unw,
      value_unweighted    = mu_unw,
      se_unweighted       = se_unw,
      estimand_weighted   = est_label_w,
      value_weighted      = mu_w,
      se_weighted         = se_w,
      n_analysis          = sum(in_S),
      sum_w               = den_w
    )
  )
  results$single_sample_mode <- single_sample_mode
  class(results) <- c("ROOT", "list")
  return(results)
}
