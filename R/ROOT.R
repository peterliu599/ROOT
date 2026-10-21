#' Rashomon Set of Optimal Trees (ROOT) for Functional Optimization
#'
#' @description
#' ROOT (Rashomon Set of Optimal Trees) is a general-purpose functional
#' optimization algorithm that learns interpretable, tree-structured binary
#' weight functions \eqn{w(X) \in \{0, 1\}}. Given a dataset \eqn{D_n} and a
#' global objective function \eqn{L(w, D_n)}, ROOT searches over the space of
#' decision trees to find weight assignments that minimize the objective function.
#'
#' @section The optimization problem:
#' ROOT solves the functional optimization problem:
#'
#' \deqn{w^* \in \arg\min_w L(D_n, w)}
#'
#' where \eqn{w: \mathbb{R}^p \to \{0, 1\}} maps a \eqn{p}-dimensional
#' covariate vector to a binary include/exclude decision. The key challenge is
#' that, unlike standard tree algorithms, the global loss
#' \eqn{L(D_n, w)} is not decomposable as a sum of losses over
#' independent subsets of the data. This means conventional greedy,
#' divide-and-conquer tree-building strategies do not apply. ROOT addresses
#' this through a randomization-based tree construction with an
#' explore-exploit strategy.
#'
#' @section How ROOT works:
#' The algorithm proceeds in several stages:
#'
#' \enumerate{
#'   \item \strong{Feature importance estimation:} Split probabilities are
#'     estimated using Ridge regression, Gradient Boosting Machine (GBM),
#'     or a user-supplied function, biasing the search toward covariates likely
#'     to be informative.
#'   \item \strong{Stochastic tree construction:} \code{num_trees} trees are
#'     grown. At each internal node, a feature is drawn according to the
#'     estimated split probabilities (or a "leaf" token is drawn, terminating
#'     the branch). Splits are made at the midpoint of the selected feature's
#'     empirical distribution. An explore-exploit strategy assigns leaf
#'     weights: with probability \code{explore_proba} a random weight is
#'     chosen; otherwise the greedy optimal weight (reducing the global
#'     objective) is used.
#'   \item \strong{Rashomon set selection:} Trees are ranked by their global
#'     objective values. The top-\code{k} trees (or all trees below a cutoff)
#'     form the Rashomon set: a collection of near-optimal but
#'     potentially different models, each providing a characterization
#'     of the optimal weight function.
#'   \item \strong{Aggregation:} Per-observation votes from the Rashomon set
#'     are combined (by default, majority vote) to produce the final weight
#'     vector \code{w_opt}.
#'   \item \strong{Characteristic tree:} A single summary decision tree is
#'     fitted to the aggregated \code{w_opt} assignments, providing a concise,
#'     interpretable description of the weight function.
#' }
#'
#' @section Generalizability mode:
#' When \code{generalizability_path = TRUE}, ROOT implements the methodology
#' of Parikh et al. (2025) for characterizing underrepresented subgroups in
#' trial-to-target generalizability analyses. In this mode:
#'
#' \itemize{
#'   \item \code{data} must contain columns \code{Y} (outcome), \code{Tr}
#'     (treatment, 0/1), and \code{S} (sample indicator, 1 = trial, 0 =
#'     target).
#'   \item ROOT internally computes transportability scores based on
#'     inverse-probability weighting (IPW), estimates the selection model
#'     \eqn{P(S = 1 \mid X)}, and constructs inverse propensity weighting
#'     (IPW) transported scores.
#'   \item The default objective minimizes the variance of the weighted
#'     target average treatment effect (WTATE) estimator. This objective
#'     accounts for both the selection odds (trial participation
#'     probability) and treatment effect heterogeneity, so that subgroups
#'     are flagged as underrepresented primarily when they both lack trial
#'     representation and exhibit effect modification.
#'   \item The output includes the unweighted sample target average
#'     treatment effect (STATE) and the WTATE with standard errors.
#' }
#'
#' See \code{\link{characterizing_underrep}} for a higher-level wrapper that
#' additionally produces a leaf-level summary table, and
#' \code{vignette("generalizability_path_example")} for a worked example.
#'
#' @section General optimization mode:
#' When \code{generalizability_path = FALSE}, ROOT operates as a general
#' functional optimizer. The user supplies any \code{data.frame} and
#' (optionally) a custom \code{global_objective_fn}. If no objective is
#' supplied, ROOT uses the default variance-based loss
#' \deqn{L(w) = \sqrt{\sum_i w_i\, vsq_i / (\sum_i w_i)^2},}
#' where \code{vsq} is taken from the data as follows:
#' \itemize{
#'   \item if \code{data} has a \code{vsq} column (per-unit variance
#'     contribution), it is used as supplied;
#'   \item otherwise, if \code{data} has a \code{v} column (per-unit score),
#'     \eqn{vsq_i = (v_i - \bar v_w)^2} is recomputed for every candidate
#'     rule, where \eqn{\bar v_w} is the mean of \code{v} over the
#'     currently included units;
#'   \item otherwise, \eqn{vsq_i = 1} is used as a fallback. The loss then
#'     equals \eqn{1/\sqrt{n_{kept}}}, which is minimized by keeping every
#'     observation, so a message is issued and users should supply
#'     \code{vsq}, \code{v}, or a \code{global_objective_fn}.
#' }
#' See \code{vignette("optimization_path_example")} for an example with a
#' user-supplied, covariance-based objective.
#'
#'
#' @param data A data.frame containing the dataset.
#'
#'   In \strong{general optimization mode} (\code{generalizability_path =
#'   FALSE}), \code{data} can contain any covariates and auxiliary columns.
#'   The user either supplies a \code{global_objective_fn} that takes a data
#'   frame with a column \code{w} and returns a scalar loss, or relies on the
#'   default objective, which reads an optional \code{vsq} or \code{v}
#'   column (see the General optimization mode section). Columns named
#'   \code{vsq} and \code{v} are auxiliary and are never used for splitting.
#'
#'   In \strong{generalizability mode} (\code{generalizability_path = TRUE}),
#'   \code{data} must contain columns \code{"Y"} (outcome), \code{"Tr"}
#'   (treatment indicator, 0/1), and \code{"S"} (sample indicator, 1 = trial,
#'   0 = target). ROOT internally constructs transportability scores and, if
#'   no custom objective is given, uses a default variance-based loss.
#' @param global_objective_fn A function with signature
#'   \code{function(D) -> numeric} scoring the entire state and minimized by
#'   ROOT. \code{D} is the working data frame containing at least a column
#'   \code{w} with the current weight assignments. If \code{NULL} (default), a
#'   variance-based objective is used.
#' @param generalizability_path \code{Logical(1)}. If \code{TRUE}, use the
#'   built-in transportability objective for trial-to-target generalizability.
#'   If \code{FALSE}, treat \code{data} as arbitrary and rely on
#'   \code{global_objective_fn}. Default \code{FALSE}.
#' @param leaf_proba A numeric tuning parameter that increases the chance
#'   a node stops splitting by selecting a synthetic \code{"leaf"} feature.
#'   Internally, the probability of choosing \code{"leaf"} is
#'   \code{leaf_proba / (1 + leaf_proba)} (assuming the
#'   covariate probabilities sum to 1). Higher values produce shallower,
#'   sparser trees. Default \code{0.25}.
#' @param seed An optional numeric seed for reproducibility.
#' @param num_trees An integer number of trees to grow. More trees explore the
#'   solution space more thoroughly. Default 10.
#' @param vote_threshold Controls how per-observation votes from the Rashomon
#'   set trees are aggregated into the final binary weight \code{w_opt}.
#'   Accepts one of:
#'   \itemize{
#'     \item A \strong{numeric} in \code{(0, 1)}: \code{w_opt = 1} when the
#'       mean fraction of Rashomon-set trees voting 1 strictly exceeds this
#'       value. An error is raised if any observation's mean vote equals the
#'       threshold exactly (tie is undefined). Default \code{2/3}.
#'     \item \code{"majority"}: equivalent to the numeric threshold \code{2/3}.
#'     \item \code{"mean"}: \code{w_opt = 1} when the mean vote strictly
#'       exceeds \code{0.5}. Errors on exact \code{0.5} ties.
#'     \item \code{"median"}: \code{w_opt = 1} when the per-row median vote
#'       strictly exceeds \code{0.5}. Errors on exact \code{0.5} ties.
#'     \item A \strong{function} with signature
#'       \code{function(votes) -> integer vector}, where \code{votes} is a
#'       numeric matrix with one row per observation and one column per
#'       Rashomon-set tree (each cell 0 or 1). Must return an integer vector
#'       of 0s and 1s of length \code{nrow(votes)}.
#'   }
#' @param explore_proba A numeric giving the exploration probability at
#'   leaves. With probability \code{explore_proba}, a leaf receives a random
#'   weight; otherwise the greedy-optimal weight is used. Default \code{0.05}.
#' @param feature_est Either \code{"Ridge"}, \code{"GBM"}, or a custom
#'   \code{function(X, y, ...) -> named numeric vector} returning nonnegative
#'   importance values with names matching columns of \code{X}. Used to bias
#'   which covariates are chosen for splitting. If the method fails, ROOT
#'   falls back to uniform feature sampling with a warning.
#' @param feature_est_args A list of additional arguments passed to a
#'   user-supplied \code{feature_est} function.
#' @param top_k_trees \code{Logical(1)}. If \code{TRUE}, select the top
#'   \code{k} trees by objective value for the Rashomon set. If \code{FALSE},
#'   use the \code{cutoff} threshold. Default \code{FALSE}.
#' @param k An integer giving the number of top trees when
#'   \code{top_k_trees = TRUE}. Default 10.
#' @param cutoff A numeric or \code{"baseline"}. Used as the Rashomon cutoff
#'   when \code{top_k_trees = FALSE}. \code{"baseline"} uses the objective at
#'   \eqn{w \equiv 1} (all weights equal to 1, i.e., no exclusions).
#' @param max_depth Maximum depth of each tree grown during the forest
#'   construction stage. A node at \code{depth == max_depth} is forced to be a
#'   leaf. Shallower trees are more interpretable but less flexible. Default
#'   \code{8}.
#' @param min_leaf_n Minimum number of observations required in a node for
#'   splitting to be attempted. If a node contains fewer than
#'   \code{min_leaf_n} observations it becomes a leaf. Default \code{2}.
#' @param max_rejects_per_node Maximum number of consecutive rejected splits
#'   (splits that do not improve the objective) allowed at a single node
#'   before the node is forced to become a leaf. This prevents infinite
#'   recursion in pathological cases. Default \code{10}.
#' @param verbose \code{Logical(1)}. If \code{TRUE}, prints the unweighted and
#'   (when available) weighted treatment effect estimates and standard errors
#'   in generalizability mode.
#' @param positivity_trim Controls pre-trimming of target observations with
#'   near-zero trial selection probability \eqn{P(S = 1 \mid X)}, addressing
#'   structural positivity violations before ROOT optimization. Only used when
#'   \code{generalizability_path = TRUE}; ignored otherwise.
#'   \itemize{
#'     \item \code{NULL} (default): no positivity trimming.
#'     \item A \strong{numeric} in \eqn{(0, 0.5)}: target observations with
#'       \eqn{\hat P(S = 1 \mid X)} below this value are removed.
#'     \item \code{"crump"}: use the data-driven threshold from Crump, Hotz,
#'       Imbens & Mitnik (2009, \emph{Biometrika}), which solves for the
#'       optimal symmetric trimming cutoff \eqn{\alpha} that minimizes the
#'       asymptotic variance of the treatment-effect estimator.
#'   }
#'   Note that trimming changes the estimand: the WTATE is now estimated over
#'   the subpopulation \eqn{\{x : \hat P(S = 1 \mid X = x) \ge \alpha\}}.
#' @param split_strategy Character; either \code{"midpoint"} (default) or
#'   \code{"random_quantile"}. Controls how the split point for a chosen
#'   feature is determined at each tree node. \code{"midpoint"} uses the
#'   classic \code{(min + max) / 2} rule; \code{"random_quantile"} draws a
#'   uniform random quantile within \code{quantile_range}, introducing
#'   controlled stochasticity that improves covariate-space coverage across
#'   the forest.
#' @param quantile_range Numeric vector of length 2 giving the lower and upper
#'   bounds of the quantile range when
#'   \code{split_strategy = "random_quantile"}. Default \code{c(0.1, 0.9)}.
#'   Ignored when \code{split_strategy = "midpoint"}.
#' @param min_kept Numeric in \eqn{[0, 1)}. Generalizability mode only:
#'   the smallest share of analysis (trial) units that a candidate rule may
#'   keep when the default objective is used. Rules keeping fewer than
#'   \code{max(2, ceiling(min_kept * n))} units receive an infinite objective
#'   value. Default \code{0.1}. Without this bound, the finite-sample standard
#'   error can be made arbitrarily small by keeping a handful of units whose
#'   scores happen to be similar, a degenerate solution that a thorough search
#'   will find. Set \code{min_kept = 0} to remove the bound. Ignored when a
#'   \code{global_objective_fn} is supplied or when
#'   \code{generalizability_path = FALSE}.
#' @param summary_max_depth Positive integer. Maximum depth of the single
#'   characteristic (summary) tree fitted to the final rule \code{w_opt} with
#'   \code{rpart}. Default \code{3}. A rule that excludes \eqn{m} separate
#'   subgroups may need depth \eqn{m} to be described exactly; larger values
#'   describe the ensemble rule more faithfully at the cost of a less compact
#'   tree. This affects only the summary tree, not \code{w_opt} or the
#'   reported estimates.
#'
#' @return An object of class \code{"ROOT"} (a list) with elements:
#'   \item{D_rash}{Data frame containing the Rashomon-set tree votes and the
#'     final aggregated weight \code{w_opt} for each observation.}
#'   \item{D_forest}{Data frame with all forest-level working columns.}
#'   \item{w_forest}{List of per-tree results from the tree-building routine.}
#'   \item{rashomon_set}{Integer vector of indices identifying which trees
#'     were selected into the Rashomon set.}
#'   \item{global_objective_fn}{The objective function used.}
#'   \item{f}{The characteristic (summary) tree fitted to \code{w_opt}, as an
#'     \code{rpart} object. \code{NULL} if all observations received the same
#'     weight or no covariates were available.}
#'   \item{testing_data}{Data frame of observations used for optimization
#'     (trial units when \code{generalizability_path = TRUE}).}
#'   \item{estimate}{(Only when \code{generalizability_path = TRUE}) A list
#'     with the unweighted (STATE) and weighted (WTATE) point estimates and
#'     their standard errors, the number of analysis observations
#'     (\code{n_analysis}), the size of the kept set (\code{n_A}), the
#'     selection-odds re-normalization factor \code{kappa_w}, and
#'     \code{se_weighted_note}.}
#'   \item{generalizability_path}{Logical flag echoing the input.}
#'   \item{positivity_trim_info}{(Only when \code{positivity_trim} is not
#'     \code{NULL}) A list with the threshold used (\code{alpha}), the number
#'     of target observations removed (\code{n_trimmed}), and the method
#'     (\code{method}).}
#'
#' @references
#' Parikh H, Ross RK, Stuart E, Rudolph KE (2025).
#' "Who Are We Missing?: A Principled Approach to Characterizing the
#' Underrepresented Population."
#' \emph{Journal of the American Statistical Association}.
#' \doi{10.1080/01621459.2025.2495319}
#'
#' Crump RK, Hotz VJ, Imbens GW, Mitnik OA (2009).
#' "Dealing with Limited Overlap in Estimation of Average Treatment Effects."
#' \emph{Biometrika}, \strong{96}(1), 187--199.
#'
#' @seealso
#' \code{\link{characterizing_underrep}} for a higher-level wrapper with
#' leaf-summary output; \code{vignette("generalizability_path_example")} for
#' the generalizability workflow;
#' \code{vignette("optimization_path_example")} for general optimization.
#'
#' @examples
#' # --- Generalizability mode ---
#' data(diabetes_data, package = "ROOT")
#' root_fit <- ROOT(
#'   data                  = diabetes_data,
#'   generalizability_path = TRUE,
#'   num_trees             = 20,
#'   top_k_trees           = TRUE,
#'   k                     = 10,
#'   seed                  = 123
#' )
#' # --- General optimization mode (custom, covariance-based objective) ---
#' set.seed(123)
#' n_assets <- 100
#'
#' # Asset features
#' volatility <- runif(n_assets, 0.05, 0.40)  # annualised volatility
#' beta       <- runif(n_assets, 0.5,  1.8)   # market beta
#' sector     <- sample(c("Tech", "Finance", "Energy", "Health"),
#'                      n_assets, replace = TRUE)
#'
#' # Simulate returns: r_i = beta_i * r_market + epsilon_i
#' market      <- rnorm(1000, 0.0005, 0.01)
#' returns_mat <- sapply(seq_len(n_assets), function(i)
#'   beta[i] * market + rnorm(1000, 0, volatility[i] / sqrt(252))
#' )
#' Sigma_hat <- stats::cov(returns_mat)   # 100 x 100 return covariance
#'
#' # Variance of the equal-weight portfolio of the included assets:
#' # L(w) = w' Sigma w / (1'w)^2. It depends on all pairwise covariances,
#' # so it cannot be split into independent per-asset contributions.
#' portfolio_variance <- function(D) {
#'   w <- D$w
#'   if (sum(w) == 0) return(Inf)
#'   sum(Sigma_hat[w == 1, w == 1]) / sum(w)^2  # = w' Sigma w / (1'w)^2
#' }
#'
#' my_data <- data.frame(
#'   vsq    = diag(Sigma_hat),   # marginal variances; guide split sampling only
#'   vol    = volatility,
#'   beta   = beta,
#'   sector = as.integer(factor(sector))
#' )
#'
#' opt_fit <- ROOT(
#'   data                = my_data,
#'   global_objective_fn = portfolio_variance,
#'   num_trees           = 20,
#'   top_k_trees         = TRUE,
#'   k                   = 10,
#'   seed                = 42
#' )
#'
#' @export
ROOT <- function(data,
                 global_objective_fn   = NULL,
                 generalizability_path = FALSE,
                 leaf_proba            = 0.25,
                 seed                  = NULL,
                 num_trees             = 10,
                 vote_threshold        = 2 / 3,
                 explore_proba         = 0.05,
                 feature_est           = "Ridge",
                 feature_est_args      = list(),
                 top_k_trees           = FALSE,
                 k                     = 10,
                 cutoff                = "baseline",
                 max_depth             = 8L,
                 min_leaf_n            = 2L,
                 max_rejects_per_node  = 10L,
                 verbose               = FALSE,
                 positivity_trim       = NULL,
                 split_strategy        = "midpoint",
                 quantile_range        = c(0.1, 0.9),
                 summary_max_depth     = 3L,
                 min_kept              = 0.1) {

  # ---- Small helpers ------------------------------------------------------
  coerce01 <- function(x, allow_na = TRUE) {
    # Logical -> {0,1}
    if (is.logical(x)) x <- as.integer(x)

    # Numeric: must be exactly 0/1 (optionally NA)
    if (is.numeric(x)) {
      bad <- !is.na(x) & !(x %in% c(0, 1))
      if (any(bad)) stop("Non 0/1 values found.", call. = FALSE)
      out <- as.integer(x)
      if (!allow_na && any(is.na(out))) stop("Non 0/1 values found.", call. = FALSE)
      return(out)
    }

    # Character/factor mapping (existing behavior)
    if (is.factor(x)) x <- as.character(x)
    x_chr <- trimws(tolower(as.character(x)))
    map1 <- c("1", "yes", "treated", "t", "true")
    map0 <- c("0", "no", "control", "c", "false")
    out <- rep(NA_integer_, length(x_chr))
    out[x_chr %in% map1] <- 1L
    out[x_chr %in% map0] <- 0L
    if (!allow_na && any(is.na(out)))
      stop("Non 0/1 values found.", call. = FALSE)
    out
  }

  .norm_feat_prob <- function(imp, X_df) {
    if (!is.numeric(imp) || is.null(names(imp))) {
      stop("Custom feature importance must be a named numeric vector.", call. = FALSE)
    }
    imp <- imp[colnames(X_df)]
    if (any(is.na(imp))) stop("Importance missing for some X_df columns.", call. = FALSE)
    if (any(imp < 0))   stop("Importances must be non-negative.", call. = FALSE)
    s <- sum(imp)
    if (!is.finite(s) || s <= 0) {
      warning(
        "Feature importance was degenerate; now using a simple random probabilistic ",
        "setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- imp / s
    }
    feat_prob
  }

  compute_split_feature <- function(X_df, vsq_vals) {
    if (ncol(X_df) < 1L) {
      stop("compute_split_feature(): no features available for weighting.", call. = FALSE)
    }

    ## 1) User-supplied feature importance function
    if (is.function(feature_est)) {
      imp <- try(
        if (is.null(seed)) {
          do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
        } else {
          withr::with_seed(
            seed + 777L,
            do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
          )
        },
        silent = TRUE
      )
      if (inherits(imp, "try-error")) {
        warning(
          "feature_est() failed; now using a simple random probabilistic setting ",
          "which is uniform feature sampling.",
          call. = FALSE
        )
        feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
        if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
        return(feat_prob)
      }
      return(.norm_feat_prob(imp, X_df))
    }

    ## 2) Ridge-based importance
    if (tolower(feature_est) == "ridge") {
      use_ridge <- ncol(X_df) > 0 &&
        any(abs(vsq_vals) > .Machine$double.eps) &&
        is.finite(sum(abs(vsq_vals)))
      if (use_ridge) {
        ridge_fit <- try(
          MASS::lm.ridge(vsq_vals ~ ., data = data.frame(X_df), lambda = 1),
          silent = TRUE
        )
        if (!inherits(ridge_fit, "try-error")) {
          beta <- as.numeric(ridge_fit$coef)
          names(beta) <- colnames(X_df)
          if (!all(beta == 0) && is.finite(sum(abs(beta)))) {
            feat_prob <- abs(beta) / sum(abs(beta))
            return(feat_prob)
          }
        }
      }
      warning(
        "Ridge feature importance failed or was degenerate; now using a simple random ",
        "probabilistic setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
      return(feat_prob)
    }

    ## 3) GBM-based importance
    df_gbm <- data.frame(X_df, vsq = vsq_vals)
    fit_gbm <- function() {
      gbm::gbm(
        vsq ~ ., data = df_gbm, distribution = "gaussian",
        n.trees = 100, interaction.depth = 3, shrinkage = 0.1,
        bag.fraction = 0.5, train.fraction = 1.0,
        n.minobsinnode = 2, verbose = FALSE
      )
    }
    gbm_fit <- try(
      if (is.null(seed)) fit_gbm() else withr::with_seed(seed + 1000L, fit_gbm()),
      silent = TRUE
    )
    if (inherits(gbm_fit, "try-error")) {
      warning(
        "GBM feature importance failed; now using a simple random probabilistic ",
        "setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
      return(feat_prob)
    }
    rel_inf <- gbm::relative.influence(gbm_fit, n.trees = 100, sort. = FALSE)
    names(rel_inf) <- colnames(X_df)
    if (sum(rel_inf) <= 0 || !is.finite(sum(rel_inf))) {
      warning(
        "GBM feature importance was degenerate; now using a simple random probabilistic ",
        "setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- rel_inf / sum(rel_inf)
    }
    feat_prob
  }

  # ---- Argument checks ----------------------------------------------------
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.logical(generalizability_path) || length(generalizability_path) != 1L)
    stop("`generalizability_path` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(leaf_proba) || leaf_proba < 0 || leaf_proba > 1)
    stop("`leaf_proba` must be between 0 and 1.", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L))
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  if (!is.numeric(num_trees) || num_trees < 1)
    stop("`num_trees` must be positive.", call. = FALSE)
  num_trees <- as.integer(num_trees)
  # Validate vote_threshold
  .valid_str <- c("majority", "mean", "median")
  if (is.numeric(vote_threshold)) {
    if (length(vote_threshold) != 1L || vote_threshold <= 0 || vote_threshold > 1)
      stop("`vote_threshold` as a numeric must be a single value in (0, 1].", call. = FALSE)
  } else if (is.character(vote_threshold)) {
    if (length(vote_threshold) != 1L || !vote_threshold %in% .valid_str)
      stop(sprintf(
        '`vote_threshold` as a string must be one of: %s.',
        paste(sprintf('"%s"', .valid_str), collapse = ", ")
      ), call. = FALSE)
  } else if (!is.function(vote_threshold)) {
    stop(
      '`vote_threshold` must be a numeric threshold in (0, 1], ',
      'one of "majority" / "mean" / "median", ',
      "or a function(votes) -> integer vector.",
      call. = FALSE
    )
  }
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
  if (!is.numeric(max_depth) || length(max_depth) != 1L || max_depth < 0)
    stop("`max_depth` must be a non-negative integer.", call. = FALSE)
  max_depth <- as.integer(max_depth)
  if (!is.numeric(min_leaf_n) || length(min_leaf_n) != 1L || min_leaf_n < 1)
    stop("`min_leaf_n` must be a positive integer.", call. = FALSE)
  min_leaf_n <- as.integer(min_leaf_n)
  if (!is.numeric(max_rejects_per_node) || length(max_rejects_per_node) != 1L ||
      max_rejects_per_node < 1)
    stop("`max_rejects_per_node` must be a positive integer.", call. = FALSE)
  max_rejects_per_node <- as.integer(max_rejects_per_node)
  # Validate positivity_trim
  if (!is.null(positivity_trim)) {
    if (!(identical(positivity_trim, "crump") ||
          (is.numeric(positivity_trim) && length(positivity_trim) == 1L &&
           positivity_trim > 0 && positivity_trim < 0.5))) {
      stop('`positivity_trim` must be NULL, "crump", or a numeric in (0, 0.5).',
           call. = FALSE)
    }
  }
  # Validate split_strategy
  if (!is.character(split_strategy) || length(split_strategy) != 1L ||
      !split_strategy %in% c("midpoint", "random_quantile")) {
    stop('`split_strategy` must be "midpoint" or "random_quantile".', call. = FALSE)
  }
  # Validate quantile_range
  if (!is.numeric(quantile_range) || length(quantile_range) != 2L ||
      quantile_range[1] < 0 || quantile_range[2] > 1 ||
      quantile_range[1] >= quantile_range[2]) {
    stop("`quantile_range` must be a numeric vector of length 2 with 0 <= lower < upper <= 1.",
         call. = FALSE)
  }
  # Validate summary_max_depth (rpart allows maxdepth up to 30)
  if (!is.numeric(summary_max_depth) || length(summary_max_depth) != 1L ||
      is.na(summary_max_depth) || summary_max_depth < 1 ||
      summary_max_depth > 30 || summary_max_depth != round(summary_max_depth)) {
    stop("`summary_max_depth` must be a single integer between 1 and 30.",
         call. = FALSE)
  }
  summary_max_depth <- as.integer(summary_max_depth)
  # Validate min_kept
  if (!is.numeric(min_kept) || length(min_kept) != 1L || is.na(min_kept) ||
      min_kept < 0 || min_kept >= 1) {
    stop("`min_kept` must be a single number in [0, 1).", call. = FALSE)
  }
  user_supplied_objective <- !is.null(global_objective_fn)
  if (!user_supplied_objective) {
    global_objective_fn <- if (isTRUE(generalizability_path)) {
      objective_generalizability
    } else {
      objective_default
    }
  }
  if (!is.function(global_objective_fn))
    stop("`global_objective_fn` must be a function(D) -> numeric.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  # ---- Step 1: generalizability_path vs general optimization --------------------
  if (generalizability_path) {
    outcome_col   <- "Y"
    treatment_col <- "Tr"
    sample_col    <- "S"

    needed <- c(outcome_col, treatment_col, sample_col)
    missing_cols <- setdiff(needed, names(data))
    if (length(missing_cols) > 0L) {
      stop(
        "With generalizability_path = TRUE, `data` must contain columns ",
        "'Y' (outcome), 'Tr' (treatment indicator), and 'S' (sample indicator). ",
        "Missing: ", paste(missing_cols, collapse = ", "), ".",
        call. = FALSE
      )
    }

    if (!is.numeric(data[[outcome_col]])) {
      stop("In generalizability_path mode, `Y` must be numeric.", call. = FALSE)
    }

    # Coerce S first (needed for subsequent validation)
    data[[sample_col]] <- coerce01(data[[sample_col]], allow_na = FALSE)
    trial_mask <- data[[sample_col]] == 1

    # Check for NA in Y among trial units (S == 1)
    if (any(is.na(data[[outcome_col]][trial_mask]))) {
      stop("Y (outcome) contains missing values among trial units (S == 1). Please handle NA before training.", call. = FALSE)
    }

    # Treatment: allow NA for target population (S == 0), but not for trial (S == 1)
    data[[treatment_col]] <- coerce01(data[[treatment_col]], allow_na = TRUE)
    if (any(is.na(data[[treatment_col]][trial_mask]))) {
      stop("Tr (treatment) contains missing values among trial units (S == 1). Please handle NA before training.", call. = FALSE)
    }

    # ---- Positivity trimming (Crump et al., 2009) ----------------------------
    positivity_trim_info <- NULL
    if (!is.null(positivity_trim)) {
      trim_result <- trim_positivity_violations(
        data          = data,
        sample_col    = sample_col,
        outcome_col   = outcome_col,
        treatment_col = treatment_col,
        threshold     = positivity_trim
      )
      if (trim_result$n_trimmed > 0L) {
        message(sprintf(
          paste0("Positivity trimming: removed %d target observation(s) with ",
                 "P(S=1|X) < %g (method = %s)."),
          trim_result$n_trimmed, trim_result$alpha,
          if (identical(positivity_trim, "crump")) "crump" else "fixed"
        ))
      }
      positivity_trim_info <- list(
        method    = if (identical(positivity_trim, "crump")) "crump" else "fixed",
        alpha     = trim_result$alpha,
        n_trimmed = trim_result$n_trimmed
      )
      data <- trim_result$data
      # Re-derive trial_mask after trimming
      trial_mask <- data[[sample_col]] == 1
    }

    scores <- compute_transport_scores(
      data     = data,
      outcome  = outcome_col,
      treatment = treatment_col,
      sample   = sample_col
    )
    v_all    <- scores$v
    vsq_all  <- scores$vsq
    sel_all  <- scores$sel_odds
    S_full  <- data[[sample_col]]

    analysis_idx <- which(S_full == 1L)
    if (length(analysis_idx) == 0L)
      stop("No observations with S == 1 (trial units) found in `data`.", call. = FALSE)

    testing_data <- data[analysis_idx, , drop = FALSE]
    v_vals       <- v_all[analysis_idx]
    vsq_vals     <- vsq_all[analysis_idx]
    sel_vals     <- sel_all[analysis_idx]
    S_vec        <- rep(1L, length(analysis_idx))
    include_vsq  <- TRUE   # kept so user-supplied objectives can still read D$vsq
  } else {
    positivity_trim_info <- NULL
    testing_data <- data
    n <- nrow(testing_data)

    if ("v" %in% names(testing_data)) {
      v_vals <- testing_data[["v"]]
    } else {
      v_vals <- rep(0, n)
    }

    sel_vals <- NULL
    if ("vsq" %in% names(testing_data)) {
      vsq_vals    <- testing_data[["vsq"]]
      include_vsq <- TRUE
    } else if ("v" %in% names(testing_data) && is.numeric(v_vals)) {
      # Only `v` supplied: do NOT pass a precomputed vsq to the objective.
      # objective_default() then centers v at the mean of the currently
      # included units on every evaluation. The fixed-mean version below is
      # used only to guide feature-importance sampling.
      mu_v        <- mean(v_vals, na.rm = TRUE)
      vsq_vals    <- (v_vals - mu_v)^2
      include_vsq <- FALSE
    } else {
      # Fallback: vsq = 1. The default objective is then 1 / sqrt(n_kept),
      # which is minimized by keeping every observation.
      vsq_vals    <- rep(1, n)
      include_vsq <- TRUE
      if (!user_supplied_objective) {
        message(
          "No `vsq` or `v` column and no `global_objective_fn` supplied; using ",
          "vsq = 1. The default objective then reduces to 1 / sqrt(n_kept), ",
          "which favors keeping all observations. Supply `vsq`, `v`, or a ",
          "custom `global_objective_fn` for a meaningful optimization."
        )
      }
    }

    if ("S" %in% names(testing_data)) {
      S_vec <- coerce01(testing_data[["S"]], allow_na = TRUE)
    } else {
      S_vec <- rep(1L, n)
    }
  }

  n <- nrow(testing_data)

  # Minimum kept-set size for the default generalizability objective: without
  # it, keeping a few units with similar scores drives the estimated SE to ~0.
  if (isTRUE(generalizability_path) && !user_supplied_objective) {
    n_min_kept <- max(2L, as.integer(ceiling(min_kept * n)))
    base_objective <- global_objective_fn
    global_objective_fn <- function(D) {
      if (sum(D$w == 1, na.rm = TRUE) < n_min_kept) return(Inf)
      base_objective(D)
    }
  }

  # ---- Step 2: Prepare covariates for splitting --------------------------
  drop_cols <- intersect(
    names(testing_data),
    c("Y", "Tr", "S", "v", "vsq", "w", "lX", "sel_odds")
  )
  X_df <- testing_data[, setdiff(names(testing_data), drop_cols), drop = FALSE]

  if (ncol(X_df) == 0L) {
    # no covariates to split on; construct leaf-only trees
    features      <- "leaf"
    split_feature <- stats::setNames(1, "leaf")
  } else {
    features      <- c("leaf", colnames(X_df))
    split_feature <- compute_split_feature(X_df, vsq_vals)
    split_feature <- c(leaf = leaf_proba, split_feature)
    split_feature <- split_feature / sum(split_feature)
    names(split_feature) <- features
  }

  # ---- Step 3: Build forest ----------------------------------------------
  # Working-state builder: columns read by the objective functions
  make_state <- function(with_w) {
    D <- data.frame(X_df, v = v_vals, stringsAsFactors = FALSE)
    if (include_vsq)       D$vsq      <- vsq_vals
    if (!is.null(sel_vals)) D$sel_odds <- sel_vals
    if (with_w)            D$w        <- rep(1, n)
    D$S <- S_vec
    D
  }
  D_forest <- make_state(with_w = FALSE)

  w_forest <- vector("list", num_trees)
  for (t_idx in seq_len(num_trees)) {
    D_tree <- make_state(with_w = TRUE)
    rownames(D_tree) <- as.character(seq_len(n))
    rownames(X_df)   <- as.character(seq_len(n))

    grow_one <- function() {
      split_node(
        split_feature        = split_feature,
        X                    = D_tree,
        D                    = D_tree,
        parent_loss          = Inf,
        depth                = 0L,
        explore_proba        = explore_proba,
        choose_feature_fn    = choose_feature,
        reduce_weight_fn     = reduce_weight,
        global_objective_fn  = global_objective_fn,
        max_depth            = max_depth,
        min_leaf_n           = min_leaf_n,
        max_rejects_per_node = max_rejects_per_node,
        split_strategy       = split_strategy,
        quantile_range       = quantile_range
      )
    }

    tree_res <- if (is.null(seed)) grow_one() else withr::with_seed(seed + t_idx, grow_one())
    D_updated <- tree_res$D
    D_forest[[paste0("w_tree_", t_idx)]] <- D_updated$w
    w_forest[[t_idx]] <- tree_res
  }

  # ---- Step 4: Rashomon set selection ------------------------------------
  obj_values   <- vapply(w_forest, function(res) res[["local objective"]], numeric(1))
  tree_indices <- seq_along(w_forest)

  if (top_k_trees) {
    if (k > num_trees) {
      warning("k > num_trees; using k = num_trees.", call. = FALSE)
      k <- num_trees
    }
    ord <- order(obj_values, na.last = NA)
    rashomon_set <- utils::head(ord, k)
  } else {
    cutoff_val <- if (identical(cutoff, "baseline")) {
      D_base <- D_forest
      D_base$w <- 1
      global_objective_fn(D_base)
    } else as.numeric(cutoff)
    if (!is.finite(cutoff_val)) cutoff_val <- Inf
    rashomon_set <- tree_indices[obj_values < cutoff_val]
  }

  if (length(rashomon_set) == 0L) {
    if (identical(cutoff, "baseline")) {
      message("No trees selected into Rashomon set.")
    } else {
      warning("No trees selected into Rashomon set.", call. = FALSE)
    }
  }
  not_in_set <- setdiff(tree_indices, rashomon_set)

  D_rash <- D_forest
  if (length(not_in_set) > 0L) {
    drop_cols2 <- paste0("w_tree_", not_in_set)
    keep_cols  <- setdiff(names(D_rash), drop_cols2)
    D_rash     <- D_rash[, keep_cols, drop = FALSE]
  }

  weight_cols <- grep("^w_tree_", names(D_rash), value = TRUE)
  D_weights   <- if (length(weight_cols) > 0L) D_rash[, weight_cols, drop = FALSE] else data.frame()
  if (ncol(D_weights) > 0L) {
    votes_mat         <- as.matrix(D_weights)  # n x |rashomon_set| matrix of 0/1 votes
    D_rash$vote_count <- rowSums(votes_mat)

    if (is.function(vote_threshold)) {
      # User-supplied aggregation function
      w_opt_raw <- tryCatch(
        vote_threshold(votes_mat),
        error = function(e) stop(
          "vote_threshold() raised an error: ", conditionMessage(e), call. = FALSE
        )
      )
      # Validate BEFORE coercing: catches non-binary values like 0.5 that
      # would silently become 0L after as.integer(), hiding the error.
      if (length(w_opt_raw) != nrow(votes_mat) ||
          !all(w_opt_raw %in% c(0, 1, 0L, 1L), na.rm = TRUE))
        stop(
          "vote_threshold() must return an integer vector of 0s and 1s ",
          "with one element per observation.", call. = FALSE
        )
      if (!is.integer(w_opt_raw)) w_opt_raw <- as.integer(w_opt_raw)
      D_rash$w_opt <- w_opt_raw
    } else if (is.character(vote_threshold)) {
      if (vote_threshold == "majority") {
        D_rash$w_opt <- as.integer(rowMeans(votes_mat) > 2 / 3)

      } else if (vote_threshold == "mean") {
        row_means <- rowMeans(votes_mat)
        tied      <- row_means == 0.5
        if (any(tied)) {
          stop(
            sum(tied), " observation(s) have an exact mean vote of 0.5 across the ",
            "Rashomon set -- the binary weight is undefined. ",
            "Consider using an odd number of trees, a different vote_threshold, ",
            "or a custom aggregation function.",
            call. = FALSE
          )
        }
        D_rash$w_opt <- as.integer(row_means > 0.5)

      } else if (vote_threshold == "median") {
        row_medians <- apply(votes_mat, 1, stats::median)
        tied        <- row_medians == 0.5
        if (any(tied)) {
          stop(
            sum(tied), " observation(s) have an exact median vote of 0.5 across the ",
            "Rashomon set -- the binary weight is undefined. ",
            "Consider using an odd number of trees, a different vote_threshold, ",
            "or a custom aggregation function.",
            call. = FALSE
          )
        }
        D_rash$w_opt <- as.integer(row_medians > 0.5)

      } else {
        stop("Unrecognised vote_threshold string: '", vote_threshold,
             "'. Must be one of \"majority\", \"mean\", or \"median\".", call. = FALSE)
      }
    } else {
      # Numeric threshold: use rowMeans (fraction of trees voting 1)
      row_means <- rowMeans(votes_mat)
      tied      <- row_means == vote_threshold
      if (any(tied)) {
        stop(
          sum(tied), " observation(s) have a mean vote exactly equal to the numeric ",
          "threshold (", vote_threshold, ") -- the binary weight is undefined. ",
          "Adjust the threshold or the number of trees.",
          call. = FALSE
        )
      }
      D_rash$w_opt <- as.integer(row_means > vote_threshold)
    }
  } else {
    D_rash$w_opt      <- integer(nrow(D_rash))
    D_rash$vote_count <- integer(nrow(D_rash))
  }

  # ---- Step 5: Summary tree ----------------------------------------------
  final_classifier <- NULL
  if (ncol(X_df) > 0L) {
    w_classes <- unique(stats::na.omit(D_rash$w_opt))

    if (length(w_classes) == 2L && all(w_classes %in% c(0L, 1L))) {
      final_classifier <- characterize_tree(X_df, D_rash$w_opt,
                                            max_depth = summary_max_depth)
    } else if (length(w_classes) == 1L) {

      # <-- This message is what your test-root.R:312 is expecting
      message("No summary tree available to plot (single-class w_opt).")

      # Keep your informative diagnostics too
      if (w_classes == 1L) {
        message("ROOT: All observations have w_opt = 1 (all kept). No subgroup distinction found.")
      } else {
        warning(
          "ROOT: All observations have w_opt = 0 (all removed). This may indicate extreme ",
          "underrepresentation or model issues. The weighted estimate will be undefined.",
          call. = FALSE
        )
      }

    } else {
      message("No summary tree available to plot (w_opt not binary with two classes or no covariates).")
    }
  }

  # ---- Step 6: Estimands (only for generalizability_path = TRUE) ----------------
  est_list <- NULL
  if (generalizability_path) {
    in_S  <- rep(TRUE, n)  # we already subset to S == 1
    v_sel <- v_vals        # IPW scores, already scaled by pi/(1-pi)
    r_sel <- sel_vals      # inverse selection odds (1 - pi(X)) / pi(X)
    w_sel <- if ("w_opt" %in% names(D_rash)) D_rash$w_opt[in_S] else rep(1L, length(v_sel))

    ## Unweighted STATE (sample target average treatment effect)
    est_label_unw <- "STATE"
    ok_unw        <- !is.na(v_sel)
    n_eff_unw     <- sum(ok_unw)
    mu_unw        <- if (n_eff_unw > 0) mean(v_sel[ok_unw]) else NA_real_
    se_unw        <- if (n_eff_unw > 1) {
      diffs2 <- (v_sel[ok_unw] - mu_unw)^2
      sqrt(sum(diffs2) / (n_eff_unw * (n_eff_unw - 1)))
    } else NA_real_

    ## Weighted WTATE
    est_label_w  <- "WTATE"
    is_binary_w  <- all(is.na(w_sel) | (w_sel %in% c(0, 1)))
    den_w_any    <- sum(w_sel, na.rm = TRUE)
    mu_w_any     <- if (den_w_any > 0) sum(w_sel * v_sel, na.rm = TRUE) / den_w_any else NA_real_

    if (is_binary_w) {
      ok_w  <- ok_unw & !is.na(w_sel) & (w_sel == 1L)
      n_A   <- sum(ok_w)
      den_w <- sum(w_sel == 1L, na.rm = TRUE)
      if (n_A > 0) {
        # Re-normalize from pi/(1-pi) to pi_w/(1-pi_w),
        # pi_w = P(S = 1 | w(X) = 1); see objective_generalizability().
        kappa_w <- wtate_kappa(r_sel, ok_w, ok_unw)
        vbar_A  <- mean(v_sel[ok_w])
        mu_w    <- kappa_w * vbar_A
        se_w <- if (n_A > 1) {
          diffs2_A <- (v_sel[ok_w] - vbar_A)^2
          kappa_w * sqrt(sum(diffs2_A) / (n_A * (n_A - 1)))
        } else NA_real_

        se_w_note <- paste0(
          "WTATE = kappa_w * vbar_A; SE = kappa_w * sqrt( sum_{A} (v_i - vbar_A)^2 / ( n_A * (n_A - 1) ) ).\n",
          "  Here A = { i : S_i = 1, w_i = 1 }, n_A = |A|, v_i are the IPW transported scores, vbar_A is\n",
          "  their mean, and kappa_w = [n_A / sum_A r_i] / [n_1 / sum_{S=1} r_i] with r_i = (1 - pi(X_i)) / pi(X_i).\n",
          "  The SE treats w and the nuisance models as fixed.\n"
        )
        if (user_supplied_objective) {
          se_w_note <- paste0(
            se_w_note,
            "\nYou supplied a custom global_objective_fn; please verify this SE matches your ",
            "estimand, or perform additional variance analysis (e.g., bootstrap or ",
            "influence-function methods)."
          )
        }
      } else {
        mu_w     <- NA_real_
        se_w     <- NA_real_
        kappa_w  <- NA_real_
        den_w    <- 0L
        n_A      <- 0L
        se_w_note <- "SE omitted: no kept observations (A = { i : w_i = 1 } is empty)."
      }
    } else {
      # nocov start
      # Defensive only: w_opt is validated as binary during vote aggregation.
      ok_nb    <- ok_unw & !is.na(w_sel)
      kappa_w  <- (sum(w_sel[ok_nb]) / sum(w_sel[ok_nb] * r_sel[ok_nb])) /
        (sum(ok_unw) / sum(r_sel[ok_unw]))
      mu_w     <- kappa_w * mu_w_any
      den_w    <- den_w_any
      n_A      <- NA_integer_
      se_w     <- NA_real_
      se_w_note <- paste0(
        "SE omitted: non-binary weights detected; the subset-mean SE above is not appropriate.",
        " If you intentionally use a custom global_objective_fn with non-binary weights, ",
        "please conduct further variance analysis tailored to that objective."
      )
      warning("Non-binary w_opt detected; omitting weighted SE.", call. = FALSE)
      # nocov end
    }

    if (verbose) {
      message(sprintf("%s (unweighted) = %.6f, SE = %.6f", est_label_unw, mu_unw, se_unw))
      if (is.finite(se_w) || (!is.na(se_w) && !is.nan(se_w))) {
        message(sprintf("%s (weighted)   = %.6f, SE = %.6f", est_label_w, mu_w, se_w))
        message(paste0("  Note: ", se_w_note))
      } else {
        message(sprintf("%s (weighted)   = %.6f", est_label_w, mu_w))
        message(paste0("  Note: ", se_w_note))
      }
    }

    est_list <- list(
      estimand_unweighted = est_label_unw,
      value_unweighted    = mu_unw,
      se_unweighted       = se_unw,
      estimand_weighted   = est_label_w,
      value_weighted      = mu_w,
      se_weighted         = se_w,
      se_weighted_note    = se_w_note,
      n_analysis          = sum(in_S),
      sum_w               = den_w,
      kappa_w             = kappa_w,
      n_A                 = if (is_binary_w) n_A else NA_integer_
    )
  }

  res <- list(
    D_rash                  = D_rash,
    D_forest                = D_forest,
    w_forest                = w_forest,
    rashomon_set            = rashomon_set,
    global_objective_fn     = global_objective_fn,
    user_supplied_objective = user_supplied_objective,
    vote_threshold          = vote_threshold,
    f                       = final_classifier,
    testing_data            = testing_data,
    estimate                = est_list,
    generalizability_path   = generalizability_path,
    positivity_trim_info    = positivity_trim_info,
    split_strategy          = split_strategy
  )
  class(res) <- c("ROOT", "list")
  res
}
