#' Characterize Underrepresented Subgroups
#'
#' @description
#' A high-level wrapper around \code{\link{ROOT}()} for identifying and
#' characterizing subgroups that are \emph{insufficiently represented} in a
#' randomized controlled trial (RCT) relative to a target population. The
#' function returns an interpretable decision tree describing which subgroups
#' should be included (\eqn{w(X) = 1}) or excluded (\eqn{w(X) = 0}) from the
#' analysis, along with the corresponding target treatment average effect estimates.
#'
#' @section What does "underrepresented" mean?:
#' In the context of generalizing treatment effects from a trial to a target
#' population, a subgroup is considered \strong{underrepresented} (or
#' \emph{insufficiently represented}) when it occupies a region of the
#' covariate space that both (a) has limited overlap between the trial and the
#' target population, and (b) exhibits heterogeneous treatment effects.
#'
#' Formally, the contribution of a unit with covariates \eqn{X = x} to the
#' variance of the target average treatment effect (TATE) estimator depends on
#' both the selection ratio \eqn{\ell(x) = P(S=1 \mid X=x) / P(S=0 \mid X=x)}
#' and the conditional average treatment effect. Subgroups where \eqn{\ell(x)}
#' is small, and conditional average treatment effect deviates from the
#' overall TATE, contribute disproportionately to estimator variance. These are
#' the subgroups that \code{characterizing_underrep()} identifies and
#' characterizes. The sample average treatment effect (SATE) is a
#' finite sample equivalent version of the TATE.
#'
#' @section The generalizability workflow:
#' When \code{generalizability_path = TRUE}, this function implements the
#' two-stage approach of Parikh et al. (2025):
#'
#' \enumerate{
#'   \item \strong{Design stage:} ROOT learns binary weights \eqn{w(X)} that
#'     minimize the variance of the weighted target average treatment effect
#'     (WTATE) estimator, subject to interpretability constraints (tree
#'     structure). The resulting decision tree characterizes which subgroups
#'     are well-represented (\eqn{w = 1}) and which are underrepresented
#'     (\eqn{w = 0}).
#'   \item \strong{Analysis stage:} The WTATE is estimated on the refined
#'     target population that excludes the underrepresented subgroups. This
#'     estimand trades some generality for greater precision and credibility.
#' }
#'
#' The key estimands are:
#' \itemize{
#'   \item \strong{SATE} (Sample Average Treatment Effect): the treatment
#'     effect for the full target population based on the trial sample, which
#'     may be imprecise if certain subgroups are underrepresented. It is a
#'     finite sample equivalent version of the TATE.
#'   \item \strong{WTATE} (Weighted Target Average Treatment Effect): the
#'     treatment effect restricted to the sufficiently represented
#'     subpopulation, estimated with lower variance.
#' }
#'
#' @section General optimization mode:
#' When \code{generalizability_path = FALSE}, this function behaves as a
#' convenience wrapper around \code{\link{ROOT}()} for arbitrary binary weight
#' optimization. The user can supply a custom objective function via
#' \code{global_objective_fn}; ROOT will learn an interpretable tree-based
#' weight function minimizing that objective. See
#' \code{vignette("optimization_path_example")} for an example.
#'
#' @section Data requirements:
#' When \code{generalizability_path = TRUE}, \code{data} must contain the
#' following standardized columns:
#' \itemize{
#'   \item \code{Y}: numeric outcome variable.
#'   \item \code{Tr}: binary treatment indicator (0 = control, 1 = treated).
#'   \item \code{S}: binary sample indicator (1 = trial/RCT, 0 = target
#'     population).
#' }
#' All remaining columns are treated as pretreatment covariates \eqn{X}
#' available for splitting.
#'
#' @param data A \code{data.frame} containing covariates and, in
#'   generalizability mode, also columns \code{Y}, \code{Tr}, and \code{S}.
#' @param global_objective_fn A function with signature
#'   \code{function(D) -> numeric} scoring the entire state and minimized by
#'   ROOT. If \code{NULL} (the default), a variance-based objective is used
#'   that targets the precision of the WTATE estimator.
#' @param generalizability_path Logical. If \code{TRUE}, runs the full
#'   generalizability analysis (design + analysis stages) and expects columns
#'   \code{Y}, \code{Tr}, and \code{S} in \code{data}. If \code{FALSE}, runs
#'   ROOT in general optimization mode. Default \code{FALSE}.
#' @param leaf_proba A \code{numeric(1)} tuning parameter that increases the
#'   chance a node stops splitting by selecting a synthetic \code{"leaf"}
#'   feature. Internally, the probability of choosing \code{"leaf"} is
#'   \code{leaf_proba / (1 + leaf_proba)} (assuming the covariate
#'   probabilities sum to 1). Higher values produce shallower, more
#'   interpretable trees. Default \code{0.25}.
#' @param seed Random seed for reproducibility.
#' @param num_trees Number of trees to grow in the ROOT forest. More trees
#'   explore the tree space more thoroughly but increase computation time.
#' @param vote_threshold Controls how Rashomon-set tree votes are aggregated
#'   into \code{w_opt}. Accepts a numeric threshold in \code{(0, 1]}, one of
#'   \code{"majority"} / \code{"mean"} / \code{"median"}, or a custom
#'   \code{function(votes) -> integer vector}. See \code{\link{ROOT}} for
#'   full details. Default \code{2/3} (majority vote).
#' @param explore_proba Exploration probability in tree growth. Controls the
#'   explore-exploit trade-off: with probability \code{explore_proba}, a leaf
#'   is assigned a random weight; otherwise, it receives the greedy optimal
#'   weight. Default \code{0.05}.
#' @param feature_est Either \code{"Ridge"}, \code{"GBM"}, or a custom
#'   feature importance function. Used to bias which covariates are selected
#'   for splitting. Default \code{"Ridge"}.
#' @param feature_est_args List of extra arguments passed to
#'   \code{feature_est} when it is a function.
#' @param top_k_trees Logical; if \code{TRUE}, selects the top \code{k} trees
#'   by objective value for the Rashomon set. If \code{FALSE}, uses
#'   \code{cutoff} instead. Default \code{FALSE}.
#' @param k Number of trees to retain when \code{top_k_trees = TRUE}. Default
#'   \code{10}.
#' @param cutoff Numeric or \code{"baseline"}. When \code{top_k_trees = FALSE},
#'   trees with objective values below this cutoff are included in the
#'   Rashomon set. \code{"baseline"} uses the objective value at
#'   \eqn{w \equiv 1} (no subgroups excluded). Default \code{"baseline"}.
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
#' @param verbose Logical; if \code{TRUE}, prints the unweighted and weighted
#'   estimands with standard errors. Default \code{FALSE}.
#'
#' @return A \code{characterizing_underrep} S3 object (a \code{list}) with:
#'   \item{root}{The \code{\link{ROOT}} object returned by \code{ROOT()}.
#'     Contains the Rashomon set, weight assignments (\code{root$D_rash$w_opt}),
#'     summary tree (\code{root$f}), and (in generalizability mode) treatment
#'     effect estimates (\code{root$estimate}).}
#'   \item{combined}{The input \code{data}.}
#'   \item{leaf_summary}{A \code{data.frame} with one row per terminal node of
#'     the summary (characteristic) tree, giving the decision rule, predicted
#'     weight, sample size, and a label indicating whether the
#'     subgroup is "Represented (keep, w = 1)" or "Under-represented
#'     (drop, w = 0)". \code{NULL} if no summary tree was produced.}
#'
#' @references
#' Parikh H, Ross RK, Stuart E, Rudolph KE (2025).
#' "Who Are We Missing?: A Principled Approach to Characterizing the
#' Underrepresented Population."
#' \emph{Journal of the American Statistical Association}.
#' \doi{10.1080/01621459.2025.2495319}
#'
#' @seealso
#' \code{\link{ROOT}} for the underlying optimization engine;
#' \code{vignette("generalizability_path_example")} for a detailed worked
#' example of the generalizability workflow;
#' \code{vignette("optimization_path_example")} for general optimization mode.
#'
#' @examples
#' \dontrun{
#' # --- Generalizability analysis ---
#' # diabetes_data has columns Y, Tr, S, and covariates
#' data(diabetes_data, package = "ROOT")
#'
#' char_fit <- characterizing_underrep(
#'   data                  = diabetes_data,
#'   generalizability_path = TRUE,
#'   num_trees             = 20,
#'   top_k_trees           = TRUE,
#'   k                     = 10,
#'   seed                  = 123
#' )
#'
#' # View the characterization tree
#' plot(char_fit)
#'
#' # Inspect which subgroups are underrepresented
#' char_fit$leaf_summary
#'
#' # Treatment effect estimates (SATE and WTATE)
#' char_fit$root$estimate
#' }
#' @export
characterizing_underrep <- function(data,
                                    global_objective_fn   = NULL,
                                    generalizability_path = FALSE,
                                    leaf_proba            = 0.25,
                                    seed                  = 123,
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
                                    verbose               = FALSE) {
  # Data frame check
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!is.logical(generalizability_path) || length(generalizability_path) != 1L) {
    stop("`generalizability_path` must be TRUE or FALSE.", call. = FALSE)
  }

  # If we are in generalizability_path mode, enforce Y / Tr / S presence
  if (isTRUE(generalizability_path)) {
    needed <- c("Y", "Tr", "S")
    missing_cols <- setdiff(needed, names(data))
    if (length(missing_cols) > 0L) {
      stop(
        "For generalizability_path = TRUE, `data` must contain columns: ",
        paste(needed, collapse = ", "),
        ". Missing: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # 1. Call ROOT with the requested path
  root_out <- ROOT(
    data                  = data,
    generalizability_path = generalizability_path,
    leaf_proba            = leaf_proba,
    seed                  = seed,
    num_trees             = num_trees,
    vote_threshold        = vote_threshold,
    explore_proba         = explore_proba,
    feature_est           = feature_est,
    feature_est_args      = feature_est_args,
    top_k_trees           = top_k_trees,
    k                     = k,
    cutoff                = cutoff,
    max_depth             = max_depth,
    min_leaf_n            = min_leaf_n,
    max_rejects_per_node  = max_rejects_per_node,
    verbose               = verbose,
    global_objective_fn   = global_objective_fn
  )

  # 2. Summarize terminal nodes (if characterization tree exists)
  leaf_summary <- NULL
  if (!is.null(root_out$f)) {
    f   <- root_out$f
    frm <- f$frame
    is_leaf <- frm$var == "<leaf>"

    paths <- rpart::path.rpart(
      f,
      nodes    = as.numeric(rownames(frm)[is_leaf]),
      print.it = FALSE
    )
    path_text <- vapply(paths, function(p) paste(p, collapse = " & "), character(1))

    ylv <- if (!is.null(f$ylevels)) f$ylevels else levels(stats::model.frame(f)$w)
    pred_class <- ylv[frm$yval[is_leaf]]

    n_in_leaf <- frm$n[is_leaf]
    total_n   <- sum(n_in_leaf)
    pct       <- if (total_n > 0) round(100 * n_in_leaf / total_n, 1) else NA_real_

    label <- ifelse(
      pred_class %in% c("0", "0.0", "0L", "FALSE"),
      "Under-represented (drop, w = 0)",
      "Represented (keep, w = 1)"
    )

    leaf_summary <- data.frame(
      leaf_id     = rownames(frm)[is_leaf],
      rule        = path_text,
      predicted_w = pred_class,
      n           = n_in_leaf,
      pct         = pct,
      label       = label,
      row.names   = NULL,
      check.names = FALSE
    )
  }

  # 3. Construct and return S3 object
  res <- list(
    root         = root_out,
    combined     = data,
    leaf_summary = leaf_summary
  )
  class(res) <- c("characterizing_underrep", "list")
  res
}
