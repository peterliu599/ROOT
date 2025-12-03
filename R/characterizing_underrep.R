#' Characterize under-represented subgroups (wraps \code{ROOT})
#'
#' Combines an RCT (\code{S = 1}) and a target dataset (\code{S = 0}), then calls \code{ROOT()} to
#' learn a weighted tree that identifies subgroups with different representation
#' in the target population compared to the trial.
#'
#' @section Abbreviations:
#' RCT means randomized clinical trial. ATE means Average Treatment Effect.
#' WTATE means Weighted Transported ATE. SE means Standard Error.
#'
#' @param DataRCT A \code{data.frame} containing the randomized clinical trial data.
#'   Must include treatment, outcome, and covariate columns.
#' @param covariateColName_RCT A \code{character} vector of covariate column names in \code{DataRCT}.
#' @param trtColName_RCT A \code{character(1)} naming the treatment column in \code{DataRCT} with values \code{0} or \code{1}.
#' @param outcomeColName_RCT A \code{character(1)} naming the outcome column in \code{DataRCT}.
#' @param DataTarget A \code{data.frame} containing the target population covariates only.
#' @param covariateColName_TargetData A \code{character} vector of covariate column names in \code{DataTarget}.
#' @param leaf_proba A \code{numeric(1)} giving the probability for the \code{"leaf"} option in \code{ROOT} tree growth. Default \code{0.25}.
#' @param seed An \code{integer(1)} seed for reproducibility. Default \code{123}.
#' @param num_trees An \code{integer(1)} number of trees to grow. Default \code{10}.
#' @param vote_threshold A \code{numeric(1)} in \code{(0.5, 1]} for the majority vote threshold. Default \code{2/3}.
#' @param explore_proba A \code{numeric(1)} exploration probability. Default \code{0.05}.
#' @param feature_est Either a \code{character(1)} in \code{c("Ridge","GBM")} or a \code{function(X, y, ...)} that returns a named nonnegative \code{numeric} vector of importances.
#' @param feature_est_args A named \code{list} of extra arguments passed to the user supplied \code{feature_est} function.
#' @param top_k_trees A \code{logical(1)}. If \code{TRUE}, selects top \code{k} trees by objective; otherwise uses \code{cutoff}. Default \code{FALSE}.
#' @param k An \code{integer(1)} number of trees used when \code{top_k_trees = TRUE}. Default \code{10}.
#' @param cutoff A \code{numeric(1)} or the value \code{"baseline"} used as the Rashomon set cutoff when \code{top_k_trees = FALSE}.
#' @param verbose A \code{logical(1)}. If \code{TRUE}, prints progress and estimand summaries. Default \code{FALSE}.
#' @param global_objective_fn A \code{function} with signature \code{function(D) -> numeric} to minimize. Default \code{objective_default}.
#' @param keep_threshold Unused; kept for backward compatibility. A \code{numeric(1)} if provided.
#' @param lX_threshold Unused; kept for backward compatibility. A \code{numeric(1)} if provided.
#'
#' @return A \code{characterizing_underrep} S3 object (a \code{list}) with components:
#'   \item{root}{The resulting \code{ROOT} object returned by \code{ROOT()}.}
#'   \item{combined}{A \code{data.frame} with the stacked RCT and target data used for analysis.}
#'   \item{leaf_summary}{A \code{data.frame} of terminal node summaries with rules, counts, percentages, and labels when a summary tree exists; otherwise \code{NULL}.}
#'
#' @references
#' Parikh, H., Ross, R. K., Stuart, E., and Rudolph, K. E. (2025).
#' Who Are We Missing?: A Principled Approach to Characterizing the Underrepresented Population.
#' \emph{Journal of the American Statistical Association}, 1–32.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(diabetes_data)
#'
#' # Split into Trial (S=1) and Target (S=0)
#' trial  <- subset(diabetes_data, S == 1)
#' target <- subset(diabetes_data, S == 0)
#'
#' # Run characterization
#' res <- characterizing_underrep(
#'   DataRCT = trial,
#'   covariateColName_RCT = c("Race_Black", "Sex_Male", "DietYes", "Age45"),
#'   trtColName_RCT = "Tr",
#'   outcomeColName_RCT = "Y",
#'   DataTarget = target,
#'   covariateColName_TargetData = c("Race_Black", "Sex_Male", "DietYes", "Age45"),
#'   seed = 123
#' )
#' }
#' @export
characterizing_underrep <- function(
    DataRCT,
    covariateColName_RCT,
    trtColName_RCT,
    outcomeColName_RCT,
    DataTarget,
    covariateColName_TargetData,
    leaf_proba = 0.25,
    seed = 123,
    num_trees = 10,
    vote_threshold = 2 / 3,
    explore_proba = 0.05,
    feature_est = "Ridge",
    feature_est_args = list(),
    top_k_trees = FALSE,
    k = 10,
    cutoff = "baseline",
    verbose = FALSE,
    global_objective_fn = objective_default,
    keep_threshold = 0.50,
    lX_threshold = NULL
) {
  # 1. Input Validation
  if (!is.data.frame(DataRCT) || !is.data.frame(DataTarget)) {
    stop("DataRCT and DataTarget must be data.frames.", call. = FALSE)
  }
  if (!all(c(trtColName_RCT, outcomeColName_RCT) %in% names(DataRCT))) {
    stop("trtColName_RCT and/or outcomeColName_RCT not found in DataRCT.", call. = FALSE)
  }
  if (!all(covariateColName_RCT %in% names(DataRCT))) {
    miss <- setdiff(covariateColName_RCT, names(DataRCT))
    stop("Missing RCT covariates in DataRCT: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!all(covariateColName_TargetData %in% names(DataTarget))) {
    miss <- setdiff(covariateColName_TargetData, names(DataTarget))
    stop("Missing target covariates in DataTarget: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # 2. Harmonize Covariates
  # We take the intersection of covariates available in both datasets
  common_covs <- intersect(covariateColName_RCT, covariateColName_TargetData)
  if (length(common_covs) == 0L) {
    stop("No overlapping covariate names between RCT and Target.", call. = FALSE)
  }

  # Warn if some requested covariates are being dropped
  if (!identical(sort(common_covs), sort(covariateColName_RCT)) ||
      !identical(sort(common_covs), sort(covariateColName_TargetData))) {
    warning("Using intersection of covariates: ", paste(common_covs, collapse = ", "))
  }

  X_rct <- DataRCT[, common_covs, drop = FALSE]
  X_tgt <- DataTarget[, common_covs, drop = FALSE]

  # 3. Construct Combined Dataset
  # S = 1 for RCT, S = 0 for Target
  combined_rct <- data.frame(
    X_rct,
    Y  = DataRCT[[outcomeColName_RCT]],
    Tr = DataRCT[[trtColName_RCT]],
    S  = 1L
  )
  combined_tgt <- data.frame(
    X_tgt,
    Y  = NA_real_,
    Tr = NA_integer_,
    S  = 0L
  )
  combined <- rbind(combined_rct, combined_tgt)
  rownames(combined) <- NULL

  # 4. Call ROOT (Logic Only)
  # We strip away all plotting arguments here.
  root_out <- ROOT(
    data = combined,
    outcome = "Y",
    treatment = "Tr",
    sample = "S",
    leaf_proba = leaf_proba,
    seed = seed,
    num_trees = num_trees,
    vote_threshold = vote_threshold,
    explore_proba = explore_proba,
    feature_est = feature_est,
    feature_est_args = feature_est_args,
    top_k_trees = top_k_trees,
    k = k,
    cutoff = cutoff,
    verbose = verbose,
    global_objective_fn = global_objective_fn
  )

  # 5. Summarize Terminal Nodes
  # Extract rules from the summary tree (if it exists)
  leaf_summary <- NULL
  if (!is.null(root_out$f)) {
    f <- root_out$f
    frm <- f$frame
    is_leaf <- frm$var == "<leaf>"

    # Get paths (rules)
    # path.rpart returns a list of paths for the specified nodes
    paths <- rpart::path.rpart(f, nodes = as.numeric(rownames(frm)[is_leaf]), print.it = FALSE)
    path_text <- vapply(paths, function(p) paste(p, collapse = " & "), character(1))

    # Predicted class label per leaf
    ylv <- if (!is.null(f$ylevels)) f$ylevels else levels(stats::model.frame(f)$w)
    pred_class <- ylv[frm$yval[is_leaf]]

    n_in_leaf <- frm$n[is_leaf]
    total_n <- sum(n_in_leaf)
    pct <- if (total_n > 0) round(100 * n_in_leaf / total_n, 1) else NA_real_

    # Human-readable label
    label <- ifelse(pred_class %in% c("0", "0.0", "0L", "FALSE"),
                    "Under-represented (drop, w=0)",
                    "Represented (keep, w=1)")

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

  # 6. Construct and Return S3 Object
  res <- list(
    root         = root_out,
    combined     = combined,
    leaf_summary = leaf_summary
  )
  class(res) <- c("characterizing_underrep", "list")
  return(res)
}
