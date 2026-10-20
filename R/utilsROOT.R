#' Compute Crump et al. (2009) optimal propensity-score trimming threshold
#'
#' Implements the data-driven rule from Corollary 1 of Crump, Hotz, Imbens &
#' Mitnik (2009, \emph{Biometrika}) for selecting the symmetric trimming
#' threshold \eqn{\alpha} that minimizes the asymptotic variance of the average
#' treatment effect estimator under homoskedasticity. The optimal set retains
#' observations with propensity scores in \eqn{[\alpha, 1 - \alpha]}.
#'
#' The threshold solves the fixed-point equation
#' \deqn{\frac{1}{\alpha(1 - \alpha)} = 2 \, E\!\left[\frac{1}{e(X)\{1 - e(X)\}}
#'   \;\middle|\; \frac{1}{e(X)\{1 - e(X)\}} \le \frac{1}{\alpha(1 - \alpha)}\right].}
#'
#' @param e_hat Numeric vector of estimated propensity scores (values in
#'   \eqn{(0, 1)}). In the generalizability context this is
#'   \eqn{\hat P(S = 1 \mid X)}.
#' @param grid_size Integer; resolution of the search grid over
#'   \eqn{(0, 0.5)}. Default 1000.
#'
#' @return A single numeric value \eqn{\alpha \in [0, 0.5]}. A value of 0
#'   indicates no trimming is needed.
#'
#' @references
#' Crump RK, Hotz VJ, Imbens GW, Mitnik OA (2009).
#' "Dealing with Limited Overlap in Estimation of Average Treatment Effects."
#' \emph{Biometrika}, \strong{96}(1), 187--199.
#'
#' @keywords internal
crump_alpha <- function(e_hat, grid_size = 1000L) {
  if (!is.numeric(e_hat) || length(e_hat) == 0L) {
    stop("`e_hat` must be a non-empty numeric vector of propensity scores.",
         call. = FALSE)
  }
  # Remove any boundary values that would produce Inf

  e_hat <- e_hat[is.finite(e_hat) & e_hat > 0 & e_hat < 1]
  if (length(e_hat) == 0L) return(0)

  g <- e_hat * (1 - e_hat)   # e(1-e) for each observation

  # Check if no trimming is needed (Corollary 1 condition):
  # sup 1/{e(1-e)} <= 2 * E[1/{e(1-e)}]
  inv_g <- 1 / g
  if (max(inv_g) <= 2 * mean(inv_g)) return(0)

  alphas <- seq(0.001, 0.499, length.out = grid_size)

  for (a in alphas) {
    threshold <- a * (1 - a)
    keep <- g >= threshold
    if (sum(keep) == 0L) break

    lhs <- 1 / threshold               # 1 / {alpha(1-alpha)}
    rhs <- 2 * mean(inv_g[keep])        # 2 * E[1/{e(1-e)} | kept]

    if (lhs <= rhs) return(a)
  }

  0.5
}


#' Trim target observations for positivity violations
#'
#' Removes target-population rows (\code{S == 0}) whose estimated trial
#' selection probability \eqn{\hat P(S = 1 \mid X)} falls below a threshold,
#' addressing structural positivity violations before ROOT optimization.
#'
#' When \code{threshold = "crump"}, the data-driven cutoff from
#' \code{\link{crump_alpha}()} is used. Otherwise a fixed numeric threshold
#' (e.g. 0.05) is applied directly to \eqn{\hat P(S = 1 \mid X)}.
#'
#' @param data A \code{data.frame} containing columns \code{S} (sample
#'   indicator, 1 = trial, 0 = target) plus covariates.
#' @param sample_col Character name of the sample indicator column.
#' @param outcome_col Character name of the outcome column (excluded from the
#'   sampling model).
#' @param treatment_col Character name of the treatment column (excluded from
#'   the sampling model).
#' @param threshold Either \code{"crump"} for data-driven trimming (see
#'   \code{\link{crump_alpha}()}), or a numeric value in \eqn{(0, 0.5)}
#'   specifying the lower bound on \eqn{P(S = 1 \mid X)} below which target
#'   observations are removed.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{data}}{The trimmed \code{data.frame}.}
#'     \item{\code{n_trimmed}}{Integer count of removed target rows.}
#'     \item{\code{alpha}}{The threshold actually used.}
#'     \item{\code{pi_s}}{Numeric vector of estimated \eqn{P(S = 1 \mid X)} for
#'       the original data (length = \code{nrow} of original data).}
#'   }
#'
#' @keywords internal
trim_positivity_violations <- function(data, sample_col, outcome_col,
                                       treatment_col, threshold = 0.05) {
  covars <- setdiff(colnames(data), c(outcome_col, treatment_col, sample_col))
  if (length(covars) == 0L) {
    return(list(data = data, n_trimmed = 0L, alpha = 0, pi_s = rep(0.5, nrow(data))))
  }

  # Fit sampling model P(S = 1 | X)
  formula_s <- stats::as.formula(paste(sample_col, "~", paste(covars, collapse = "+")))
  model_s   <- stats::glm(formula_s, data = data, family = stats::binomial())
  pi_s      <- stats::predict(model_s, type = "response")

  # Determine the effective alpha

  if (identical(threshold, "crump")) {
    alpha <- crump_alpha(pi_s)
  } else {
    alpha <- as.numeric(threshold)
  }

  if (alpha <= 0) {
    return(list(data = data, n_trimmed = 0L, alpha = 0, pi_s = pi_s))
  }

  # Trim target rows with P(S=1|X) < alpha
  is_target <- data[[sample_col]] == 0
  drop_mask <- is_target & (pi_s < alpha)

  list(
    data      = data[!drop_mask, , drop = FALSE],
    n_trimmed = sum(drop_mask),
    alpha     = alpha,
    pi_s      = pi_s
  )
}


#' Summarize a ROOT fit
#'
#' Provides a readable summary of a \code{ROOT} object, including:
#' \enumerate{
#'   \item the summary characterization tree \code{f},
#'   \item whether the user supplied a custom \code{global_objective_fn} (Yes/No), and
#'   \item in generalizability mode (\code{generalizability_path = TRUE}), the
#'         unweighted and weighted estimands with their standard errors.
#' }
#'
#' @section Abbreviations:
#' \describe{
#'   \item{ATE}{Average treatment effect.}
#'   \item{RCT}{Randomized controlled trial.}
#'   \item{SE}{Standard error.}
#'   \item{TATE}{Target average treatment effect.}
#'   \item{WTATE}{Weighted target average treatment effect.}
#'   \item{STATE}{Sample target average treatment effect.}
#' }
#'
#' When \code{generalizability_path = TRUE}, the unweighted estimand corresponds
#' to a STATE-type quantity and the weighted estimand to a WTATE-type
#' quantity for the transported target population. When \code{generalizability_path = FALSE},
#' ROOT is used for general functional optimization and no causal labels
#' are imposed; the summary focuses on the tree and diagnostics.
#'
#' @section Diagnostics:
#' The summary also reports:
#' \itemize{
#'   \item the number of trees grown,
#'   \item the size of the Rashomon set,
#'   \item the percentage of observations with ensemble vote \code{w_opt == 1}.
#' }
#'
#' @param object A \code{"ROOT"} S3 object returned by \code{ROOT()}.
#' @param ... Currently unused and included for S3 compatibility.
#'
#' @return \code{object} returned invisibly. Printed output is for inspection.
#'
#' @method summary ROOT
#'
#' @examples
#' ROOT.output = ROOT(diabetes_data,generalizability_path = TRUE, seed = 123)
#' summary(ROOT.output)
#'
#' @export
summary.ROOT <- function(object, ...) {
  x <- object

  cat("ROOT object\n")
  cat("  Generalizability mode:", isTRUE(x$generalizability_path), "\n\n")

  ## Summary tree
  cat("Summary classifier (f):\n")
  if (!is.null(x$f)) {
    print(x$f)
  } else {
    cat("  <no summary tree available>\n")
  }
  cat("\n")

  ## Testing data
  #cat("Testing data (head):\n")
  #if (!is.null(x$testing_data)) {
  #  print(utils::head(x$testing_data))
  #} else {
  #  cat("  <no testing_data stored>\n")
  #}
  #cat("\n")

  ## Objective function
  cat("Global objective function:\n")
  if (isTRUE(x$user_supplied_objective)) {
    cat("  User-supplied: Yes\n")
  } else {
    cat("  User-supplied: No (default objective used)\n")
  }
  cat("\n")

  ## Estimands (only in generalizability mode)
  if (isTRUE(x$generalizability_path) && !is.null(x$estimate)) {
    est <- x$estimate
    cat("Estimand summary (generalizability mode):\n")
    cat("  Unweighted ", est$estimand_unweighted,
        " = ", est$value_unweighted,
        ", SE = ", est$se_unweighted, "\n", sep = "")
    cat("  Weighted   ", est$estimand_weighted,
        " = ", est$value_weighted,
        ", SE = ", est$se_weighted, "\n", sep = "")
    cat("\n")
    #cat("  Note: ", est$se_weighted_note, "\n", sep = "")
  }

  ## Core diagnostics
  cat("Diagnostics:\n")
  # number of trees
  if (!is.null(x$w_forest)) {
    cat("  Number of trees grown: ", length(x$w_forest), "\n", sep = "")
  } else {
    cat("  Number of trees grown: <unknown>\n")
  }

  # Rashomon size
  if (!is.null(x$rashomon_set)) {
    cat("  Rashomon set size: ", length(x$rashomon_set), "\n", sep = "")
  } else {
    cat("  Rashomon set size: <unknown>\n")
  }

  # Split strategy
  if (!is.null(x$split_strategy)) {
    cat("  Split strategy: ", x$split_strategy, "\n", sep = "")
  }

  # Positivity trimming
  if (!is.null(x$positivity_trim_info)) {
    pti <- x$positivity_trim_info
    cat(sprintf("  Positivity trimming: %s (alpha = %g, %d target obs. removed)\n",
                pti$method, pti$alpha, pti$n_trimmed))
  }

  # % with w_opt == 1
  if (!is.null(x$D_rash) && "w_opt" %in% names(x$D_rash)) {
    pct_keep <- mean(x$D_rash$w_opt == 1, na.rm = TRUE) * 100
    cat(sprintf("  %% observations with w_opt == 1: %.1f%%\n", pct_keep))
  } else {
    cat("  % observations with w_opt == 1: <not available>\n")
  }
  cat("\n")

  invisible(x)
}

#' Print a ROOT fit
#'
#' Provides a human-readable brief summary of a \code{ROOT} object, including:
#' \enumerate{
#'   \item the summary characterization tree \code{f},
#'   \item in generalizability mode (\code{generalizability_path = TRUE}), the
#'         unweighted and weighted estimands with their standard errors
#'         and an explanatory note for the weighted standard error (SE).
#' }
#'
#' @section Abbreviations:
#' \describe{
#'   \item{ATE}{Average treatment effect.}
#'   \item{RCT}{Randomized controlled trial.}
#'   \item{SE}{Standard error.}
#'   \item{TATE}{Target average treatment effect.}
#'   \item{WTATE}{Weighted target average treatment effect.}
#'   \item{STATE}{Sample target average treatment effect.}
#' }
#'
#' When \code{generalizability_path = TRUE}, the unweighted estimand corresponds
#' to a STATE-type quantity and the weighted estimand to a WTATE-type
#' quantity for the transported target population. When \code{generalizability_path = FALSE},
#' ROOT is used for general functional optimization and no causal labels
#' are imposed.
#'
#' @param x A \code{"ROOT"} S3 object returned by \code{ROOT()}.
#' @param ... Currently unused and included for S3 compatibility.
#'
#' @return \code{object} returned invisibly. Printed output is for inspection.
#'
#' @method print ROOT
#'
#' @examples
#' ROOT.output = ROOT(diabetes_data,generalizability_path = TRUE, seed = 123)
#' print(ROOT.output)
#'
#' @export
print.ROOT <- function(x, ...) {
  cat("ROOT object\n")
  cat("  Generalizability mode:", isTRUE(x$generalizability_path), "\n\n")

  ## Summary tree
  cat("Summary classifier (f):\n")
  if (!is.null(x$f)) {
    print(x$f)
  } else {
    cat("  <no summary tree available>\n")
  }
  cat("\n")

  ## Estimands (only in generalizability mode)
  if (isTRUE(x$generalizability_path) && !is.null(x$estimate)) {
    est <- x$estimate
    cat("Estimand summary (generalizability mode):\n")
    cat("  Unweighted ", est$estimand_unweighted,
        " = ", est$value_unweighted,
        ", SE = ", est$se_unweighted, "\n", sep = "")
    cat("  Weighted   ", est$estimand_weighted,
        " = ", est$value_weighted,
        ", SE = ", est$se_weighted, "\n", sep = "")
    #cat("  Note: ", est$se_weighted_note, "\n", sep = "")
  }

  invisible(x)
}

#' Plot the ROOT summary tree
#'
#' Visualizes the decision tree that characterizes the weighted subgroup
#' (the weight function \eqn{w(d)} in \code{\{0,1\}}) identified by \code{ROOT()},
#' using \code{rpart.plot::prp()}.
#'
#' @param x A \code{"ROOT"} S3 object returned by \code{ROOT()} with
#'   \code{x$f} an \code{rpart} object representing the summary / characterization tree.
#' @param ... Additional arguments passed to \code{rpart.plot::prp()}.
#'
#' @return No return value; the plot is drawn to the active graphics device.
#' @importFrom rpart.plot prp
#'
#' @examples
#' ROOT.output = ROOT(diabetes_data,generalizability_path = TRUE, seed = 123)
#' plot(ROOT.output)
#'
#' @export
plot.ROOT <- function(x, ...) {
  if (is.null(x$f) || !inherits(x$f, "rpart")) {
    message("No summary tree available to plot (x$f is NULL or not an 'rpart' object).")
    return(invisible(NULL))
  }

  # Default arguments (can be overridden by ...)
  args <- list(...)
  if (!"type" %in% names(args))          args$type <- 2
  if (!"extra" %in% names(args))         args$extra <- 109
  if (!"under" %in% names(args))         args$under <- TRUE
  if (!"faclen" %in% names(args))        args$faclen <- 0
  if (!"tweak" %in% names(args))         args$tweak <- 1.1
  if (!"fallen.leaves" %in% names(args)) args$fallen.leaves <- TRUE
  if (!"box.palette" %in% names(args))   args$box.palette <- c("pink", "palegreen3")
  if (!"shadow.col" %in% names(args))    args$shadow.col <- "gray"
  if (!"branch.lty" %in% names(args))    args$branch.lty <- 3
  if (!"main" %in% names(args))          args$main <- "Final Characterized Tree from Rashomon Set"

  args$x <- x$f
  do.call(rpart.plot::prp, args)
}
