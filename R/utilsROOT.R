#' Detect single-sample mode in a ROOT fit
#'
#' Returns \code{TRUE} if a \code{ROOT} object represents single-sample (trial-only)
#' mode (ATE in RCT/Weighted ATE in RCT), and \code{FALSE} otherwise (TATE/WTATE).
#'
#' @param object A \code{ROOT} object (list with component \code{D_forest}).
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.root_is_single_sample <- function(object) {
  # Prefer the explicit flag if present
  if (!is.null(object$single_sample_mode)) return(isTRUE(object$single_sample_mode))
  # Fallback: infer from lX only (do NOT check S==1, since testing_data is S==1 by design)
  lx <- if ("lX" %in% names(object$D_forest)) object$D_forest$lX else NA_real_
  all(is.na(lx))
}

#' Extract covariate names used by ROOT
#'
#' Returns the names of covariate columns in \code{object$D_forest}, excluding internal
#' columns such as \code{v}, \code{vsq}, \code{S}, \code{lX}, and any \code{w_tree_*}.
#'
#' @param object A \code{ROOT} object.
#' @return Character vector of column names.
#' @keywords internal
#' @noRd
.root_covariate_names <- function(object) {
  wt_cols <- grep("^w_tree_", names(object$D_forest), value = TRUE)
  setdiff(names(object$D_forest), c("v","vsq","S","lX", wt_cols))
}

#' Compute the baseline loss for a ROOT fit
#'
#' Computes the baseline loss (no selection) as \eqn{\sqrt{\sum vsq}/n^2} using
#' \code{object$D_forest$vsq} and the number of rows \code{n}.
#'
#' @param object A \code{ROOT} object.
#' @return Numeric scalar (baseline loss).
#' @keywords internal
#' @noRd
.root_baseline_loss <- function(object) {
  n <- nrow(object$D_forest)
  sqrt(sum(object$D_forest$vsq, na.rm = TRUE) / (n^2))
}

#' Get objectives of the selected Rashomon trees
#'
#' Returns the vector of objective values for trees included in the Rashomon set.
#' If none are selected, returns \code{numeric(0)}.
#'
#' @param object A \code{ROOT} object.
#' @return Numeric vector of objectives (possibly empty).
#' @keywords internal
#' @noRd
.root_selected_objectives <- function(object) {
  vals <- vapply(object$w_forest, function(t) { val <- unname(t[["local objective"]]); if (is.null(val)) NA_real_ else val },
                 numeric(1))
  if (length(object$rashomon_set) == 0) return(numeric(0))
  vals[object$rashomon_set]
}


#' Summarize a ROOT fit
#'
#' Summarizes a \code{ROOT} object by reporting the primary estimands and key
#' model diagnostics. The first lines report:
#' \enumerate{
#'   \item the \strong{unweighted} estimate (ATE in RCT for single-sample or TATE for two-sample)
#'         and its \strong{standard error (SE)};
#'   \item the \strong{weighted} estimate (Weighted ATE in RCT or WTATE using \code{w_opt}) and its
#'         \strong{SE} \emph{when the default objective is used}; if a custom
#'         \code{global_objective_fn} was supplied in \code{ROOT()}, the weighted SE is omitted.
#' }
#' Subsequent lines describe the estimand type, number of trees, size of the Rashomon set,
#' presence of a summary tree, covariate count, observation count, baseline loss,
#' selected-tree losses, and the proportion kept by \code{w_opt}.
#'
#' @param object A \code{ROOT} object returned by \code{ROOT()}.
#' @param ... Unused; included for S3 compatibility.
#'
#' @return The input \code{object}, invisibly. Printed output is a human-readable summary.
#'
#' @details
#' This method prefers pre-computed estimates in \code{object$estimate} (as stored by \code{ROOT()}).
#' If unavailable, it recomputes:
#' \itemize{
#'   \item Unweighted effect as \eqn{\bar v}{mean(v)} over the analysis set
#'         (all rows in single-sample; \code{S == 1} in two-sample);
#'   \item Unweighted SE as \eqn{\sqrt{\frac{1}{n}\sum_{i=1}^n (v_i - \bar v)^2}}{sqrt(var(v)/n)}
#'         on the same set;
#'   \item Weighted effect as \eqn{\frac{\sum_i w_i v_i}{\sum_i w_i}}{sum(w*v)/sum(w)},
#'         where \code{w = w\_opt} (binary);
#'   \item Weighted SE as \eqn{\frac{\sqrt{\sum_i w_i (v_i - \bar v_w)^2}}{\sum_i w_i}}{sqrt(sum(w*(v-mu_w)^2))/sum(w)},
#'         printed only when the default objective was used.
#' }
#'
#' @method summary ROOT
#' @export
summary.ROOT <- function(object, ...) {
  if (!inherits(object, "ROOT")) stop("Not a ROOT object.")

  # Analysis set (single-sample: all rows; two-sample: S==1)
  single <- .root_is_single_sample(object)
  in_S <- if (single) rep(TRUE, nrow(object$D_forest)) else (object$D_forest$S == 1L)

  # Prefer precomputed estimates from ROOT()
  if (!is.null(object$estimate)) {
    eu <- object$estimate

    # Unweighted line
    if (!is.null(eu$estimand_unweighted) && !is.null(eu$value_unweighted) && !is.null(eu$se_unweighted)) {
      cat(sprintf("%s (unweighted) = %.6f, SE = %.6f\n",
                  eu$estimand_unweighted, eu$value_unweighted, eu$se_unweighted))
    }

    # Weighted line with conditional SE
    if (!is.null(eu$estimand_weighted) && !is.null(eu$value_weighted)) {
      if (!is.null(eu$se_weighted) && is.finite(eu$se_weighted)) {
        cat(sprintf("%s (weighted)   = %.6f, SE = %.6f\n",
                    eu$estimand_weighted, eu$value_weighted, eu$se_weighted))
      } else {
        note <- if (!is.null(eu$se_weighted_note)) paste0(" (", eu$se_weighted_note, ")") else ""
        cat(sprintf("%s (weighted)   = %.6f%s\n",
                    eu$estimand_weighted, eu$value_weighted, note))
      }
    }

  } else {
    # Fallback recompute (SE-only); respect whether default objective was used when possible
    label_unw <- if (single) "ATE in RCT" else "TATE"
    label_w   <- if (single) "Weighted ATE in RCT" else "WTATE"

    v <- object$D_forest$v[in_S]
    w <- if ("w_opt" %in% names(object$D_rash)) object$D_rash$w_opt[in_S] else rep(1L, length(v))

    mu_unw <- mean(v, na.rm = TRUE)
    n_eff_unw <- sum(!is.na(v))
    se_unw <- if (n_eff_unw > 1) sqrt(stats::var(v, na.rm = TRUE) / n_eff_unw) else NA_real_
    cat(sprintf("%s (unweighted) = %.6f, SE = %.6f\n", label_unw, mu_unw, se_unw))

    den_w <- sum(w, na.rm = TRUE)
    if (isTRUE(den_w > 0)) {
      mu_w <- sum(w * v, na.rm = TRUE) / den_w
      is_default_obj <- if (!is.null(object$objective_is_default)) {
        isTRUE(object$objective_is_default)
      } else {
        TRUE  # assume default if flag absent (back-compat)
      }
      if (is_default_obj) {
        se_w <- sqrt(sum(w * (v - mu_w)^2, na.rm = TRUE)) / den_w
        cat(sprintf("%s (weighted)   = %.6f, SE = %.6f\n", label_w, mu_w, se_w))
      } else {
        cat(sprintf("%s (weighted)   = %.6f (SE omitted: custom global_objective_fn; please compute your own SE)\n",
                    label_w, mu_w))
      }
    } else {
      cat(sprintf("%s (weighted)   = NA (no kept observations)\n", label_w))
    }
  }

  # Body — diagnostics
  cat("ROOT object\n")
  estimand <- if (single) "ATE in RCT (single-sample)" else "TATE/PATE (transported ATE)"
  cat("  Estimand:       ", estimand, "\n", sep = "")

  wt_cols <- grep("^w_tree_", names(object$D_forest), value = TRUE)
  cat("  Trees grown:    ", length(wt_cols), "\n", sep = "")
  cat("  Rashomon size:  ", length(object$rashomon_set), "\n", sep = "")

  if (!is.null(object$f)) cat("  Summary tree:    present (rpart)\n") else cat("  Summary tree:    none\n")

  covs <- .root_covariate_names(object)
  cat("  Covariates:     ", length(covs), " (", paste(utils::head(covs, 4), collapse = ", "),
      if (length(covs) > 4) ", ..." else "", ")\n", sep = "")
  cat("  Observations:   ", nrow(object$testing_data), " (analysis set for trees)\n", sep = "")

  base <- .root_baseline_loss(object)
  sel  <- .root_selected_objectives(object)
  cat("  Baseline loss:  ", sprintf("%.5f", base), "\n", sep = "")
  if (length(sel)) {
    cat("  Selected loss:  min/median = ",
        sprintf("%.5f/%.5f", min(sel, na.rm = TRUE), stats::median(sel, na.rm = TRUE)), "\n", sep = "")
  } else {
    cat("  Selected loss:  (no trees selected)\n")
  }

  if ("w_opt" %in% names(object$D_rash)) {
    keep_rate <- mean(object$D_rash$w_opt[in_S] %in% c(1L, 1), na.rm = TRUE)
    cat("  Kept (w_opt=1): ", sprintf("%.1f%%", 100 * keep_rate), "\n", sep = "")
  }

  invisible(object)
}
