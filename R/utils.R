#' Check for missing values in training data
#'
#' Ensures there are no \code{NA} values in any of the relevant columns of the dataset.
#'
#' @param data A \code{data.frame}.
#' @param cols A \code{character} vector of column names to check in \code{data}.
#'
#' @return Invisibly returns \code{TRUE} if no \code{NA} is found. Otherwise an error is thrown.
#' @keywords internal
#' @noRd
check_no_na <- function(data, cols) {
  data_name <- deparse(substitute(data))
  for (col in cols) {
    if (anyNA(data[[col]])) {
      stop(sprintf("Data `%s` column '%s' contains missing values. Please handle NA before training.", data_name, col),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Generic objective interface
#'
#' Default objective for general optimization mode: the standard error of a
#' weighted mean of per-unit scores.
#'
#' Computes \code{sqrt(sum(w_i * vsq_i) / (sum(w_i))^2)}. If \code{D} has a
#' \code{vsq} column it is used as supplied. Otherwise, if \code{D} has a
#' \code{v} column, \code{vsq_i = (v_i - vbar_w)^2} is recomputed at every
#' call, where \code{vbar_w} is the mean of \code{v} over the currently
#' included units, so the criterion is a genuine (w-dependent) variance.
#' The goal is to minimize the value. The generalizability path uses
#' \code{objective_generalizability} instead.
#' Supply your own \code{function(D) -> numeric} to use a different objective.
#'
#' @param D A \code{data.frame} with at least numeric columns \code{vsq} and \code{w}.
#'
#' @return A \code{numeric(1)} objective value. Returns \code{Inf} when undefined.
#' @keywords internal
#' @noRd
objective_default <- function(D) {
  if (!("w" %in% names(D))) {
    stop("objective_default() expects a column `w` in D.", call. = FALSE)
  }
  w <- D$w
  if (all(is.na(w)) || sum(w, na.rm = TRUE) <= 0) return(Inf)

  if ("vsq" %in% names(D) && is.numeric(D$vsq) && any(is.finite(D$vsq))) {
    vsq <- D$vsq
  } else if ("v" %in% names(D) && is.numeric(D$v) && any(is.finite(D$v))) {
    v  <- D$v
    mu <- stats::weighted.mean(v, w = w, na.rm = TRUE)
    vsq <- (v - mu)^2
  } else {
    # Fallback when no variance proxy is available:
    # Standard error (SE) of mean of a Bernoulli with p approx mean(w>0) and n = #kept.
    n_keep <- sum(w > 0, na.rm = TRUE)
    if (n_keep <= 1) return(Inf)
    p <- mean(w > 0, na.rm = TRUE)
    return(sqrt(p * (1 - p) / n_keep))
  }

  num <- sum(w * vsq, na.rm = TRUE)
  den <- sum(w, na.rm = TRUE)^2
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(Inf)
  sqrt(num / den)
}

#' Default objective for generalizability mode: SE of the IPW WTATE estimator
#'
#' Finite-sample standard error of the IPW estimator of the weighted target
#' average treatment effect (Parikh et al. 2025; Appendix A of the ROOT paper),
#' \deqn{L(w) = \kappa_w \sqrt{\sum_{i \in A} (v_i - \bar v_A)^2} / n_A,}
#' where \eqn{A = \{i : S_i = 1, w_i = 1\}}, \eqn{n_A = |A|}, \eqn{\bar v_A}
#' is the mean of \code{v} over \eqn{A}, and
#' \deqn{\kappa_w = \frac{n_A / \sum_{i \in A} r_i}{n_1 / \sum_{S_i = 1} r_i}}
#' re-normalizes the selection odds from the full population,
#' \eqn{\pi/(1-\pi)}, to the refined population,
#' \eqn{\pi_w/(1-\pi_w)} with \eqn{\pi_w = P(S = 1 \mid w(X) = 1)}.
#' Here \eqn{r_i = (1-\hat\pi(X_i))/\hat\pi(X_i)} (column \code{sel_odds}),
#' and \eqn{\sum_{i \in A} r_i} estimates the number of target units with
#' \eqn{w = 1}, because \eqn{E[S\,w(X)\,r(X)] = E[(1-S)\,w(X)]}.
#' \eqn{\kappa_w = 1} when \eqn{w \equiv 1}.
#' This equals \eqn{\sqrt{\hat V^{ipw,w}/n_1}}.
#'
#' @param D A \code{data.frame} with numeric columns \code{w}, \code{v}
#'   (from \code{compute_transport_scores}, trial rows only) and \code{sel_odds}.
#' @return A \code{numeric(1)} objective value; \code{Inf} when undefined.
#' @keywords internal
#' @noRd
objective_generalizability <- function(D) {
  if (!all(c("w", "v", "sel_odds") %in% names(D))) {
    stop("objective_generalizability() expects columns `w`, `v`, and `sel_odds` in D.",
         call. = FALSE)
  }
  ok <- !is.na(D$v) & is.finite(D$v) & !is.na(D$sel_odds)
  w  <- D$w
  A  <- ok & !is.na(w) & w == 1
  n_A <- sum(A)
  if (n_A < 2L) return(Inf)
  kappa <- wtate_kappa(D$sel_odds, A, ok)
  v_A <- D$v[A]
  val <- kappa * sqrt(sum((v_A - mean(v_A))^2)) / n_A
  if (!is.finite(val)) return(Inf)
  val
}

#' Selection-odds re-normalization for a refined population
#'
#' Returns \eqn{\kappa_w = (n_A / \sum_{A} r_i) / (n_1 / \sum_{trial} r_i)},
#' the factor converting the full-population IPW scores \code{v} to the
#' refined population \eqn{w = 1}.
#'
#' @param sel_odds Numeric vector of \eqn{r_i} for trial units.
#' @param A Logical vector: units with \eqn{w = 1}.
#' @param ok Logical vector: usable trial units.
#' @return \code{numeric(1)}.
#' @keywords internal
#' @noRd
wtate_kappa <- function(sel_odds, A, ok = rep(TRUE, length(sel_odds))) {
  (sum(A) / sum(sel_odds[A])) / (sum(ok) / sum(sel_odds[ok]))
}

#' Helper to evaluate the objective after a hypothetical local change
#'
#' Evaluates \code{global_objective_fn} on a temporary copy of \code{D} after setting
#' \code{w = val} for the rows selected by \code{indices}.
#'
#' @param val A \code{numeric(1)} that must be either \code{0} or \code{1}.
#' @param indices An \code{integer} vector of row indices or a \code{character} vector of row names that receive \code{val}.
#' @param D A \code{data.frame} used by \code{global_objective_fn}. Must contain columns \code{w} and \code{vsq}.
#' @param global_objective_fn A \code{function} with signature \code{function(D) -> numeric}.
#'
#' @return A \code{numeric(1)} objective value after the hypothetical change.
#' @keywords internal
#' @noRd
objective_if <- function(val, indices, D, global_objective_fn) {
  stopifnot(is.function(global_objective_fn), length(val) == 1, val %in% c(0,1))
  rows <- integer(0)
  if (length(indices)) {
    rows <- if (is.numeric(indices)) as.integer(indices)
    else which(rownames(D) %in% as.character(indices))
  }
  Dtmp <- D
  if (length(rows)) Dtmp[rows, "w"] <- val
  global_objective_fn(Dtmp)
}

#' Backward and fast path micro evaluator adaptor
#'
#' Wraps a global objective \code{global_objective_fn(D)} into a splitter compatible
#' loss function \code{loss_fn(val, indices, D)} by calling \code{\link{objective_if}}
#' on a temporary copy of \code{D}.
#'
#' @param global_objective_fn A \code{function} with signature \code{function(D) -> numeric}
#'   that returns a scalar to be minimized. For example \code{\link{objective_default}}.
#'
#' @return A \code{function} \code{loss_fn(val, indices, D)} suitable for use in
#'   \code{\link{ROOT}} and \code{\link{split_node}}. It sets \code{w = val} on \code{indices}
#'   without mutation of the original \code{D} and then returns \code{global_objective_fn(D)}.
#' @keywords internal
#' @noRd
loss_from_objective <- function(global_objective_fn) {
  force(global_objective_fn)
  function(val, indices, D) objective_if(val, indices, D, global_objective_fn)
}
