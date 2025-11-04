#' Summarize a characterizing_underrep fit
#'
#' Prints the \code{ROOT} summary (un/weighted estimates with standard deviations) and a brief
#' overview of terminal rules from the annotated summary tree, if available.
#'
#' @param object A \code{characterizing_underrep} object.
#' @param ... Unused; included for S3 compatibility.
#' @return \code{object}, invisibly.
#' @details Delegates core statistics to \code{summary(object$root)}; previews up to
#'   five terminal rules when a summary tree exists, and reports plot availability.
#' @method summary characterizing_underrep
#' @export
summary.characterizing_underrep <- function(object, ...) {
  if (!inherits(object, "characterizing_underrep")) stop("Not a characterizing_underrep object.")
  cat("characterizing_underrep object\n")
  cat("  --- ROOT summary ---\n")
  summary(object$root)  # uses your summary.ROOT()

  # Leaf summary
  if (is.null(object$leaf_summary)) {
    cat("  Leaf summary:    none (no summarized tree)\n")
  } else {
    cat("  Leaf summary:    ", nrow(object$leaf_summary), " terminal nodes\n", sep = "")
    prev <- utils::head(object$leaf_summary, 10)
    prev$rule <- substr(prev$rule, 1, 60)
    print(prev, row.names = FALSE)
    if (nrow(object$leaf_summary) > 5) cat("  ...\n")
  }

  # Plots
  cat("  Plots available:\n")
  cat("    * ROOT tree:          ", if (!is.null(object$tree_plot_root))     "yes" else "no", "\n", sep = "")
  cat("    * Underrep tree:      ", if (!is.null(object$tree_plot_underrep)) "yes" else "no", "\n", sep = "")

  invisible(object)
}
