#' Summarize a characterizing_underrep fit
#'
#' Summarizes the \code{ROOT} summary which includes unweighted and (when in
#' generalization mode) weighted estimates with standard errors, as reported by
#' \code{summary.ROOT()}. Provides a brief overview of terminal rules from the
#' annotated summary tree when available.
#'
#' @section Abbreviations:
#' ATE means Average Treatment Effect. RCT means Randomized Controlled Trial.
#' SE means Standard Error. TATE means Transported ATE. WTATE means Weighted TATE.
#' WATE means Weighted ATE. PATE means Population ATE.
#'
#' @param object A \code{characterizing_underrep} S3 object. Expected components include
#'   \code{root} which is a \code{ROOT} object (summarized by \code{summary.ROOT()})
#'   and may contain \code{f} which is an \code{rpart} object for the summary tree,
#'   and \code{leaf_summary} which is a \code{data.frame} with one row per terminal
#'   node and may include a \code{rule} column of type \code{character}.
#' @param ... Currently unused. Included for S3 compatibility.
#'
#' @return \code{object} returned invisibly. Printed output is a human readable summary.
#'
#' @details
#' Delegates core statistics and estimands to \code{summary(object$root)}.
#' Previews up to ten terminal rules when a summary tree exists.
#'
#' @method summary characterizing_underrep
#' @examples
#' \dontrun{
#' char.output = characterizing_underrep(diabetes_data,generalizability_path = TRUE, seed = 123)
#' summary(char.output)
#' }
#' @export
summary.characterizing_underrep <- function(object, ...) {
  if (!inherits(object, "characterizing_underrep")) {
    stop("Not a characterizing_underrep object.")
  }

  cat("characterizing_underrep object\n")
  cat("  --- ROOT summary ---\n")
  summary(object$root)  # uses summary.ROOT()

  # Leaf summary
  if (is.null(object$leaf_summary)) {
    cat("  Leaf summary:    none (no summarized tree)\n")
  } else {
    cat("  Leaf summary:    ", nrow(object$leaf_summary),
        " terminal nodes\n", sep = "")
    prev <- utils::head(object$leaf_summary, 10L)
    if ("rule" %in% names(prev)) {
      prev$rule <- substr(prev$rule, 1L, 60L)
    }
    print(prev, row.names = FALSE)
    if (nrow(object$leaf_summary) > 10L) cat("  ...\n")
  }
  invisible(object)
}


#' Print a characterizing_underrep fit
#'
#' Print the \code{ROOT} summary which includes unweighted and (when in
#' generalization mode) weighted estimates with standard errors, as reported by
#' \code{summary.ROOT()}.
#'
#' @section Abbreviations:
#' ATE means Average Treatment Effect. RCT means Randomized Controlled Trial.
#' SE means Standard Error. TATE means Transported ATE. WTATE means Weighted TATE.
#' WATE means Weighted ATE. PATE means Population ATE.
#'
#' @param x A \code{characterizing_underrep} S3 object. Expected components include
#'   \code{root} which is a \code{ROOT} object (summarized by \code{print.ROOT()})
#'   and may contain \code{f} which is an \code{rpart} object for the summary tree,
#'   and \code{leaf_summary} which is a \code{data.frame} with one row per terminal
#'   node and may include a \code{rule} column of type \code{character}.
#' @param ... Currently unused. Included for S3 compatibility.
#'
#' @return \code{object} returned invisibly. Printed output is a human readable brief summary.
#'
#' @details
#' Delegates core statistics and estimands to \code{print(x$root)}.
#'
#' @method print characterizing_underrep
#' @examples
#' \dontrun{
#' char.output = characterizing_underrep(diabetes_data,generalizability_path = TRUE, seed = 123)
#' print(char.output)
#' }
#' @export
print.characterizing_underrep <- function(x, ...) {
  if (!inherits(x, "characterizing_underrep")) {
    stop("Not a characterizing_underrep object.")
  }

  cat("characterizing_underrep object\n")
  cat("  --- ROOT brief summary ---\n")
  print(x$root)  # uses print.ROOT()

  invisible(x)
}


#' Plot Under represented Population Characterization
#'
#' Visualizes the decision tree derived from the \code{ROOT} analysis. Highlights
#' which subgroups are represented where \code{w = 1} versus underrepresented
#' where \code{w = 0} in generalization mode, or simply \code{w(x)} in \code{\{0,1\}}
#' in general optimization mode.
#'
#' @param x A \code{characterizing_underrep} S3 object with \code{x$root$f} present as an
#'   \code{rpart} object for the summary or characterization tree.
#' @param main Character string for the plot title. Default is
#'   \code{"Underrepresented Population Characterization"}.
#' @param cex.main Numeric scaling factor for the title text size. Default is \code{1.2}.
#' @param ... Additional arguments passed to \code{rpart.plot::prp()}.
#'
#' @return \code{NULL}. The plot is drawn to the active graphics device.
#'
#' @importFrom rpart.plot prp
#' @importFrom graphics par legend title
#' @export
plot.characterizing_underrep <- function(x,
                                         main = "Final Characterized Tree from Rashomon Set",
                                         cex.main = 1.2,
                                         ...) {
  # --- 1. Safety Checks ---
  if (is.null(x$root$f)) {
    message("No summary tree available to plot (possibly no covariates or tree failed to grow).")
    return(invisible(NULL))
  }

  f   <- x$root$f
  frm <- f$frame
  is_leaf <- frm$var == "<leaf>"

  # --- 2. Identify the 'Represented' (w=1) Class Index ---
  ylevels <- attr(f, "ylevels")
  if (is.null(ylevels)) ylevels <- f$ylevels

  pos_class_idx <- 2L # default assumption

  if (!is.null(ylevels)) {
    hits <- which(tolower(ylevels) %in% c("1", "w=1", "yes", "represented", "keep", "true"))
    if (length(hits) > 0L) {
      pos_class_idx <- hits[1L]
    } else {
      pos_class_idx <- length(ylevels)
    }
  }

  # --- 3. Determine status (represented vs not) from yval ---
  if (!is.null(ylevels)) {
    is_represented <- (frm$yval == pos_class_idx)
  } else {
    # regression fallback
    is_represented <- (frm$yval >= 0.5)
  }

  # --- 4. Define Colors ---
  col_keep <- "#4E79A7" # Blue   (w = 1)
  col_drop <- "#F28E2B" # Orange (w = 0)

  box_col <- rep("white", nrow(frm))
  box_col[is_leaf] <- ifelse(is_represented[is_leaf], col_keep, col_drop)

  # --- 5. Custom labeling: show % in each leaf ---
  total_n <- frm$n[1L]

  node_fun <- function(x, labs, digits, varlen) {
    rows <- x$frame
    out  <- character(nrow(rows))
    pct  <- if (total_n > 0) rows$n / total_n else 0

    for (i in seq_len(nrow(rows))) {
      if (rows$var[i] == "<leaf>") {
        pct_str <- sprintf("%.0f%%", 100 * pct[i])
        out[i] <- pct_str
      } else {
        out[i] <- labs[i]
      }
    }
    out
  }

  # --- 6. Plot arguments (bigger defaults) ---
  args <- list(...)
  # Remove main from args if present - we handle title separately
  args$main <- NULL

  if (is.null(args$type))          args$type <- 2
  if (is.null(args$extra))         args$extra <- 0
  if (is.null(args$under))         args$under <- TRUE
  if (is.null(args$tweak))         args$tweak <- 1.3
  if (is.null(args$fallen.leaves)) args$fallen.leaves <- TRUE
  if (is.null(args$shadow.col))    args$shadow.col <- "gray"

  args$x             <- f
  args$box.col       <- box_col
  args$node.fun      <- node_fun
  args$split.box.col <- NA

  # --- 7. Draw plot, title, legend (legend at bottom) ---
  op_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op_par), add = TRUE)

  # Add room for a bottom legend (tune 2.5 up/down if needed)
  graphics::par(mar = op_par$mar + c(2.5, 0, 0, 0))

  do.call(rpart.plot::prp, args)

  if (!is.null(main) && nzchar(main)) {
    suppressWarnings(graphics::title(main = main, cex.main = cex.main))
  }

  # Legend text depends on generalizability_path
  is_gen <- isTRUE(x$root$generalizability_path)
  legend_labels <- if (is_gen) {
    c("w(x) = 1 Sufficiently represented",
      "w(x) = 0 Underrepresented")
  } else {
    c("w(x) = 1",
      "w(x) = 0")
  }

  # Put legend below plot
  graphics::par(xpd = NA)
  suppressWarnings(graphics::legend(
    x      = "bottom",
    legend = legend_labels,
    fill   = c(col_keep, col_drop),
    border = NA,
    bty    = "n",
    horiz  = TRUE,
    cex    = 1.0,
    inset  = c(0, -0.04)  # move legend downward; adjust if needed
  ))

  invisible(NULL)
}
