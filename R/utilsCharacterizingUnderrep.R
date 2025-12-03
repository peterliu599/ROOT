#' Summarize a characterizing_underrep fit
#'
#' Summarizes the \code{ROOT} summary which includes unweighted and weighted estimates
#' with standard errors. The weighted SE is omitted when a custom \code{global_objective_fn}
#' was used in \code{ROOT()}. Provides a brief overview of terminal rules from the
#' annotated summary tree when available.
#'
#' @section Abbreviations:
#' ATE means Average Treatment Effect. RCT means Randomized Controlled Trial.
#' SE means Standard Error. TATE means Transported ATE. WTATE means Weighted TATE.
#' WATE means Weighted ATE. PATE means Population ATE.
#'
#' @param object A \code{characterizing_underrep} S3 object. Expected components include
#'   \code{root} which is a \code{ROOT} object and may contain \code{f} which is an
#'   \code{rpart} object for the summary tree, and \code{leaf_summary} which is a
#'   \code{data.frame} with one row per terminal node and may include a \code{rule}
#'   column of type \code{character}.
#' @param ... Currently unused. Included for S3 compatibility.
#'
#' @return \code{object} returned invisibly. Printed output is a human readable summary.
#'
#' @details Delegates core statistics to \code{summary(object$root)}. Previews up to
#'   ten terminal rules when a summary tree exists and reports plot availability.
#'
#' @method summary characterizing_underrep
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
#'
#' # View Summary
#' summary(res)
#' }
#' @export
summary.characterizing_underrep <- function(object, ...) {
  if (!inherits(object, "characterizing_underrep")) stop("Not a characterizing_underrep object.")

  cat("characterizing_underrep object\n")
  cat("  --- ROOT summary ---\n")
  summary(object$root)  # uses summary.ROOT() defined above

  # Leaf summary
  if (is.null(object$leaf_summary)) {
    cat("  Leaf summary:    none (no summarized tree)\n")
  } else {
    cat("  Leaf summary:    ", nrow(object$leaf_summary), " terminal nodes\n", sep = "")
    prev <- utils::head(object$leaf_summary, 10)
    # Trim very long rule strings for console
    if ("rule" %in% names(prev)) prev$rule <- substr(prev$rule, 1, 60)
    print(prev, row.names = FALSE)
    if (nrow(object$leaf_summary) > 10) cat("  ...\n")
  }
  return(invisible(object))
}



#' Print a characterizing_underrep fit
#'
#' Prints the \code{ROOT} summary which includes unweighted and weighted estimates
#' with standard errors. The weighted SE is omitted when a custom \code{global_objective_fn}
#' was used in \code{ROOT()}. Provides a shorter overview of results.
#'
#' @section Abbreviations:
#' ATE means Average Treatment Effect. RCT means Randomized Controlled Trial.
#' SE means Standard Error. TATE means Transported ATE. WTATE means Weighted TATE.
#' WATE means Weighted ATE. PATE means Population ATE.
#'
#' @param x A \code{characterizing_underrep} S3 object. Expected components include
#'   \code{root} which is a \code{ROOT} object. The \code{ROOT} object may contain
#'   \code{D_forest} which is data frame like, \code{D_rash} which may include \code{w_opt},
#'   and \code{rashomon_set} which is an \code{integer} vector of selected tree indices.
#' @param ... Currently unused. Included for S3 compatibility.
#'
#' @return \code{x} returned invisibly. Printed output is a human readable summary.
#'
#' @details Delegates core statistics, number of trees grown and Rashomon size,
#'   and the percentage of observations with ensemble vote \code{w_opt == 1} to
#'   \code{print(x$root)}.
#'
#' @method print characterizing_underrep
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
#'
#' # Print core results
#' print(res)
#' }
#' @export
print.characterizing_underrep <- function(x, ...) {
  if (!inherits(x, "characterizing_underrep")) stop("Not a characterizing_underrep object.")

  cat("characterizing_underrep object\n")
  cat("  --- ROOT summary ---\n")
  invisible(print(x$root))  # uses print.ROOT()
}



#' Plot Under represented Population Characterization
#'
#' Visualizes the decision tree derived from the \code{ROOT} analysis. Highlights
#' which subgroups are represented where \code{w = 1} versus underrepresented where \code{w = 0}.
#'
#' @section Abbreviations:
#' ATE means Average Treatment Effect. RCT means Randomized Controlled Trial.
#' SE means Standard Error. TATE means Transported ATE. WTATE means Weighted TATE.
#' WATE means Weighted ATE. PATE means Population ATE.
#'
#' @param x A \code{characterizing_underrep} S3 object with \code{x$root$f} present as an
#'   \code{rpart} object for the summary or characterization tree. The \code{rpart} model
#'   should have classification levels that allow identification of the represented class.
#' @param ... Additional arguments passed to \code{rpart.plot::prp()}.
#'
#' @return \code{NULL}. The plot is drawn to the active graphics device.
#'
#' @importFrom rpart.plot prp
#' @importFrom graphics par legend
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
#'
#' # Plot the annotated tree
#' plot(res)
#' }
#' @export
plot.characterizing_underrep <- function(x, ...) {
  # --- 1. Safety Checks ---
  if (is.null(x$root$f)) {
    message("No summary tree available to plot (possibly no covariates or tree failed to grow).")
    return(invisible(NULL))
  }

  f   <- x$root$f
  frm <- f$frame
  is_leaf <- frm$var == "<leaf>"

  # --- 2. Identify the 'Represented' (w=1) Class Index ---
  # We need to know if 'yval=1' means Represented or if 'yval=2' means Represented.
  ylevels <- attr(f, "ylevels")
  if (is.null(ylevels)) ylevels <- f$ylevels

  pos_class_idx <- 2 # Default assumption: 2nd level is '1' (e.g. levels are "0", "1")

  if (!is.null(ylevels)) {
    # Robustly find the level named "1", "w=1", "yes", etc.
    # This handles cases where levels might be c("1", "0") or c("No", "Yes")
    hits <- which(tolower(ylevels) %in% c("1", "w=1", "yes", "represented", "keep", "true"))
    if (length(hits) > 0) {
      pos_class_idx <- hits[1]
    } else {
      # Fallback: use the last level (standard R behavior for 0 vs 1)
      pos_class_idx <- length(ylevels)
    }
  }

  # --- 3. Determine Status using rpart's predicted class (yval) ---
  # This guarantees the color matches the tree's actual decision (0 or 1).
  if (!is.null(ylevels)) {
    # Classification: Node is represented if its predicted class index matches the positive index
    is_represented <- (frm$yval == pos_class_idx)
  } else {
    # Regression: Node is represented if prediction >= 0.5
    is_represented <- (frm$yval >= 0.5)
  }

  # --- 4. Define Colors ---
  col_keep <- "#4E79A7" # Blue   (Represented)
  col_drop <- "#F28E2B" # Orange (Underrepresented)

  # Create a color vector for all nodes (initially white)
  box_col <- rep("white", nrow(frm))

  # Color only the leaves based on status
  box_col[is_leaf] <- ifelse(is_represented[is_leaf], col_keep, col_drop)

  # --- 5. Custom Labeling Function ---
  total_n <- frm$n[1]

  node_fun <- function(x, labs, digits, varlen) {
    rows <- x$frame
    out  <- character(nrow(rows))
    pct  <- if (total_n > 0) rows$n / total_n else 0

    for (i in seq_len(nrow(rows))) {
      if (rows$var[i] == "<leaf>") {
        # Re-check representation status for this specific node
        node_yval <- rows$yval[i]

        # Determine if this node is "Represented"
        node_is_rep <- FALSE
        if (!is.null(ylevels)) {
          node_is_rep <- (node_yval == pos_class_idx)
        } else {
          node_is_rep <- (node_yval >= 0.5)
        }

        pct_str <- sprintf("%.0f%%", 100 * pct[i])

        # Apply the requested labels
        if (node_is_rep) {
          out[i] <- paste0("REPRESENTED\n", pct_str)
        } else {
          out[i] <- paste0("UNDERREPRESENTED\n", pct_str)
        }
      } else {
        # Internal nodes get standard split labels
        out[i] <- labs[i]
      }
    }
    out
  }

  # --- 6. Plot Arguments ---
  args <- list(...)
  # Set defaults if not provided by user
  if (is.null(args$type))          args$type <- 2
  if (is.null(args$extra))         args$extra <- 0
  if (is.null(args$under))         args$under <- TRUE
  if (is.null(args$tweak))         args$tweak <- 1.1
  if (is.null(args$fallen.leaves)) args$fallen.leaves <- TRUE
  if (is.null(args$shadow.col))    args$shadow.col <- "gray"

  # Force our calculated values
  args$x <- f
  args$box.col <- box_col
  args$node.fun <- node_fun
  args$split.box.col <- NA # Transparent background for split labels

  # Draw the plot
  do.call(rpart.plot::prp, args)

  # Legend
  op <- graphics::par(xpd = NA); on.exit(graphics::par(op), add = TRUE)
  graphics::legend(
    "topleft",
    legend = c("w(x) = 1 Represented", "w(x) = 0 Underrepresented"),
    fill   = c(col_keep, col_drop),
    border = NA, bty = "n", cex = 0.9
  )
}
