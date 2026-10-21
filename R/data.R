#' Simulated diabetes trial and target data
#'
#' A simulated dataset for illustrating ROOT in generalizability mode. It
#' mimics a randomized trial of a glucose-lowering drug (about 2,000
#' participants, \code{S = 1}) and a target population of adults with type 2
#' diabetes (about 8,000 people, \code{S = 0}). Treatment and outcome are
#' observed only in the trial.
#'
#' The design contains all four combinations of trial under-sampling and
#' treatment-effect modification: chronic kidney disease (\code{CKD}), heart
#' failure (\code{HeartFailure}), insulin use (\code{Insulin}), and age 70 or
#' older are under-sampled \emph{and} modify the treatment effect; \code{Race_Black} and \code{Rural}
#' are under-sampled but do \emph{not} modify the effect; baseline
#' \code{HbA1c} mildly modifies the effect but is well sampled;
#' \code{Sex_Male}, \code{BMI}, and \code{DietYes} neither modify the effect
#' nor are under-sampled.
#'
#' @section Data-generating process:
#' Trial participation follows
#' \deqn{\mathrm{logit}\, P(S = 1 \mid X) = -0.9 - 1.25\,\mathrm{CKD} -
#'   0.8\,\mathrm{HF} - 0.4\,\mathrm{Insulin} - 0.018(\mathrm{Age} - 58) -
#'   0.5\,\mathrm{Black} - 0.35\,\mathrm{Rural},}
#' where HF denotes \code{HeartFailure}. Within the trial, \code{Tr} is
#' randomized with probability 0.5. The conditional average treatment effect
#' is
#' \deqn{\tau(X) = 1 + 0.1(\mathrm{HbA1c} - 8) - 2.0\,\mathrm{CKD} -
#'   2.4\,\mathrm{HF} + 1.5\,\mathrm{Insulin} +
#'   2.0\,\mathbf{1}\{\mathrm{Age} \ge 70\},}
#' and the outcome is \eqn{Y = \mu_0(X) + \mathrm{Tr}\,\tau(X) + \varepsilon},
#' with \eqn{\mu_0(X) = 0.2(\mathrm{HbA1c} - 8) + 0.1\,\mathrm{Insulin} -
#' 0.2\,\mathrm{Black} + 0.1\,\mathrm{Rural} - 0.01(\mathrm{BMI} - 31) +
#' 0.2\,\mathrm{DietYes} + 0.1\,\mathrm{HF}} and
#' \eqn{\varepsilon \sim N(0, 0.5^2)}. The effects are illustrative.
#' The full simulation code is in \code{data-raw/diabetes_data.R} in the
#' package source.
#'
#' @format A \code{data.frame} with 10,000 rows and 13 columns:
#' \describe{
#'   \item{Age}{Age in years (30 to 85).}
#'   \item{Sex_Male}{Indicator in \code{0}/\code{1}: male.}
#'   \item{Race_Black}{Indicator in \code{0}/\code{1}: race is Black.}
#'   \item{Rural}{Indicator in \code{0}/\code{1}: lives in a rural area.}
#'   \item{BMI}{Body mass index (kg/m^2).}
#'   \item{HbA1c}{Baseline HbA1c (percent).}
#'   \item{Insulin}{Indicator in \code{0}/\code{1}: currently uses insulin.}
#'   \item{CKD}{Indicator in \code{0}/\code{1}: chronic kidney disease.}
#'   \item{HeartFailure}{Indicator in \code{0}/\code{1}: heart failure.}
#'   \item{DietYes}{Indicator in \code{0}/\code{1}: enrolled in a diet program.}
#'   \item{S}{Sample indicator: \code{1} = trial, \code{0} = target.}
#'   \item{Tr}{Treatment assignment in \code{0}/\code{1} (trial only; \code{NA} in the target).}
#'   \item{Y}{Reduction in HbA1c (percentage points) at 6 months (trial only; \code{NA} in the target).}
#' }
#'
#' @usage data(diabetes_data)
#' @name diabetes_data
#' @docType data
#' @keywords datasets
#' @aliases diabetes_data
NULL
