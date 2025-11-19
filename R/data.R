#' Simulated Diabetes Trial and Target Population Data
#'
#' A simulated dataset combining a randomized clinical trial (RCT) and a target population data.
#' Used to demonstrate the functionality of transporting treatment effects to underrepresented populations.
#'
#' @format A data frame with 2500 rows and 7 variables:
#' \describe{
#'   \item{Race_Black}{Binary covariate: 1 if Black, 0 otherwise.}
#'   \item{Sex_Male}{Binary covariate: 1 if Male, 0 otherwise.}
#'   \item{DietYes}{Binary covariate: 1 if managing diet, 0 otherwise.}
#'   \item{Age45}{Binary covariate: 1 if age >= 45, 0 otherwise.}
#'   \item{Y}{Observed outcome (blood sugar levels).}
#'   \item{Tr}{Treatment assignment (1 = Treatment 1, 0 = Treatment 2). NA for Target population.}
#'   \item{S}{Sample indicator: 1 if in RCT, 0 if in Target population.}
#' }
#' @usage data(diabetes_data)
"simulated_diabetes_data"
