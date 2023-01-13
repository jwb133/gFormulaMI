#' Simulated fully observed data frame
#'
#' A simulated observational study data frame with no missing data.
#'
#' @format ## `simDataFullyObs`
#' A data frame with 5000 rows and the following variables:
#' \describe{
#'   \item{l0}{Continuous baseline confounder}
#'   \item{a0}{Binary baseline treatment}
#'   \item{l1}{Continuous confounder at time 1}
#'   \item{a1}{Binary treatment at time 1}
#'   \item{l2}{Continuous confounder at time 2}
#'   \item{a2}{Binary treatment at time 2}
#'   \item{y}{Continuous final outcome}
#'   ...
#' }
"simDataFullyObs"
