#' Pool estimates and variances obtained by analysing multiple synthetic datasets
#'
#' This function pools estimates and variances which have been obtained by
#' analysing multiple synthetic imputations (e.g. created used [gFormulaImpute])
#' using the method developed by Raghunathan et al 2003.
#'
#' The only argument to `syntheticPool` is a set of model fits obtained by running
#' an analysis on an imputed dataset collection of class `mids`, as created
#' for example using the `mice` function in the `mice` package.
#'
#' The function returns a table containing the overall parameter estimates, the within, between and total imputation
#' variances, 95% confidence intervals, and p-values testing the null hypothesis
#' that the corresponding parameters equal zero.
#'
#' It is possible for the variance estimator developed by Raghunathan et al 2003 to
#' be negative. In this case `syntheticPool` stops and informs you to re-impute
#' using a larger number of imputations `M` and/or `nSim`.
#'
#' The development of the `gFormulaMI` package was supported by a grant from the UK
#' Medical Research Council (MR/T023953/1).
#'
#' @param fits Collection of model fits produced by a call of the form
#'  `with(imps, lm(y~regime))` where `imps` is a collection of imputed datasets
#'  of class `mids`.
#'
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @references Raghunathan TE, Reiter JP, Rubin DB. 2003. Multiple imputation for statistical
#'  disclosure limitation. Journal of Official Statistics, 19(1), p.1-16.
#'
#' @export
#'
#' @examples
#' set.seed(7626)
#' #impute synthetic datasets under two regimes of interest using gFormulaImpute
#' imps <- gFormulaImpute(data=simDataFullyObs,M=10,
#'                         trtVars=c("a0","a1","a2"),
#'                         trtRegimes=list(c(0,0,0),c(1,1,1)))
#' #fit linear model to final outcome with regime as covariate
#' fits <- with(imps, lm(y~factor(regime)))
#' #pool results using Raghunathan et al 2003 rules
#' syntheticPool(fits)
#'
#' @returns A matrix containing the pooled results.
syntheticPool <- function(fits) {

  M <- length(fits$analyses)
  p <- length(stats::coefficients(fits$analyses[[1]]))
  ests <- array(0, dim=c(M,p))
  within_vars <- array(0, dim=c(M,p))

  for (i in 1:M) {
    ests[i,] <- stats::coefficients(fits$analyses[[i]])
    within_vars[i,] <- diag(stats::vcov(fits$analyses[[i]]))
  }

  overall_ests <- colMeans(ests)
  v <- colMeans(within_vars)
  b <- diag(stats::var(ests))
  total_var <- (1+1/M)*b - v

  #check that all variances are positive
  if (any(total_var<=0)) {
    stop("Some parameters have estimated total variances which are not positive. Re-run the imputation process using larger M and/or nSim values.")
  }

  #degrees of freedom
  df <- (M-1)*(1-(M*v)/((M+1)*b))^2

  resTable <- array(0, dim=c(p,8))
  resTable[,1] <- overall_ests
  resTable[,2] <- v
  resTable[,3] <- b
  resTable[,4] <- total_var
  resTable[,5] <- df
  #95% confidence interval
  resTable[,6] <- overall_ests - stats::qt(0.975,df=df)*sqrt(total_var)
  resTable[,7] <- overall_ests + stats::qt(0.975,df=df)*sqrt(total_var)
  #two sided p-value
  resTable[,8] <- 2*stats::pt(-abs(overall_ests/sqrt(total_var)),df=df)

  colnames(resTable) <- c("Estimate", "Within", "Between", "Total", "df",
                          "95% CI L", "95% CI U", "p")
  rownames(resTable) <- names(stats::coefficients(fits$analyses[[1]]))
  resTable
}
