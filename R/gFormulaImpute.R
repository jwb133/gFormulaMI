#' G-formula multiple imputation
#'
#' `gFormulaImpute` creates multiple imputed synthetic datasets of longitudinal histories under
#' specified treatment regimes of interest, based on the G-formula.
#'
#' `gFormulaImpute` creates multiple imputed synthetic datasets of longitudinal histories under
#' specified treatment regimes of interest, based on the G-formula, as described by Bartlett et al (2025)
#' \doi{10.1177/09622802251316971}. Specifically, to
#' the observed data frame, an additional `nSim` rows are added in which all variables are set to
#' missing, except the time-varying treatment variables. The latter are set to the values
#' as specified in the `trtRegimes` argument. If multiple treatment regimes are specified,
#' `nSim` rows are added for each of the specified treatment regimes.
#'
#' `gFormulaImpute` uses the `mice` package to impute the potential outcome values of the
#' time-varying confounders and outcome in the synthetic datasets. Imputation is performed
#' sequentially from left to right in the data frame. As such, the variables must be ordered
#' in time in the input data frame, with the time-varying confounders at each time followed
#' by the corresponding treatment variable at that time.
#'
#' For the data argument, `gFormulaImpute` expects either a fully observed (complete) data frame,
#' or else a set of multiple imputation stored in an object of class `mids` (from the `mice` package).
#'
#' Unlike with Rubin's regular multiple imputation pooling rules, it is possible
#' for the pooling rules developed by Raghunathan et al (2003) to give negative
#' variance estimates. The probability of this occurring is reduced by increasing
#' `M` and/or `nSim`.
#'
#' `gFormulaImpute` returns an object of class `mids`. This can be analysed using the same
#' methods that imputed datasets from `mice` can be analysed with (see examples). However, Rubin's standard
#' pooling rules are not valid for analysis of the synthetic datasets. Instead, the synthetic
#' variance estimator of Raghunathan et al (2003) must be used, as implemented in the
#' [syntheticPool] function.
#'
#' The development of the `gFormulaMI` package was supported by a grant from the UK
#' Medical Research Council (MR/T023953/1).
#'
#' @param data The observed data frame
#' @param M The number of imputed datasets to generate
#' @param trtVars A vector of variable names indicating the time-varying treatment variables
#' @param trtRegimes A vector specifying the treatment regime of interest, or a list of
#' vectors specifying the treatment regimes of interest
#' @param nSim The number of individuals to simulate in each imputed dataset. Defaults to
#' number of individuals in observed data
#' @param micePrintFlag TRUE/FALSE specifying whether the output from the call(s) to mice
#' should be printed
#' @param silent TRUE/FALSE indicating whether to print output to console (FALSE) or not (TRUE)
#' @param method An optional method argument to pass to mice. If specified, this should be
#' a named vector of methods, indicating which imputation method to use for each
#' time-varying confounder variable and final outcome. If not specified, the default
#' is to impute continuous variables using normal linear regression (norm), binary variables using
#' logistic regression (logreg), polytomous regression for unordered factors and
#' proportional odds model for ordered factors
#' @param predictorMatrix An optional predictor matrix to specify which variables to use
#' as predictors in the imputation models. The default is to impute sequentially, i.e. impute
#' using all variables to the left of the variable being imputed as covariates
#' @param missingDataCheck TRUE/FALSE indicating whether `gFormulaMI` checks, when
#' passed a regular data frame, whether there any missing values.
#'
#' @returns an S3 object of class mids (multiply imputed dataset)
#'
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @references Bartlett JW, Olarte Parra C, Granger E, Keogh RH., van Zwet EW
#' and Daniel RM, 2025. G-formula with multiple imputation for causal inference
#' with incomplete data. Statistical Methods in Medical Research.
#'
#' Raghunathan TE, Reiter JP, Rubin DB. 2003. Multiple imputation for statistical
#'  disclosure limitation. Journal of Official Statistics, 19(1), p.1-16.
#'
#' @export
#'
#' @examples
#' set.seed(7626)
#' #impute synthetic datasets under two regimes of interest
#' imps <- gFormulaImpute(data=simDataFullyObs,M=10,
#'                         trtVars=c("a0","a1","a2"),
#'                         trtRegimes=list(c(0,0,0),c(1,1,1)))
#' #fit linear model to final outcome with regime as covariate
#' fits <- with(imps, lm(y~factor(regime)))
#' #pool results using Raghunathan et al 2003 rules
#' syntheticPool(fits)
#'
#'
gFormulaImpute <- function(data, M=50, trtVars, trtRegimes,
                           nSim=NULL, micePrintFlag=FALSE,silent=FALSE,
                           method=NULL,predictorMatrix=NULL,
                           missingDataCheck=TRUE) {

  if (inherits(data, "mids")) {
    missingData <- TRUE
    if (silent==FALSE) {
      print("Input data is a mice created multiple imputation object.")
    }
    if (data$m!=M) {
      if (silent==FALSE) {
        print("Value passed to M being ignored.")
        print(paste("Number of synthetic imputations to be generated set to",data$m, "as in mids object passed to gFormulaImpute."))
      }
    }
    M <- data$m
    firstImp <- mice::complete(data,1)
  } else if (inherits(data,"data.frame"))  {
    missingData <- FALSE
    if (silent==FALSE) {
      print("Input data is a regular data frame.")
    }
    #check there are no missing values, unless user has turned off this check
    if (missingDataCheck==TRUE) {
      if (sum(is.na(data))>0) {
        stop("Missing values detected - please multiply impute these and pass a mids type object as input.")
      }
    }
  } else {
    stop("Input dataset should either be a data frame or a mids object created by mice.")
  }

  timePoints <- length(trtVars)

  if (typeof(trtRegimes)=="list") {
    numRegimes <- length(trtRegimes)
    #check trtRegimes elements each have length equal to timePoints
    for (i in 1:numRegimes) {
      if (length(trtRegimes[[i]])!=timePoints) {
        stop("Length of one of the specified regimes does not match number of treatment variables.")
      }
    }
  } else {
    numRegimes <- 1
    #check trtRegimes elements each have length equal to timePoints
    if (length(trtRegimes)!=timePoints) {
      stop("Length of treatment regime does not match number of treatment variables.")
    }
  }

  #check treatment variables are in data frame
  if (missingData==FALSE) {
    dataVars <- colnames(data)
  } else {
    dataVars <- colnames(firstImp)
  }
  if (all(trtVars %in% dataVars)==FALSE) {
    stop("Some of the treatment variables you specified are not in the data frame.")
  }

  if (missingData==TRUE) {
    n <- nrow(firstImp)
  } else {
    n <- nrow(data)
  }

  if (is.null(nSim)) {
    nSim <- n
  }

  #create blank dataset with treatment indicators set as per desired regime
  if (missingData==TRUE) {
    syntheticDataBlank <- data.frame(matrix(NA,nrow=nSim*numRegimes,
                                            ncol=ncol(firstImp)))
    colnames(syntheticDataBlank) <- colnames(firstImp)
  } else {
    syntheticDataBlank <- data.frame(matrix(NA,nrow=nSim*numRegimes,ncol=ncol(data)))
    colnames(syntheticDataBlank) <- colnames(data)
  }

  #create new variable which in the end will indicate which treatment regime
  #the row corresponds to
  syntheticDataBlank$regime <- as.factor(0)

  if (numRegimes==1) {
    syntheticDataBlank$regime <- as.factor(1)
    for (i in 1:timePoints) {
      #set treatment indicator according to specified regime
      syntheticDataBlank[1:n,trtVars[i]] <- trtRegimes[i]
    }
  } else {
    syntheticDataBlank$regime <- as.factor(rep(1:numRegimes,each=nSim))
    for (j in 1:numRegimes) {
      for (i in 1:timePoints) {
        syntheticDataBlank[((j-1)*nSim+1):(j*nSim), trtVars[i]] <- trtRegimes[[j]][i]
      }
    }
  }

  #set up predictor matrix for mice, exploiting monotone pattern
  predMat <- mice::make.predictorMatrix(syntheticDataBlank)
  predMat[,] <- 1*(lower.tri(predMat))

  if (missingData==FALSE) {
    data$regime <- as.factor(0)
    inputData <- rbind(data,syntheticDataBlank)

    if (is.null(method)) {
      method <- mice::make.method(data=inputData,defaultMethod = c("norm", "logreg", "polyreg","polr"))
    } else {
      #add on an empty imputation method for the new variable regime
      method <- c(method,regime="")
    }

    if (is.null(predictorMatrix)) {
      predictorMatrix <- predMat
    } else {
      #need to append user provided predictor matrix with extra row and column
      #corresponding to new variable regime
      predMat[1:nrow(predictorMatrix),1:ncol(predictorMatrix)] <- predictorMatrix
      predictorMatrix <- predMat
    }

    imps <- mice::mice(data=inputData,
               method=method,
               predictorMatrix = predictorMatrix,m=M,maxit=1,
               printFlag = micePrintFlag)

    if (silent==FALSE) {
      print("Variables imputed using:")
      print(imps$method)
      print("Predictor matrix is set to:")
      print(imps$predictorMatrix)
    }

    #remove original data from imputations
    #but first declare regime=NULL to avoid CRAN note
    regime <- NULL
    returnImps <- mice::filter(imps, regime!=0)

  } else {

    #create empty 'long' dataframe that will store all the imputations
    imputedDatasetsLong <- data.frame(matrix(NA,nrow=M*nSim*numRegimes,
                      ncol=ncol(firstImp)))
    imputedDatasetsLong$regime <- factor(NA, levels=c(as.character(1:numRegimes)))
    colnames(imputedDatasetsLong) <- colnames(syntheticDataBlank)
    imputedDatasetsLong$.imp <- rep(1:M,each=nSim*numRegimes)
    imputedDatasetsLong$.id <- rep(1:(nSim*numRegimes), times=M)

    for (i in 1:M) {
      #impute data in synthetic part using mice, with m=1
      inputData <- mice::complete(data, action=i)
      inputData$regime <- as.factor(0)
      inputData <- rbind(inputData,syntheticDataBlank)

      if (i==1) {
        if (is.null(method)) {
          method <- mice::make.method(data=inputData,defaultMethod = c("norm", "logreg", "polyreg","polr"))
        } else {
          #add on an empty imputation method for the new variable regime
          method <- c(method,regime="")
        }
        if (is.null(predictorMatrix)) {
          predictorMatrix <- predMat
        } else {
          #need to append user provided predictor matrix with extra row and column
          #corresponding to new variable regime
          predMat[1:nrow(predictorMatrix),1:ncol(predictorMatrix)] <- predictorMatrix
          predictorMatrix <- predMat
        }
      }

      imps <- mice::mice(data=inputData,
                   method=method,
                   predictorMatrix = predictorMatrix,m=1,maxit=1,
                   printFlag = micePrintFlag)

      if (i==1) {
        if (silent==FALSE) {
          print("Variables imputed using:")
          print(imps$method)
          print("Predictor matrix is set to:")
          print(imps$predictorMatrix)
        }
      }

      #prepare single imputation for copying to imputeDatasetsLong
      imputedDataset <- mice::complete(imps,action=1)
      #remove original data (regime=0)
      imputedDataset <- imputedDataset[imputedDataset$regime!=0,]
      imputedDataset$regime <- droplevels(imputedDataset$regime)
      #copy single imputation into long dataframe
      imputedDatasetsLong[((i-1)*(nSim*numRegimes)+1):(i*(nSim*numRegimes)),1:ncol(inputData)] <- imputedDataset
    }

    #put 'original' data at top
    imputedDatasetsLong <- rbind(cbind(syntheticDataBlank,.imp=0,.id=1:(nSim*numRegimes)),
                                 imputedDatasetsLong)
    #turn back into a mids object
    returnImps <- mice::as.mids(imputedDatasetsLong)
    #copy over predictor matrix  used into imps
    returnImps$predictorMatrix <- predMat
    returnImps$method <- imps$method
  }
  #return the imputations
  returnImps


}

