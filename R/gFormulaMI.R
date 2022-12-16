#timePoints - number of follow-up timepoints (including baseline)
#trtRegimes - list of treatment regimes of interest


gFormulaImpute <- function(data, M=50, trtVarStem, timePoints, trtRegimes,
                           micePrintFlag=FALSE, nSim=NULL) {

  if (class(data)=="mids") {
    missingData <- TRUE
    print("Input data is a mice created multiple imputation object.")
  } else if (class(data)=="data.frame")  {
    missingData <- FALSE
    print("Input data is a regular data frame.")
  } else {
    print("Input dataset should either be a data frame or a mids object created by mice.")
    stop()
  }

  if (typeof(trtRegimes)=="list") {
    numRegimes <- length(trtRegimes)
    #check trtRegimes elements each have length equal to timePoints
  } else {
    numRegimes <- 1
    #check trtRegimes elements each have length equal to timePoints+1
    if (length(trtRegimes)!=timePoints) {
      stop("Length of treatment regime does not match number of timepoints.")
    }
  }

  n <- nrow(data)
  if (is.null(nSim)) {
    nSim <- n
  }

  #create new variable which in the end will indicate which treatment regime
  #the row corresponds to
  data$regime <- as.factor(0)

  #create blank dataset with treatment indicators set as per desired regime
  syntheticDataBlank <- data.frame(matrix(NA,nrow=nSim*numRegimes,ncol=ncol(data)))
  colnames(syntheticDataBlank) <- colnames(data)
  #create vector of treatment variables
  trtVarNames <- paste0(trtVarStem,0:(timePoints-1))

  if (numRegimes==1) {
    syntheticDataBlank$regime <- as.factor(1)
    for (i in 1:timePoints) {
      #set treatment indicator according to specified regime
      syntheticDataBlank[1:n,trtVarNames[i]] <- trtRegimes[i]
    }
  } else {
    syntheticDataBlank$regime <- as.factor(rep(1:numRegimes,each=nSim))
    for (j in 1:numRegimes) {
      for (i in 1:timePoints) {
        syntheticDataBlank[((j-1)*nSim+1):(j*nSim), trtVarNames[i]] <- trtRegimes[[j]][i]
      }
    }
  }
  #set up predictor matrix for mice
  predMat <- mice::make.predictorMatrix(data)
  predMat[,] <- lower.tri(predMat)

  if (missingData==FALSE) {
    inputData <- rbind(data,syntheticDataBlank)
    print("Imputing potential outcomes using sequential models.
          Please check variables are ordered in time as desired!")

    imps <- mice::mice(data=inputData,
               defaultMethod = c("norm", "logreg", "polyreg","polr"),
               predictorMatrix = predMat,m=M,maxit=1,
               printFlag = micePrintFlag)
  } else {

    #now impute potential outcomes
    imputedDatasets <- vector(mode = "list", length = M)
    for (i in 1:M) {
      #impute data in synthetic part using mice, with m=1
      inputData <- rbind(complete(intermediateImps,action=i),
                         syntheticDataBlank)

      imps <- mice(data=inputData,
                   defaultMethod = c("norm", "logreg", "polyreg","polr"),
                   predictorMatrix = predMat,m=1,maxit=1,
                   printFlag = FALSE)

      imputedDatasets[[i]] <- complete(imps, action=1)
      #extract just the synthetic part
      imputedDatasets[[i]] <- imputedDatasets[[i]][(n+1):(3*n),]
    }
  }

  #remove original data from imputations
  imps <- mice::complete(data=imps,action="long",include=TRUE)
  imps <- imps[imps$regime!=0,]
  #turn back into a mids object
  imps <- suppressWarnings(mice::as.mids(imps))
}

gFormulaAnalyse <- function(fits) {

  M <- length(fits$analyses)
  p <- length(fits$analyses[[1]]$coefficients)
  ests <- array(0, dim=c(M,p))
  within_vars <- array(0, dim=c(M,p))

  for (i in 1:M) {
      ests[i,] <- coefficients(fits$analyses[[i]])
      within_vars[i,] <- diag(vcov(fits$analyses[[i]]))
  }

  overall_ests <- colMeans(ests)
  v <- colMeans(within_vars)
  b <- diag(var(ests))
  total_var <- (1+1/M)*b - v
  #degrees of freedom
  df <- (M-1)*(1-(M*v)/((M+1)*b))

  resTable <- array(0, dim=c(p,8))
  resTable[,1] <- overall_ests
  resTable[,2] <- v
  resTable[,3] <- b
  resTable[,4] <- total_var
  resTable[,5] <- df
  #95% confidence interval
  resTable[,6] <- overall_ests - qt(0.975,df=df)*sqrt(total_var)
  resTable[,7] <- overall_ests + qt(0.975,df=df)*sqrt(total_var)
  #two sided p-value
  resTable[,8] <- 2*pt(-abs(overall_ests/sqrt(total_var)),df=df)

  colnames(resTable) <- c("Estimate", "Within", "Between", "Total", "df",
                          "95% CI L", "95% CI U", "p")
  rownames(resTable) <- names(coefficients(fits$analyses[[1]]))
  resTable
}
