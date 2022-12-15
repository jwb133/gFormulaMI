#timePoints - number of follow-up timepoints (not including baseline)
#trtRegime - list of regimes of interest, but in first interest one regime


gFormulaImpute <- function(data, M=50, trtVarStem, timePoints, trtRegime,
                           micePrintFlag=FALSE) {

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

  #check trtRegime elements each have length equal to timePoints

  n <- nrow(data)

  #create new variable which in the end will indicate which treatment regime
  #the row corresponds to
  data$regime <- as.factor(0)

  #create blank dataset with treatment indicators set as per desired regime
  syntheticDataBlank <- data
  #set everything to missing
  syntheticDataBlank[,] <- NA
  #create vector of treatment variables
  trtVarNames <- paste0(trtVarStem,0:timePoints)
  syntheticDataBlank$regime <- as.factor(1)
  for (i in 1:(timePoints+1)) {
    #set treatment indicator according to specified regime
    syntheticDataBlank[1:n,trtVarNames[i]] <- 0
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
  imps <- mice::as.mids(imps)
}
#
# gFormulaAnalyse <- function(imps,analysis) {
#
#   M <- max(imps$.imp)
#   fit <- do.call(what=analysis, argswith(imps[imps$.imp==1,], expr=analysis)
#   p <- length(coef(fit))
#   ests <- array(0, dim=c(M,p))
#   within_vars <- array(0, dim=c(M,p))
#
#   for (i in 1:M) {
#       fit <- with(imps[imps$.imp==i,], expr=analysis)
#       ests[i,] <- coef(fit)
#       within_vars[i,] <- diag(vcov(fit))
#   }
#
#   v <- colMeans(within_vars)
#   b <- var(ests)
#   total_var <- (1+1/M)*b - v
#   print(total_var)
# }
