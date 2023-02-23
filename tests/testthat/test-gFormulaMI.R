#tests when passing a complete dataset to syntheticPool

test_that("gFormulaImpute and syntheticPool run when input data frame has no missing data", {
  expect_error({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataFullyObs,M=10,trtVars=c("a0","a1","a2"), trtRegimes=c(0,0,0))
    fits <- with(impRes, lm(y~1))
    syntheticPool(fits)
  }, NA)
})

test_that("Check 95% confidence interval includes true value, one regime 0,0,0", {
  expect_equal({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVars=c("a0","a1","a2"), trtRegimes=c(0,0,0))
    fits <- with(impRes, lm(y~1))
    res <- syntheticPool(fits)
    #true mean under regime is 0
    1*((res[1,6]<0) & (res[1,7]>0))
  }, 1)
})

test_that("Check 95% confidence interval includes true value, one regime 1,1,1", {
  expect_equal({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVars=c("a0","a1","a2"),trtRegimes=c(1,1,1))
    fits <- with(impRes, lm(y~1))
    res <- syntheticPool(fits)
    #true mean under regime is 3
    1*((res[1,6]<3) & (res[1,7]>3))
  }, 1)
})

test_that("gFormulaImpute runs with two regimes", {
  expect_error({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataFullyObs,M=10,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)))
    fits <- with(impRes, lm(y~factor(regime)))
    syntheticPool(fits)
  }, NA)
})

test_that("Check 95% confidence interval includes true value, two regimes", {
  expect_equal({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)))
    fits <- with(impRes, lm(y~factor(regime)))
    res <- syntheticPool(fits)
    #true mean of regime 1 is 3, and true mean of regime 0 is 0
    1*((res[2,6]<3) & (res[2,7]>3))
  }, 1)
})

test_that("If passed a regular data frame with missing values, it should error", {
  expect_error({
    modData <- simDataFullyObs
    modData$l1[1:5] <- NA
    impRes <- gFormulaImpute(data=modData,M=50,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)))
  })
})

test_that("Catching mismatch in treatment variables", {
  expect_error({
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0),c(1,1)))
  })
})

test_that("Catching mismatch in treatment variables", {
  expect_error({
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVars=c("a0","a1"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)))
  })
})

test_that("gFormulaImpute runs when manually specifying method vector", {
  expect_error({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataFullyObs,M=10,trtVars=c("a0","a1","a2"),
                             trtRegimes=c(0,0,0),
                             method=c("norm","", "norm","","norm", "", "norm"))
  }, NA)
})

test_that("When passed custom method, returned object method matches", {
  expect_equal({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataFullyObs,M=10,trtVars=c("a0","a1","a2"),
                             trtRegimes=c(0,0,0),
                             method=c("pmm","", "norm","","norm", "", "norm"))
    as.vector(impRes$method)
    #note value to compare against below includes extra """ for the added variable regime
  }, c("pmm","", "norm","","norm", "", "norm", ""))
})

test_that("gFormulaImpute runs when manually specifying predictor matrix", {
  expect_error({
    set.seed(7626)
    temp <- simDataFullyObs
    temp$regime <- factor(0)
    myPredMat <- mice::make.predictorMatrix(temp)
    myPredMat[,] <- lower.tri(myPredMat)
    myPredMat["l2",c("l0","a0")] <- 0
    impRes <- gFormulaImpute(data=simDataFullyObs,M=10,trtVars=c("a0","a1","a2"),
                             trtRegimes=c(0,0,0),
                             predictorMatrix=myPredMat)
  }, NA)
})

#tests when passing a mids object of multiple imputations

test_that("Synthetic imputation runs when passed a mids object", {
  expect_error({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataMisImps,M=10,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)))
  },NA)
})

test_that("Check 95% confidence interval includes true value,
          two regimes, missing data imputations as input", {
  expect_equal({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataMisImps,M=50,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)))
    fits <- with(impRes, lm(y~factor(regime)))
    res <- syntheticPool(fits)
    #true mean of regime 1 is 3, and true mean of regime 0 is 0
    1*((res[2,6]<3) & (res[2,7]>3))
  }, 1)
})

test_that("When passed custom method, returned object method matches (passing mids object)", {
  expect_equal({
    set.seed(7626)
    impRes <- gFormulaImpute(data=simDataMisImps,M=50,trtVars=c("a0","a1","a2"),
                             trtRegimes=list(c(0,0,0),c(1,1,1)),
                             method=c("pmm","", "norm","","norm", "", "norm"))
    as.vector(impRes$method)
    #note value to compare against below includes extra """ for the added variable regime
  },c("pmm","", "norm","","norm", "", "norm", ""))
})

test_that("Check that syntheticPool catches non-positive variances", {
  expect_error({
    set.seed(766)
    temp <- simDataFullyObs[1:10,]
    impRes <- gFormulaImpute(data=simDataFullyObs,M=2,trtVars=c("a0","a1","a2"), trtRegimes=c(0,0,0))
    fits <- with(impRes, lm(y~1))
    syntheticPool(fits)
  })
})

#simulation studies to assess unbiasedness and coverage

test_that("Check 95% CI coverage", {
  skip_on_cran()
  expect_equal({
    expit <- function(x) exp(x)/(1+exp(x))
    set.seed(7626)
    nSim <- 1000
    n <- 500
    ciRegime1 <- array(0, dim=c(nSim,2))
    ciContrast <- array(0, dim=c(nSim,2))

    for (i in 1:nSim) {

      l0 <- rnorm(n)
      a0 <- 1*(runif(n)<expit(l0))
      l1 <- l0+a0+rnorm(n)
      a1 <- 1*(runif(n)<expit(l1+a0))
      l2 <- l1+a1+rnorm(n)
      a2 <- 1*(runif(n)<expit(l2+a1))
      y <- l2+a2+rnorm(n)

      simData <- data.frame(l0=l0,a0=a0,l1=l1,a1=a1,l2=l2,a2=a2,y=y)

      impRes <- gFormulaImpute(data=simData,M=50,trtVars=c("a0","a1","a2"),
                               trtRegimes=list(c(0,0,0),c(1,1,1)))
      fits <- with(impRes, lm(y~factor(regime)))
      res <- syntheticPool(fits)
      ciRegime1[i,] <- c(res[1,6], res[1,7])
      ciContrast[i,] <- c(res[2,6], res[2,7])
    }
    ciRegime1Included <- sum((ciRegime1[,1]<0) & (ciRegime1[,2]>0))
    ciContrastIncluded <- sum((ciContrast[,1]<3) & (ciContrast[,2]>3))
    prop.test(ciRegime1Included,nSim)
    ciCov <- prop.test(ciContrastIncluded,nSim)
    #check 95% CI for 95% CI coverage for contrast of regimes includes 95%
    1*((ciCov$conf.int[1]<0.95) & (ciCov$conf.int[2]>0.95))
  }, 1)
})
