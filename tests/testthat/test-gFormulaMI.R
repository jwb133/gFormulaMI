test_that("gFormulaImpute and gFormulaAnalyse runs when input data frame has no missing data", {
  expect_error({
    impRes <- gFormulaImpute(data=simDataFullyObs,M=10,trtVarStem="a", timePoints=2, trtRegime=c(0,0))
    fits <- with(impRes, lm(y~1))
    gFormulaAnalyse(fits)
  }, NA)
})

test_that("Check 95% confidence interval includes true value, one regime", {
  expect_equal({
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVarStem="a", timePoints=2, trtRegime=c(0,0))
    fits <- with(impRes, lm(y~1))
    res <- gFormulaAnalyse(fits)
    #true mean under regime is 0
    1*((res[1,6]<0) & (res[1,7]>0))
  }, 1)
})

test_that("gFormulaImpute runs with two regimes", {
  expect_error({
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVarStem="a",timePoints=2,
                             trtRegime=list(c(0,0),c(1,1)))
    fits <- with(impRes, lm(y~factor(regime)))
    gFormulaAnalyse(fits)
  }, NA)
})

test_that("Check 95% confidence interval includes true value, two regimes", {
  expect_equal({
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVarStem="a",timePoints=2,
                             trtRegime=list(c(0,0),c(1,1)))
    fits <- with(impRes, lm(y~factor(regime)))
    res <- gFormulaAnalyse(fits)
    #true mean of regime 1 is 3, and true mean of regime 0 is 0
    1*((res[2,6]<3) & (res[2,7]>3))
  }, 1)
})

test_that("If passed a regular data frame with missing values, it should error", {
  expect_error({
    simDataFullyObs$l1[1:5] <- NA
    impRes <- gFormulaImpute(data=simDataFullyObs,M=50,trtVarStem="a",timePoints=2,
                             trtRegime=list(c(0,0),c(1,1)))
  })
})
