n <- 500

expit <- function(x) exp(x)/(1+exp(x))


l0 <- rnorm(n)
a0 <- 1*(runif(n)<expit(l0))
l1 <- l0+a0+rnorm(n)
a1 <- 1*(runif(n)<expit(l1+a0))
l2 <- l1+a1+rnorm(n)
a2 <- 1*(runif(n)<expit(l2+a1))
y <- l2+a2+rnorm(n)

obsData <- data.frame(l0=l0,a0=a0,l1=l1,a1=a1,l2=l2,a2=a2,y=y)

rm(l0,a0,l1,a1,l2,a2,y)

imps <- gFormulaImpute(obsData,M=50,trtVarStem="a",timePoints=2,trtRegime=c(0,0,0))

fit <- with(imps, lm(y[regime!=0]~l0[regime!=0]))
