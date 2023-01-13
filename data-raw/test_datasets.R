#simulated dataset with some missing data, already imputed by mice
simDataMis <- simDataFullyObs
simDataMis$l1[1:5] <- NA
set.seed(7247)
simDataMisImps <- mice::mice(simDataMis,m=50)
usethis::use_data(simDataMisImps, internal = TRUE, overwrite=TRUE)





