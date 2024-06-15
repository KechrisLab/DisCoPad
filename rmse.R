#Root-Mean-Squared-Error calculations
#Inputs
#data: a vector with estimated values
#true: the "true" value of the estimand

#Output: RMSE

rmse = function(data, true){
  RMSE = sqrt(mean((data - true)^2))
  return(RMSE)
}