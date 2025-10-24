#' Preprocess the tropical cyclone dataset for model training.
#' 
#' @param dta = the raw data frame from the envData dataset.
#'
#' @returns a list of 2 elements: the preprocessed data matrix with
#' an 1s intercept column, and a vector of binary responses

preprocess = function(dta){
  ## encode starting time
  start.time = Sys.time()
  
  ## encode missing data from 'basin' as an actual category
  #dta$basin = sapply(dta$basin, function(cat){if(is.na(cat)){"NA"}else{cat}})
  dta$basin[is.na(dta$basin)] = "NA" # better readability and efficiency
  
  # might be a good idea to keep basin as a factor for modelling purposes!
  dta$basin = as.factor(dta$basin)
  
  # encode categories in 'basin' using one hot encoding
  dta.encoded = cbind(dta, model.matrix(~ 0 + basin, data = dta))
  
  ## get data frame of predictors only: we can remove the original 'basin' column
  dta.predictors = subset(dta.encoded, select = -c(basin, TC_genesis))
  
  ## get data frame of the response only
  dta.response = dta.encoded$TC_genesis
  
  ## format data in a analysis amendable way
  X = dta.predictors
  y = dta.response
  
  X = data.matrix(X)
  # add intercept
  X = cbind(matrix(1, nrow=nrow(X), ncol=1), X)
  y = as.matrix(y)
  
  ## encode ending time and report execution time
  end.time = Sys.time()
  time.diff = end.time - start.time
  cat("Overall time taken (minutes): ", as.numeric(time.diff, units = "mins"))
  
  return(list(predictors = X,
              response = y))
}