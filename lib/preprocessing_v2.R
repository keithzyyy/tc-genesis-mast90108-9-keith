

#' Preprocess the tropical cyclone dataset for model training with
#' FULL rank encoding and NO intercept. 
#'
#' 
#' @param dta = the raw data frame from the envData dataset.
#' @param fix.NA = if NAs already rectified as a proper factor variable. 
#' Primary use case for FALSE is for stratify sampling.
#'
#' @returns a list of 2 elements: the preprocessed data matrix with
#' an 1s intercept column, and a vector of binary responses

preprocess_v2 = function(dta,
                         fix.NA=TRUE){
  ## encode starting time
  start.time = Sys.time()
  
  if (fix.NA == TRUE){
    ## encode missing data from 'basin' as an actual category
    #dta$basin = sapply(dta$basin, function(cat){if(is.na(cat)){"NA"}else{cat}})
    dta$basin[is.na(dta$basin)] = "NA" # better readability and efficiency
    
    # might be a good idea to keep basin as a factor for modelling purposes!
    dta$basin = as.factor(dta$basin)
  }

  
  # encode categories in 'basin' using full rank encoding with NO intercept
  dta.encoded = cbind(dta, model.matrix(~ 0 + basin, data = dta))
  
  ## get data frame of predictors only: we can remove the original 'basin' column
  dta.predictors = subset(dta.encoded, select = -c(basin, TC_genesis))
  
  ## get data frame of the response only
  dta.response = dta.encoded$TC_genesis
  
  ## format data in a analysis amendable way
  X = dta.predictors
  y = dta.response
  
  X = data.matrix(X)
  # nuh-uh! full rank encoding means no intercept!
  #X = cbind(matrix(1, nrow=nrow(X), ncol=1), X)
  y = as.matrix(y)
  
  ## encode ending time and report execution time
  end.time = Sys.time()
  time.diff = end.time - start.time
  cat("Overall time taken (minutes): ", as.numeric(time.diff, units = "mins"))
  
  return(list(predictors = X,
              response = y))
}