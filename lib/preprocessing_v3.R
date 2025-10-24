#' Fixes NA value as an actual category "NA" in
#' the basin variable, and turn the basin variable
#' as a factor.
#' 
#' @param dta = the raw data frame from the envData dataset.
#' 
#'
#' @returns corrected data. 
fix_NA_value_and_refactor_basin = function(dta){
  
  dta$basin[is.na(dta$basin)] = "NA" # better readability and efficiency
  
  # might be a good idea to keep basin as a factor for modelling purposes!
  dta$basin = as.factor(dta$basin)
  
  return(dta)
  
}

#' Preprocess the tropical cyclone dataset for model training.
#' One hot encodes (dummy encoding into k-1 dummy variables out of 
#' k levels) the basin variable.
#' 
#' Primary use case for fix.NA=FALSE is for stratified sampling, since
#' stratified variable is 'basin', which is the variable that contains
#' NA values.
#' 
#' @param dta = the raw data frame from the envData dataset.
#' @param basin.ref = reference category that one will use for 'basin'
#' @param fix.NA = whether one wants to convert NA to category "NA" or not.
#'
#' @returns a list of 2 elements: the preprocessed data matrix with
#' an 1s intercept column, and a vector of binary responses

preprocess_v3 = function(dta,
                         basin.ref='NA',
                         fix.NA=TRUE){
  ## encode starting time
  start.time = Sys.time()
  
  if (fix.NA==TRUE){
    # encode missing data from 'basin' as an actual category
    # might be a good idea to keep basin as a factor for modelling purposes!
    dta = fix_NA_value_and_refactor_basin(dta)
    
  }
  
  # define reference level for basin
  dta$basin = relevel(dta$basin, basin.ref)
  
  # encode categories in 'basin' using one hot encoding.
  # exclude the reference category from being encoded 
  
  dta.encoded = cbind(dta, model.matrix(~ 1 + basin,
                                        data = dta))
  
  ## get data frame of predictors only: we can remove the original 'basin' column
  dta.predictors = subset(dta.encoded, select = -c(basin, TC_genesis))
  
  ## get data frame of the response only
  dta.response = dta.encoded$TC_genesis
  
  ## format data in a analysis amendable way
  X = dta.predictors
  y = dta.response
  
  X = data.matrix(X)
  y = as.matrix(y)
  
  ## encode ending time and report execution time
  end.time = Sys.time()
  time.diff = end.time - start.time
  cat("Overall time taken (minutes): ", as.numeric(time.diff, units = "mins"))
  
  return(list(predictors = X,
              response = y))
}




