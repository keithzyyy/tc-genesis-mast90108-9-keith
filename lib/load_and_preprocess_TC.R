X_TC_FP = here("data", "curated","X_TC_genesis.rds")
Y_TC_FP = here("data", "curated","y_TC_genesis.rds")
X_TC_EXIST = file.exists(X_TC_FP)
Y_TC_EXIST = file.exists(Y_TC_FP)

if (X_TC_EXIST & Y_TC_EXIST){
  cat("Preprocessed (X,y) already exist. Importing now.\n")
  X = readRDS(X_TC_FP)
  y = readRDS(Y_TC_FP)
} else {
  cat("Preprocessing the TC genesis dataset....\n")
  data_path = here("data", "envDataset_12h_10x10_with_mask.csv")
  # import the data
  data = read.csv(data_path)
  # fix basin NA category first (NA missing to "NA")
  data.NA.fixed = fix_NA_value_and_refactor_basin(data)
  X = data.NA.fixed[, ! colnames(data.NA.fixed) %in% ( c("basin", "TC_genesis") )]
  y = data.NA.fixed$TC_genesis
  
  saveRDS(X, here("data", "curated","X_TC_genesis.rds"))
  saveRDS(y, here("data", "curated","y_TC_genesis.rds"))
}