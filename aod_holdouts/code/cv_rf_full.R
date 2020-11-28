# Run on full observed data cross-validation folds
# Predict 10 spatially clustered cross-validation folds with random forest
# (using ranger's implementation)

# Cluster code for personal library
# mylib = "/home/bkiani2/rlib36/"
# .libPaths(mylib)

library(dplyr)
library(ranger)
library(sessioninfo)
sessionInfo()
session_info()

fit_cv_rf_full <- function(days = 182:212, 
                           prefix_read = "./rawrds/split_fcv_",
                           prefix_save = "./aod_holdouts/cvdata_full/rf_fcv_",
                           mtry = 7, min.node.size = 5, num.trees = 1000) {
  
  
  for (day in days) {
  
    print(day)  
    # Read in training data
    train_list <- readRDS(paste0(prefix_read, day, ".rds"))
    x <- train_list$df
    folds <- train_list$folds
    rm(train_list)    
    
    x_basic <- x %>%
      dplyr::select(cmaq_id, foldID, cmaq_x, cmaq_y, cv_nndist, aod_value)
    x_basic$rf_cv <- NA
    
    x2 <- x %>%
      dplyr::select(-year, -day, -lon, -lat, -cmaq_id, 
                    -starts_with("narr"), -pm25_value,
                    -foldID, -cv_nndist)
    
    rm(x)
    
    num_folds <- length(folds)
    
    for (fi in seq(num_folds)) {
      
      print(paste0("Fold ", fi))
      
      x_train <- x2[folds[[fi]][[1]], ]
      x_cv <- x2[folds[[fi]][[2]], ]
      
      # Ensure that the fold is correctly encoded
      testthat::expect_true(
        dplyr::setequal(which(x_basic$foldID == fi), 
                        folds[[fi]][[2]])
      )
      
      set.seed(day * fi)
      ranger_fit <- ranger(aod_value ~ ., x_train, 
                           write.forest = TRUE, # Need for prediction
                           importance = "none", # Not recording importance
                           num.trees = num.trees, 
                           mtry = mtry,
                           min.node.size = min.node.size)
      
      ranger_predict <- predict(ranger_fit, data=x_cv)$predictions
      
      # Save the predictions correctly
      x_basic$rf_cv[folds[[fi]][[2]]] <- ranger_predict
      
      testthat::expect_equal(x_basic$cmaq_x[folds[[fi]][[2]]], 
                             x_cv$cmaq_x)
      testthat::expect_equal(x_basic$cmaq_y[folds[[fi]][[2]]], 
                             x_cv$cmaq_y)
      testthat::expect_equal(x_basic$aod_value[folds[[fi]][[2]]], 
                             x_cv$aod_value)
      
      rm(x_train, x_cv, ranger_fit, ranger_predict)
      
    }
    
    testthat::expect_equal(0, sum(is.na(x_basic$rf_cv)))
    
    # Save x_basic
    saveRDS(x_basic, 
            paste0(prefix_save, day, ".rds"))
    
    rm(x_basic, x2, num_folds, folds, fi)
    
  }
  
  return(list("Success"))
}

# debugonce(fit_cv_rf_full)
# fit_cv_rf_full(days = 182L,
#                prefix_read = "./rawrds/split_fcv_",
#                prefix_save = "./aod_holdouts/cvdata_full/rf_fcv_",
#                mtry = 7, min.node.size = 5, num.trees = 50)


fit_cv_rf_full(days = 182:212, 
               prefix_read = "./rawrds/split_fcv_",
               prefix_save = "./aod_holdouts/cvdata_full/rf_fcv_",
               mtry = 7, min.node.size = 5, num.trees = 1000)

gc()