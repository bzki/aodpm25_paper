library(dplyr)
library(ranger)
library(sessioninfo)
sessionInfo()
session_info()

# Functions to calculate R2, MSE
getR2 <- function(true, pred) {
  return(summary( lm( true ~ pred ) )$r.squared)
}
getMSE <- function(true, pred) {
  return(mean( (true - pred)^2 ))
}
getRMSE <- function(true, pred) {
  return(sqrt(mean( (true - pred)^2 )))
}
getIntSlope <- function(true, pred) {
  return(as.vector(coef(lm(true ~ pred))))
}

fit_test_rf <- function(days = 182:212, 
                        prefix_read = "./rawrds/split_4_",
                        prefix_save = "./aod_holdouts/testdata/rf_test_4_",
                        mtry = 7, min.node.size = 5, num.trees = 2000) {
  
  
  summary_names <- c("day", "train_size", "test_size", 
                     "R2", "RMSE", "MSE", "Int", "Slope")
  summary_info <- matrix(NA, nrow = length(days), ncol = length(summary_names))
  colnames(summary_info) <- summary_names
  summary_info <- as.data.frame(summary_info)

  for (day in days) {
  
    print(day)  
    # Read in full data
    x <- readRDS(paste0(prefix_read, day, ".rds"))
    testthat::expect_identical(x$day[1], day)
    
    df_save <- x %>%
      filter(split == "Test") %>%
      dplyr::select(cmaq_id, cmaq_x, cmaq_y, aod_value)
    df_save$rf_predict <- NA
    
    x2 <- x %>%
      dplyr::select(-year, -day, -lon, -lat, -cmaq_id, 
                    -starts_with("narr"), -pm25_value)
    
    df_train <- x2 %>%
      filter(split != "Test") %>%
      dplyr::select(-split)
    
    df_test <- x2 %>%
      filter(split == "Test") %>%
      dplyr::select(-split)
    
    y_test <- df_test$aod_value
    
    testthat::expect_equal(df_test$cmaq_x, df_save$cmaq_x)
    testthat::expect_equal(df_test$cmaq_y, df_save$cmaq_y)
    
    set.seed(day)
    ranger_fit <- ranger(aod_value ~ ., df_train, 
                         write.forest = TRUE, # Need for prediction
                         importance = "none", # Not recording importance
                         num.trees = num.trees, 
                         mtry = mtry,
                         min.node.size = min.node.size)
      
    ranger_predict <- predict(ranger_fit, data=df_test)$predictions
    
    # Save the predictions correctly
    df_save$rf_predict <- ranger_predict
    
    saveRDS(df_save, 
            paste0(prefix_save, day, ".rds"))
    
    day_index <- which(days == day)
    summary_info$day[day_index] <- day
    summary_info$train_size[day_index] <- nrow(df_train)
    summary_info$test_size[day_index] <- nrow(df_test)
    summary_info$R2[day_index] <- getR2(true = y_test, pred = ranger_predict)
    summary_info$RMSE[day_index] <- getRMSE(true = y_test, pred = ranger_predict)
    summary_info$MSE[day_index] <- getMSE(true = y_test, pred = ranger_predict)
    summary_info[day_index, c("Int", "Slope")] <- 
      getIntSlope(true = y_test, pred = ranger_predict)
      
    rm(x, x2, df_train, df_test, ranger_fit, ranger_predict, df_save, y_test)
      
    }

  return(summary_info)
}

# debugonce(fit_test_rf)
# test_results <-
#   fit_test_rf(days = 182:183, 
#               prefix_read = "./rawrds/split_4_",
#               prefix_save = "./aod_holdouts/testdata/rf_test_4_",
#               mtry = 7, min.node.size = 5, num.trees = 2000)

rf_results <-
  fit_test_rf(days = 182:212, 
              prefix_read = "./rawrds/split_4_",
              prefix_save = "./aod_holdouts/testdata/rf_test_4_",
              mtry = 7, min.node.size = 5, num.trees = 2000)

saveRDS(rf_results, "./aod_holdouts/testdata/rf_test_4_results.rds")
gc()
