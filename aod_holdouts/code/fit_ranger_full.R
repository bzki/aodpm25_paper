# Fit random forest to full training data
# Then predict on the test datasets and record MSE and other metrics
# Picking mtry = 7 and nodesize = 5

library(dplyr)
library(ranger)
# library(fields)
library(sessioninfo)
sessionInfo()
session_info()


fit_ranger <- function(days = 182:212, NUMTREES = 1000,
                       mtry, nodesize) {
  
  # Placeholder for predictions, results
  rf_matrix_names <- c("day", "mtry", "nodesize", "NUMTREES", 
                       # Metrics for understanding performance
                       "R2_oob", "MSE_oob", 
                       # Quantile of predictions
                       "p_min", "p_p01", "p_p02", "p_p98", "p_p99", "p_max",
                       # Quantile of full predictions
                       "f_min", "f_p01", "f_p02", "f_p98", "f_p99", "f_max",
                       # Fitting times from ranger
                       "fit_system", "fit_user", "fit_elapsed",
                       # Additional information on holdout/validation
                       "train_size", "holdout_size")
  rf_matrix <- matrix(NA, nrow = length(days), ncol = length(rf_matrix_names))
  colnames(rf_matrix) <- rf_matrix_names
  rf_matrix <- as.data.frame(rf_matrix)
  rf_matrix[, "day"] <- days
  
  # Other information to save
  # - Full and partial predictions 
  # For full predictions, ideally we have some way of identifying test data
  ppredict_list = vector("list", length = length(days))
  names(ppredict_list) = paste0("day", days)
  fpredict_list = ppredict_list

  for (day in days) {
    
    print(day)
    day_name = paste0("day", day)
    day_index = which(days == day)
    
    x <- readRDS(paste0("./rawrds/n4_day", day, ".rds"))
    testthat::expect_identical(x$day[1], day)

    full_id <- x$cmaq_id
    holdout_index <- which(is.na(x$aod_value))
    train_id <- full_id[-holdout_index]
    test_id <- full_id[holdout_index] 

    x2 <- x %>%
      dplyr::select(-year, -day, -lon, -lat, -cmaq_id, 
                    -starts_with("narr"), -pm25_value)
    
    df_train <- x2[-holdout_index, ]
    df_test <- x2[holdout_index, ]


    xy_train <- as.matrix(df_train[, c("cmaq_x", "cmaq_y")])
    xy_test <- as.matrix(df_test[, c("cmaq_x", "cmaq_y")])
    
    y_train <- as.vector(df_train$aod_value)
    y_test <- as.vector(df_test$aod_value)

    df_full <- x2
    xy_full <- as.matrix(x2[, c("cmaq_x", "cmaq_y")])
    y_full <- as.vector(df_full$aod_value)

    # Save training/testing size
    rf_matrix[day_index, "holdout_size"] <- length(y_test)
    rf_matrix[day_index, "train_size"] <- length(y_train)
    
    # Create prediction matrix for day
    ppredict_colnames <- c("cmaq_x", "cmaq_y", "aod_value", "ranger_predict")
    ppredict_mat <- matrix(NA, nrow = nrow(xy_test), 
                           ncol = length(ppredict_colnames))
    colnames(ppredict_mat) <- ppredict_colnames
    ppredict_mat <- as.data.frame(ppredict_mat)
    ppredict_mat[, c("cmaq_x", "cmaq_y")] <- xy_test
    ppredict_mat[, c("aod_value")] <- y_test
    
    fpredict_colnames <- c("cmaq_id", "cmaq_x", "cmaq_y", 
                           "aod_value", "holdout", "ranger_predict")
    fpredict_mat <- matrix(NA, nrow = nrow(xy_full), 
                           ncol = length(fpredict_colnames))
    colnames(fpredict_mat) <- fpredict_colnames
    fpredict_mat <- as.data.frame(fpredict_mat)

    fpredict_mat[, c("cmaq_x", "cmaq_y")] <- xy_full
    fpredict_mat[, "aod_value"] <- y_full
    fpredict_mat[, "cmaq_id"] <- full_id
    fpredict_mat[, "holdout"] <- NA # NA for not observed
    fpredict_mat[holdout_index, "holdout"] <- 1L # 1 for testing data
    fpredict_mat[-holdout_index, "holdout"] <- 0L # 0 for training data
    
    rf_matrix[day_index, "mtry"] <- mtry
    rf_matrix[day_index, "nodesize"] <- nodesize
    rf_matrix[day_index, "NUMTREES"] <- NUMTREES
    
    # Use seed for reproducibility
    set.seed(11223344)

    ranger_time <- 
    system.time(
      ranger_fit <- ranger(aod_value ~ ., df_train, 
                           write.forest = TRUE, # Need for prediction
                           importance = "none", # Not recording importance
                           num.trees = NUMTREES, 
                           mtry = mtry,
                           min.node.size = nodesize)
    )
    # Save time
    rf_matrix[day_index, c("fit_system", "fit_user", "fit_elapsed")] <- 
      as.vector(ranger_time[c(1,2,3)])
    # Get predictions and record R2 and MSE
    ranger_predict <- predict(ranger_fit, data=df_test)$predictions
    rf_matrix[day_index, "MSE_oob"] <- ranger_fit$prediction.error
    rf_matrix[day_index, "R2_oob"] <- ranger_fit$r.squared

    ranger_fpredict <- predict(ranger_fit, data=df_full)$predictions

    # Record quantiles of predictions on unobserved data
    rf_matrix[day_index, c("p_min", "p_p01", "p_p02", "p_p98", "p_p99", "p_max")] <- 
      c(min(ranger_predict),
        quantile(ranger_predict, probs = c(0.01, 0.02, 0.98, 0.99)),
        max(ranger_predict))

    rf_matrix[day_index, c("f_min", "f_p01", "f_p02", "f_p98", "f_p99", "f_max")] <- 
      c(min(ranger_fpredict),
        quantile(ranger_fpredict, probs = c(0.01, 0.02, 0.98, 0.99)),
        max(ranger_fpredict))

    # Save predictions for mapping
    ppredict_mat[, "ranger_predict"] <- ranger_predict
    
    # Generate a full prediction and save
    fpredict_mat[, "ranger_predict"] <- ranger_fpredict
        
    # Save prediction matrices
    ppredict_list[[day_name]] <- ppredict_mat
    fpredict_list[[day_name]] <- fpredict_mat
    
    rm(ppredict_mat, fpredict_mat, ranger_fit, ranger_predict, ranger_fpredict,
       x, x2, xy_full, xy_test, xy_train, df_full, df_train, df_test,
       y_train, y_test, y_full,
       day_name, day_index, test_id, train_id, 
       holdout_index, full_id)

  }
  
  ranger_return <- list(ranger_matrix = rf_matrix,
                        ppredict_list = ppredict_list,
                        fpredict_list = fpredict_list)
  return(ranger_return)
  
}



# Test run
# debugonce(fit_ranger)
# test_1 <- fit_ranger(days = 182:183, NUMTREES = 50, mtry = 7, nodesize = 5)

rf_results <- fit_ranger(days = 182:212, 
                         NUMTREES = 2000, mtry = 7, nodesize = 5)
saveRDS(rf_results, "./aod_holdouts/finalpred/rf_final.rds")
rm(rf_results)
gc()