# Use single parameter choice for all days, fit to training, predict on test
# Including GEOS-Chem and elevation (standardized) as predictors

# Cluster code for personal library:
# mylib = "/home/bkiani2/rlib36/"
# .libPaths(mylib)

library(LatticeKrig)
library(spam64)
library(fields)
library(sessioninfo)
library(dplyr)

session_info()
sessionInfo()

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

# days = 182:212
# prefix_read = "./rawrds/split_4_"
# prefix_save = "./aod_holdouts/testdata/test_4_"
# NC.buffer = 5 
# a.wght = 12
# nu = 0.1
# nlevel = 5
# NC = 15

fit_test_lk <- function(days = 182:212, 
                        prefix_read = "./rawrds/split_4_",
                        prefix_save = "./aod_holdouts/testdata/lk_test_4_",
                        NC.buffer = 5, a.wght = 12, nu = 0.1, 
                        nlevel = 5, NC = 15) {
  
  # This function will save predictions on the test data set
  # -- save betas, likelihood estimates as well, for diagnostics
  # Save as you go rather than all at once at the end
  
  coef.estimates = vector("list", length(days))
  names(coef.estimates) <- paste0("day", days)
  loglik.estimates <- coef.estimates
  summary_names <- c("day", "train_size", "test_size", 
                     "R2", "RMSE", "MSE", "Int", "Slope")
  summary_info <- matrix(NA, nrow = length(days), ncol = length(summary_names))
  colnames(summary_info) <- summary_names
  summary_info <- as.data.frame(summary_info)

  # Use this for LKrigSetup
  xy_full <- 
    as.matrix(
      readRDS("./rawrds/n4_day182.rds") %>%
      dplyr::select(cmaq_x, cmaq_y) %>%
      mutate(x_km = cmaq_x/1000, y_km = cmaq_y/1000) %>%
      dplyr::select(-cmaq_x, -cmaq_y)
    )
  
  # Use the full training dataset for setting boundaries
  # -- this lk_info does not change based on fold or day
  lk_info <- LKrigSetup(x = xy_full, 
                        nlevel = nlevel,
                        NC = NC,
                        NC.buffer = NC.buffer, 
                        nu = nu,
                        a.wght = a.wght)
  
  rm(xy_full)
  
  for (day in days) {

    print(day)

    # Read in full data
    x <- readRDS(paste0(prefix_read, day, ".rds"))
    testthat::expect_identical(x$day[1], day)
    
    # Covariates used for mean model
    # spatial drift (cmaq_x, cmaq_y) included by default
    # More: gc_aod + cmaq_x*cmaq_y 
    #       + cmaq_x * gc_aod + cmaq_y * gc_aod
    #       + elevation (normalized)
    
    x_simple <- x %>%
      dplyr::select(cmaq_id, split, 
                    aod_value, cmaq_x, cmaq_y, gc_aod, elev) %>%
      # For numerical stability, putting coordinates in terms of KM
      mutate(x_km = cmaq_x/1000, y_km = cmaq_y/1000) %>%
      dplyr::select(-cmaq_x, -cmaq_y) %>%
      mutate(
        xy_int = x_km * y_km,
        x_gc_aod = x_km*gc_aod,
        y_gc_aod = y_km*gc_aod,
        elev_norm = (elev - mean(elev))/sd(elev)
      ) %>%
      dplyr::select(-elev)
    
    # Covariate names should exclude x/y -- included by default in LK fits
    cov_names <- 
      names(x_simple)[!(names(x_simple) %in% 
                          c("x_km", "y_km", "aod_value", 
                            "split", "cmaq_id"))]
    
    # Create y/xy/Z matrices for training and test data
    df_train <- x_simple %>% filter(split != "Test")
    df_test <- x_simple %>% filter(split == "Test")
    y_train <- as.vector(df_train$aod_value)
    y_test <- as.vector(df_test$aod_value)
    xy_train <- as.matrix(df_train[, c("x_km", "y_km")])
    xy_test <- as.matrix(df_test[, c("x_km", "y_km")])
    Z_train <- as.matrix(df_train[, cov_names])
    Z_test <- as.matrix(df_test[, cov_names])
    
    
    df_save <- df_test %>%
      dplyr::select(cmaq_id, x_km, y_km, aod_value)
    df_save$lk_predict <- NA

    # Fit LatticeKrig
    lk_fit <- LatticeKrig(x = xy_train, y = y_train, 
                          Z = Z_train, LKinfo = lk_info)
    
    # Predict
    lk_pred <-  as.vector(predict.LKrig(lk_fit, xnew = xy_test, Znew = Z_test))
    
    # Save the predictions correctly
    df_save$lk_predict <- lk_pred
      
    # Save test predictions
    saveRDS(df_save, 
            paste0(prefix_save, day, ".rds"))
    
    day_name <- paste0("day", day)
    coef.estimates[[day_name]] <- lk_fit$d.coef
    loglik.estimates[[day_name]] <- lk_fit$MLE$summary
    
    day_index <- which(days == day)
    summary_info$day[day_index] <- day
    summary_info$train_size[day_index] <- nrow(df_train)
    summary_info$test_size[day_index] <- nrow(df_test)
    summary_info$R2[day_index] <- getR2(true = y_test, pred = lk_pred)
    summary_info$RMSE[day_index] <- getRMSE(true = y_test, pred = lk_pred)
    summary_info$MSE[day_index] <- getMSE(true = y_test, pred = lk_pred)
    summary_info[day_index, c("Int", "Slope")] <- 
      getIntSlope(true = y_test, pred = lk_pred)
    
    
    rm(df_save, y_train, y_test, xy_train, xy_test, Z_train, Z_test,
       lk_fit, lk_pred, x_simple, x, df_test, df_train)
    
  }
  
  rm(lk_info)
  lk_results <- list(summary_info = summary_info,
                     coef.estimates = coef.estimates,
                     loglik.estimates = loglik.estimates)
  return(lk_results)
  
}



#### Fit ####

lk_results <- 
  fit_test_lk(days = 182:212, 
              prefix_read = "./rawrds/split_4_",
              prefix_save = "./aod_holdouts/testdata/lk_test_4_",
              NC.buffer = 5, a.wght = 12, nu = 0.1, 
              nlevel = 5, NC = 15)

saveRDS(lk_results, "./aod_holdouts/testdata/lk_test_4_results.rds")
gc()