# Choose specific parameters for LatticeKrig
# Fit on each day and predict full-map
# nu = 0.1, a.wght = 12.0, nlevel = 5, NC = 15

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

fit_lk <- function(days = 182:212,
                   NC.buffer = 5, a.wght, nu, nlevel, NC) {
  
  # Placeholder for predictions, results
  lk_matrix_names <- c("day", "nlevel", "NC", "nu", "a.wght",
                       "lnproflike", 
                       # Quantile of predictions on full predictions
                       "f_min", "f_p01", "f_p02", "f_p98", "f_p99", "f_max",
                       # Quantile of predictions on non-observed data,
                       "p_min", "p_p01", "p_p02", "p_p98", "p_p99", "p_max",
                       # Training observed max/min
                       "train_min", "train_max", 
                       # Beta estimates
                       "beta_0", "beta_x", "beta_y", "beta_gc", "beta_xy",
                       "beta_xgc", "beta_ygc", "beta_elev",
                       # Fitting times from LatticeKrig
                       "fit_system", "fit_user", "fit_elapsed",
                       # Additional information on holdout/validation/basis count
                       "train_size", "holdout_size", "basis_count")

  lk_matrix = matrix(NA, nrow = length(days), 
                      ncol = length(lk_matrix_names))
  # class(lk_matrix) = "numeric"
  colnames(lk_matrix) <- lk_matrix_names
  lk_matrix <- as.data.frame(lk_matrix)
  lk_matrix$day <- days

  ppredict_list = vector("list", length = length(days))
  names(ppredict_list) = paste0("day", days)
  fpredict_list = ppredict_list
  
  # Use the full map for the basis function placement
  # -- does not change day to day
  xy_all <- as.matrix(readRDS(paste0("./rawrds/n4_day", 182, ".rds")) %>%
    dplyr::select(cmaq_x, cmaq_y) %>%
    mutate(x_km = cmaq_x/1000, y_km = cmaq_y/1000) %>%
    dplyr::select(-cmaq_x, -cmaq_y))
  
  lk_info <- LKrigSetup(x = xy_all, 
                        nlevel = nlevel,
                        NC = NC,
                        NC.buffer = NC.buffer, 
                        nu = nu,
                        a.wght = a.wght)
  
  rm(xy_all)
  
  grid_info <- c(lk_info$latticeInfo$delta)
  
  for (day in days) {

    day_name = paste0("day", day)
    day_index <- which(lk_matrix$day == day)
    print(day)
    
    # Read in full data
    x <- readRDS(paste0("./rawrds/n4_day", day, ".rds"))
    full_id <- x$cmaq_id
    testthat::expect_identical(x$day[1], day)
    holdout_index <- which(is.na(x$aod_value))
    train_id <- full_id[-holdout_index]
    test_id <- full_id[holdout_index] 

    # Covariates used for mean model
    # spatial drift (cmaq_x, cmaq_y) included by default
    # More: gc_aod + cmaq_x*cmaq_y 
    #       + cmaq_x * gc_aod + cmaq_y * gc_aod
    #       + elevation (standardized)
    x <- x %>%
      dplyr::select(aod_value, cmaq_x, cmaq_y, gc_aod, elev) %>%
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
    # Exclude aod_value, cmaq_x, cmaq_y (coords included already)
    cov_names <- names(x)[!(names(x) %in% c("x_km", "y_km", "aod_value"))]

    # Create training set using holes with particular radius
    df_train <- x[-holdout_index, ]
    df_test <- x[holdout_index, ]
    
    xy_train <- as.matrix(df_train %>% dplyr::select(x_km, y_km))
    Z_train <- as.matrix(df_train[, cov_names])
    y_train <- as.vector(df_train$aod_value)
    
    xy_test <- as.matrix(df_test %>% dplyr::select(x_km, y_km))
    Z_test <- as.matrix(df_test[, cov_names])
    y_test <- as.vector(df_test$aod_value)

    xy_full <- as.matrix(x %>% dplyr::select(x_km, y_km))
    Z_full <- as.matrix(x[, cov_names])
    y_full <- as.vector(x$aod_value)
    
    lk_matrix[day_index , "train_min"] <- min(y_train)
    lk_matrix[day_index , "train_max"] <- max(y_train)
    lk_matrix[day_index, "holdout_size"] <- length(y_test)
    lk_matrix[day_index, "train_size"] <- length(y_train)
    
    # Create prediction matrix for day
    ppredict_colnames <- c("x_km", "y_km", "aod_value", "lk_predict")
    ppredict_mat <- matrix(NA, nrow = nrow(xy_test), 
                           ncol = length(ppredict_colnames))
    colnames(ppredict_mat) <- ppredict_colnames
    ppredict_mat <- as.data.frame(ppredict_mat)
    # class(ppredict_mat) <- "numeric"
    ppredict_mat[, c("x_km", "y_km")] <- xy_test
    ppredict_mat[, c("aod_value")] <- y_test
    ppredict_mat$day <- day
    
    fpredict_colnames <- c("cmaq_id", "x_km", "y_km", 
                           "aod_value", "holdout", "lk_predict")
    fpredict_mat <- matrix(NA, nrow = nrow(xy_full), 
                           ncol = length(fpredict_colnames))
    colnames(fpredict_mat) <- fpredict_colnames
    fpredict_mat <- as.data.frame(fpredict_mat)
    # class(fpredict_mat) <- "numeric"
    # Changing from matrix to data.frame to allow for different data types
    # -- should not change any of the code otherwise
    fpredict_mat[, c("x_km", "y_km")] <- xy_full
    fpredict_mat[, "aod_value"] <- y_full
    fpredict_mat[, "cmaq_id"] <- full_id
    fpredict_mat[, "holdout"] <- NA 
    fpredict_mat[holdout_index, "holdout"] <- 1L # 1 for unobserved data
    fpredict_mat[-holdout_index, "holdout"] <- 0L # 0 for training data
    fpredict_mat$day <- day
    
    lk_matrix[day_index, "nu"] <- nu
    lk_matrix[day_index, "a.wght"] <- a.wght
    lk_matrix[day_index, "NC"] <- NC
    lk_matrix[day_index, "nlevel"] <- nlevel
    
    lk_matrix[day_index, "basis_count"] <- lk_info$latticeInfo$m
    
    
    # There may be singularity or numerical issues, hence tryCatch
    lk_time <- system.time(
      lk_fit <- tryCatch(
        LatticeKrig(x = xy_train, y = y_train, 
                    Z = Z_train, LKinfo = lk_info),
        error = function(e) {
          print(e)
          return(NULL)
        }
      )
    )
    
    if (!is.null(lk_fit)) {
      
      # Get predictions
      lk_pred_time <- system.time(
        lk_pred <-  as.vector(
          predict.LKrig(lk_fit,
                        xnew = xy_test,
                        Znew = Z_test))
      )
      
      lk_fullpred <- 
          as.vector(predict.LKrig(lk_fit, xnew = xy_full, Znew = Z_full))
      
      # Save time
      lk_matrix[day_index, c("fit_system", "fit_user", "fit_elapsed")] <- 
        as.vector(lk_time[c(1,2,3)])
      
      lk_matrix[day_index, c("beta_0", "beta_x", "beta_y", 
                             "beta_gc", "beta_xy", "beta_xgc", "beta_ygc",
                             "beta_elev")] <-
        as.vector(lk_fit$d.coef)
      
 
      # Record quantiles of predictions on validation data
      lk_matrix[day_index, c("p_min", "p_p01", "p_p02", "p_p98", "p_p99", "p_max")] <- 
        c(min(lk_pred),
          quantile(lk_pred, probs = c(0.01, 0.02, 0.98, 0.99)),
          max(lk_pred))

       lk_matrix[day_index, c("f_min", "f_p01", "f_p02", "f_p98", "f_p99", "f_max")] <- 
        c(min(lk_fullpred),
          quantile(lk_fullpred, probs = c(0.01, 0.02, 0.98, 0.99)),
          max(lk_fullpred))

      
      # Save predictions for mapping
      ppredict_mat[, "lk_predict"] <- lk_pred
      
      # Generate a full prediction and save
      fpredict_mat[, "lk_predict"] <- lk_fullpred
      
    }  else {
      print(paste0("On day ", day, "there was an error"))
    }
      
    # Save prediction matrices
    ppredict_list[[day_name]] <- ppredict_mat
    fpredict_list[[day_name]] <- fpredict_mat
    
    rm(x, cov_names, day_index, day_name, 
       df_test, df_train, fpredict_mat, ppredict_mat,
       xy_train, Z_train, y_train, 
       xy_test, Z_test, y_test, 
       xy_full, Z_full, y_full,
       lk_time, lk_fit, lk_pred, lk_pred_time, lk_fullpred,
       holdout_index, full_id, test_id, train_id
    )
    
  }
  
  lk_return_list <- list(lk_matrix = lk_matrix, 
                         grid_info = grid_info,
                         ppredict_list = ppredict_list,
                         fpredict_list = fpredict_list)
  return(lk_return_list)
}

# Test run
# debugonce(fit_lk)
# test_1 <- fit_lk(days = 182L,
#                  NC.buffer = 5,
#                  a.wght = 12, nu = 0.1, nlevel = 5, NC = 15)


lk_final <- fit_lk(days = 182:212,
                   NC.buffer = 5,  
                   a.wght = 12, nu = 0.1, nlevel = 5, NC = 15)
saveRDS(lk_final, "./aod_holdouts/finalpred/lk_final_elev.rds")
rm(lk_final)
gc()


