# Run on full observed data cross-validation folds
# Predict 10 spatially clustered cross-validation folds with LatticeKrig

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


fit_cv_lk_full <- function(days = 182:212, 
                           prefix_read = "./rawrds/split_fcv_",
                           prefix_save = "./aod_holdouts/cvdata_full/lk_fcv_",
                           NC.buffer = 5, a.wght = 12, nu = 0.1, 
                           nlevel = 5, NC = 15) {
        
  
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

    # Read in training data
    train_list <- readRDS(paste0(prefix_read, day, ".rds"))
    x <- train_list$df
    folds <- train_list$folds
    rm(train_list)

    # Covariates used for mean model
    # spatial drift (cmaq_x, cmaq_y) included by default
    # More: gc_aod + cmaq_x*cmaq_y 
    #       + cmaq_x * gc_aod + cmaq_y * gc_aod
    #       + elevation (normalized)
    
    x_simple <- x %>%
      dplyr::select(cmaq_id, foldID, cv_nndist,
                    aod_value, cmaq_x, cmaq_y, gc_aod, elev) %>%
      # For numerical stability, putting coordinates in terms of KM
      mutate(x_km = cmaq_x/1000, y_km = cmaq_y/1000,
             cv_nndist_km = cv_nndist/1000) %>%
      dplyr::select(-cmaq_x, -cmaq_y, -cv_nndist) %>%
      mutate(
        xy_int = x_km * y_km,
        x_gc_aod = x_km*gc_aod,
        y_gc_aod = y_km*gc_aod,
        elev_norm = (elev - mean(elev))/sd(elev)
      ) %>%
      dplyr::select(-elev)
    
    x_basic <- x_simple %>%
      dplyr::select(cmaq_id, foldID, x_km, y_km, cv_nndist_km, aod_value)
    x_basic$lk_cv <- NA
    
    # Covariate names should exclud x/y -- included by default in LK fits
    cov_names <- 
      names(x_simple)[!(names(x_simple) %in% 
                          c("x_km", "y_km", "aod_value", 
                            "foldID", "cmaq_id", "cv_nndist_km"))]
    
    y_value <- x_simple$aod_value
    xy_matrix <- as.matrix(x_simple[, c("x_km", "y_km")])
    Z_matrix <- as.matrix(x_simple[, cov_names])
    
    rm(x_simple, x)
    
    num_folds <- length(folds)
    
    for (fi in seq(num_folds)) {
      
      print(paste0("Fold ", fi))
      
      y_train <- y_value[folds[[fi]][[1]]]
      xy_train <- xy_matrix[folds[[fi]][[1]], ]
      Z_train <- Z_matrix[folds[[fi]][[1]], ]
      
      xy_cv <- xy_matrix[folds[[fi]][[2]], ]
      Z_cv <- Z_matrix[folds[[fi]][[2]], ]
      
      # Ensure that the fold is correctly encoded
      testthat::expect_true(
        dplyr::setequal(which(x_basic$foldID == fi), 
                        folds[[fi]][[2]])
      )
      
      # Fit LatticeKrig
      lk_fit <- LatticeKrig(x = xy_train, y = y_train, 
                            Z = Z_train, LKinfo = lk_info)
      
      # Predict
      lk_pred <-  as.vector(predict.LKrig(lk_fit, xnew = xy_cv, Znew = Z_cv))
      
      # Save the predictions correctly
      x_basic$lk_cv[folds[[fi]][[2]]] <- lk_pred
      
      rm(y_train, xy_train, Z_train, xy_cv, Z_cv, lk_pred, lk_fit)
      
    }
    
    testthat::expect_equal(0, sum(is.na(x_basic$lk_cv)))
    
    # Save x_basic
    saveRDS(x_basic, 
            paste0(prefix_save, day, ".rds"))
    
    rm(x_basic, y_value, xy_matrix, Z_matrix, num_folds, folds, fi)
    
  }
  
  rm(lk_info)
  return(list("Success"))
}



#### Fit ####
# debugonce(fit_cv_lk_full)
# fit_cv_lk_full(days = 182L, 
#                prefix_read = "./rawrds/split_fcv_",
#                prefix_save = "./aod_holdouts/cvdata_full/lk_fcv_",
#                NC.buffer = 5, a.wght = 12, nu = 0.1, 
#                nlevel = 5, NC = 15)

fit_cv_lk_full(days = 182:212, 
               prefix_read = "./rawrds/split_fcv_",
               prefix_save = "./aod_holdouts/cvdata_full/lk_fcv_",
               NC.buffer = 5, a.wght = 12, nu = 0.1, 
               nlevel = 5, NC = 15)

gc()