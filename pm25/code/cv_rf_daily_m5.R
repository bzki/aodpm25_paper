# Daily model -- additional model M5
# See cv_rf_daily.R for the original details.
# M5 does the following:
# -- train solely on areas where AOD is actually observed
# -- predict everywhere by combining observed + imputed AOD

# Notes
# - Be careful about calculating convolution PM2.5 here
# - Create models 5a/5b/5c -- focus is on (b) model.
# - A single set of features considered only
# - Restrict to nodesize = 5 and four different mtry values (4, 8, 12, 16)

library(dplyr)
library(ranger)
library(fields)
sessionInfo()
source("./functions/make_convolution.R", echo = T, max.deparse.length = Inf)

if (!dir.exists("./pm25/cvresults/")) {
  dir.create("./pm25/cvresults/")
}

# Load daily PM2.5 data and full AOD predictions
df_list <- readRDS("./pm25/daily_cv.rds")$df_list
x_list <- readRDS("./pm25/daily_cv.rds")
cv_fold_list <- x_list$cv_folds
b1_fold_list <- x_list$b1_folds
b2_fold_list <- x_list$b2_folds
rm(x_list)
pred_list <- readRDS("./aod_holdouts/finalpred/allpreds.rds")

# Each model may suggest a slightly different mtry
# The models run are somewhat complicated by this fact. 
mtry_options <- c(4, 8, 12, 16)
node_options <- c(5)
param <- expand.grid(mtry = mtry_options, min.node.size = node_options) %>%
  arrange(mtry, min.node.size) %>%
  mutate(pnum = seq(nrow(.)))
rm(mtry_options, node_options)
NUM.TREES = 2000

# Output
days <- 182:212
cvdf_list = vector("list", length(days))
names(cvdf_list) <- paste0("day", days)
convcv_list = convb1_list = convb2_list = b1df_list = b2df_list = cvdf_list

for (day in days) {
  print(day)
  day_name <- paste0("day", day)
  
  cv_folds <- cv_fold_list[[day_name]]
  b1_folds <- b1_fold_list[[day_name]]
  b2_folds <- b2_fold_list[[day_name]]
  
  # Read in observed PM2.5 points, link with predictions
  df <- df_list[[day_name]] %>%
    dplyr::select(-year, -day, -aod_missing) %>%
    left_join(pred_list[[day_name]] %>% select(cmaq_id, sl3_pred),
              by = c("cmaq_id" = "cmaq_id")) %>%
    rename(aod_imputed = sl3_pred) %>%
    # Create a missing indicator and an AOD/GC Combined variable
    mutate(
      aod_missing = ifelse(is.na(aod_value), 1, 0),
      aod_gc_combine = ifelse(is.na(aod_value), gc_aod, aod_value),
      x_km = cmaq_x/1000,
      y_km = cmaq_y/1000
    )
  
  # Convolution layer for PM2.5 -- separately for each fold
  conv_cv <- 
    make_convolution(df_input = df, x_name = "x_km", y_name = "y_km", 
                     var_name = "pm25_value", fold_list = cv_folds)
  colnames(conv_cv) <- paste0("pm25_f", seq(length(cv_folds)))
  conv_b1 <- 
    make_convolution(df_input = df, x_name = "x_km", y_name = "y_km", 
                     var_name = "pm25_value", fold_list = b1_folds)
  colnames(conv_b1) <- paste0("pm25_f", seq(length(b1_folds)))
  conv_b2 <- 
    make_convolution(df_input = df, x_name = "x_km", y_name = "y_km", 
                     var_name = "pm25_value", fold_list = b2_folds)
  colnames(conv_b2) <- paste0("pm25_f", seq(length(b2_folds)))
  # Save these for potential re-analysis/checking
  convcv_list[[day_name]] <- conv_cv
  convb1_list[[day_name]] <- conv_b1
  convb2_list[[day_name]] <- conv_b2
  
  # Create data.frame in which we will save all predictions
  # -- either one overall or one for each of cv/b1/b2
  df_save_cv <- df %>% select(cmaq_id, cmaq_x, cmaq_y, pm25_value, aod_value,
                              aod_missing, aod_imputed, gc_aod, fold_cv) %>%
    mutate(day = day)
  df_save_b1 <- df %>% select(cmaq_id, cmaq_x, cmaq_y, pm25_value, aod_value,
                              aod_missing, aod_imputed, gc_aod, fold_b1) %>%
    mutate(day = day)
  df_save_b2 <- df %>% select(cmaq_id, cmaq_x, cmaq_y, pm25_value, aod_value,
                              aod_missing, aod_imputed, gc_aod, fold_b2) %>%
    mutate(day = day)
  
  # Create simplified data.frame
  df_5a <- df %>% select(-x_km, -y_km, -fold_cv, -fold_b1, -fold_b2,
                        -aod_gc_combine, -aod_missing)
  print(names(df_5a))
  
  for (i in seq(nrow(param))) {
    print(i)
    df_save_cv[, paste0("p5a_", i)] <- NA
    df_save_cv[, paste0("p5b_", i)] <- NA
    df_save_cv[, paste0("p5c_", i)] <- NA
    
    for (fi in seq(length(cv_folds))) {
      print(paste("fold number", fi)) 
      
      # Create training subsets for models
      train_5a <- df_5a[cv_folds[[fi]][[1]], ] %>%
        filter(!is.na(aod_value)) %>%
        select(-aod_imputed)
        
      test_5a <- df_5a[cv_folds[[fi]][[2]], ] %>%
        mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
        select(-aod_imputed)
      
      pm25_conv <- conv_cv[, fi]
      df_5b <- cbind(df_5a, pm25_conv)
      
      train_5b <- df_5b[cv_folds[[fi]][[1]], ] %>%
        filter(!is.na(aod_value)) %>%
        select(-aod_imputed)
      
      test_5b <- df_5b[cv_folds[[fi]][[2]], ] %>%
        mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
        select(-aod_imputed)
      
      testthat::expect_identical(names(train_5b), names(test_5b))
      testthat::expect_identical(names(train_5a), names(test_5a))
      testthat::expect_identical(nrow(train_5a), nrow(train_5b))
      
      set.seed(day + fi*5)
      fit_5a <- ranger(pm25_value ~ ., data = train_5a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_cv[cv_folds[[fi]][[2]], paste0("p5a_", i)] <- 
        predict(fit_5a, test_5a)$predictions
      rm(fit_5a)
      
      set.seed(day + fi*5)
      fit_5b <- ranger(pm25_value ~ ., data = train_5b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_cv[cv_folds[[fi]][[2]], paste0("p5b_", i)] <- 
        predict(fit_5b, test_5b)$predictions
      rm(fit_5b)
      
      # Uses same training data as 5b
      set.seed(day + fi*5)
      fit_5c <- ranger(pm25_value ~ ., data = train_5b, 
                       always.split.variables = c("cmaq_x", "cmaq_y"),
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_cv[cv_folds[[fi]][[2]], paste0("p5c_", i)] <- 
        predict(fit_5c, test_5b)$predictions
      rm(fit_5c)
      rm(df_5b, train_5a, train_5b, test_5a, test_5b, pm25_conv)
      
    }
  }
  rm(i, fi)
  
  for (i in seq(nrow(param))) {
    print(i)
    df_save_b1[, paste0("p5a_", i)] <- NA
    df_save_b1[, paste0("p5b_", i)] <- NA
    df_save_b1[, paste0("p5c_", i)] <- NA
    
    for (fi in seq(length(b1_folds))) {
      print(paste("fold number", fi)) 
      
      # Create training subsets for models
      train_5a <- df_5a[b1_folds[[fi]][[1]], ] %>%
        filter(!is.na(aod_value)) %>%
        select(-aod_imputed)
      
      test_5a <- df_5a[b1_folds[[fi]][[2]], ] %>%
        mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
        select(-aod_imputed)
      
      pm25_conv <- conv_b1[, fi]
      df_5b <- cbind(df_5a, pm25_conv)
      
      train_5b <- df_5b[b1_folds[[fi]][[1]], ] %>%
        filter(!is.na(aod_value)) %>%
        select(-aod_imputed)
      
      test_5b <- df_5b[b1_folds[[fi]][[2]], ] %>%
        mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
        select(-aod_imputed)
      
      testthat::expect_identical(names(train_5b), names(test_5b))
      testthat::expect_identical(names(train_5a), names(test_5a))
      testthat::expect_identical(nrow(train_5a), nrow(train_5b))
      
      set.seed(day + fi*5)
      fit_5a <- ranger(pm25_value ~ ., data = train_5a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_b1[b1_folds[[fi]][[2]], paste0("p5a_", i)] <- 
        predict(fit_5a, test_5a)$predictions
      rm(fit_5a)
      
      set.seed(day + fi*5)
      fit_5b <- ranger(pm25_value ~ ., data = train_5b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_b1[b1_folds[[fi]][[2]], paste0("p5b_", i)] <- 
        predict(fit_5b, test_5b)$predictions
      rm(fit_5b)
      
      # Uses same training data as 5b
      set.seed(day + fi*5)
      fit_5c <- ranger(pm25_value ~ ., data = train_5b, 
                       always.split.variables = c("cmaq_x", "cmaq_y"),
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_b1[b1_folds[[fi]][[2]], paste0("p5c_", i)] <- 
        predict(fit_5c, test_5b)$predictions
      rm(fit_5c)
      rm(df_5b, train_5a, train_5b, test_5a, test_5b, pm25_conv)
      
    }
  }
  rm(i, fi) 
  
  for (i in seq(nrow(param))) {
    print(i)
    df_save_b2[, paste0("p5a_", i)] <- NA
    df_save_b2[, paste0("p5b_", i)] <- NA
    df_save_b2[, paste0("p5c_", i)] <- NA
    
    for (fi in seq(length(b2_folds))) {
      print(paste("fold number", fi)) 
      
      # Create training subsets for models
      train_5a <- df_5a[b2_folds[[fi]][[1]], ] %>%
        filter(!is.na(aod_value)) %>%
        select(-aod_imputed)
      
      test_5a <- df_5a[b2_folds[[fi]][[2]], ] %>%
        mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
        select(-aod_imputed)
      
      pm25_conv <- conv_b2[, fi]
      df_5b <- cbind(df_5a, pm25_conv)
      
      train_5b <- df_5b[b2_folds[[fi]][[1]], ] %>%
        filter(!is.na(aod_value)) %>%
        select(-aod_imputed)
      
      test_5b <- df_5b[b2_folds[[fi]][[2]], ] %>%
        mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
        select(-aod_imputed)
      
      testthat::expect_identical(names(train_5b), names(test_5b))
      testthat::expect_identical(names(train_5a), names(test_5a))
      testthat::expect_identical(nrow(train_5a), nrow(train_5b))
      
      set.seed(day + fi*5)
      fit_5a <- ranger(pm25_value ~ ., data = train_5a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_b2[b2_folds[[fi]][[2]], paste0("p5a_", i)] <- 
        predict(fit_5a, test_5a)$predictions
      rm(fit_5a)
      
      set.seed(day + fi*5)
      fit_5b <- ranger(pm25_value ~ ., data = train_5b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_b2[b2_folds[[fi]][[2]], paste0("p5b_", i)] <- 
        predict(fit_5b, test_5b)$predictions
      rm(fit_5b)
      
      # Uses same training data as 5b
      set.seed(day + fi*5)
      fit_5c <- ranger(pm25_value ~ ., data = train_5b, 
                       always.split.variables = c("cmaq_x", "cmaq_y"),
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      df_save_b2[b2_folds[[fi]][[2]], paste0("p5c_", i)] <- 
        predict(fit_5c, test_5b)$predictions
      rm(fit_5c)
      rm(df_5b, train_5a, train_5b, test_5a, test_5b, pm25_conv)
      
    }
  }
  rm(i, fi) 
  
  cvdf_list[[day_name]] <- df_save_cv
  b1df_list[[day_name]] <- df_save_b1
  b2df_list[[day_name]] <- df_save_b2
  
  rm(cv_folds, b1_folds, b2_folds, df, 
     df_save_cv, df_save_b1, df_save_b2,
     df_5a, conv_cv, conv_b1, conv_b2
  )
  gc()
}

comb_list <- list(
  cvdf_list = cvdf_list,
  b1df_list = b1df_list,
  b2df_list = b2df_list,
  param = param,
  convcv_list = convcv_list,
  convb1_list = convb1_list,
  convb2_list = convb2_list
)

saveRDS(comb_list, "./pm25/cvresults/cv_rf_results_m5.rds")
