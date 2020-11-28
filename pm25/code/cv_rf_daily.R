## 4 models for daily random forest  models ##
# M1: Without AOD and without GC
# M2: Without AOD but with GC
# M3: Imputed AOD (and with GC)
# M4: AOD observed + GC where AOD missing
# (a) models do not include convolution layer for PM2.5
# (b): Include convolution layer for PM2.5 (Weighted average of PM2.5)
# The M(b) models are of main interest. 
# nodesize is limited to 5 -- previous experience suggests very little changes
# mtry values of 4, 8, 12, 16 are used. 
# 2000 trees are fit. 

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
  df_s <- df %>% select(-x_km, -y_km, -fold_cv, -fold_b1, -fold_b2)
  
  # Create data split
  # No AOD, NO GC
  df_1a <- df_s %>%
    select(-cmaq_id, -aod_value, -aod_imputed, -gc_aod, -aod_gc_combine)
  print(names(df_1a))
  
  # Without AOD but with GC,
  df_2a <- df_s %>%
    select(-cmaq_id, -aod_value, -aod_imputed, -aod_gc_combine)
  print(names(df_2a))
  
  # Imputed AOD (leaving in GC)
  df_3a <- df_s %>%
    select(-cmaq_id, -aod_value, -aod_gc_combine)
  print(names(df_3a))
  
  # AOD observed + GC where AOD missing
  df_4a <- df_s %>%
    select(-cmaq_id, -aod_value, -aod_imputed, -gc_aod)
  print(names(df_4a))
  
  for (i in seq(nrow(param))) {
    print(i)
    df_save_cv[, paste0("p", seq(4), "a_", i)] <- NA
    df_save_cv[, paste0("p", seq(4), "b_", i)] <- NA
    
    for (fi in seq(length(cv_folds))) {
      
      print(paste("fold number", fi)) 
      
      # Create training subsets for models
      train_1a <- df_1a[cv_folds[[fi]][[1]], ]
      train_2a <- df_2a[cv_folds[[fi]][[1]], ]
      train_3a <- df_3a[cv_folds[[fi]][[1]], ]
      train_4a <- df_4a[cv_folds[[fi]][[1]], ]
  
      test_1a <- df_1a[cv_folds[[fi]][[2]], ]
      test_2a <- df_2a[cv_folds[[fi]][[2]], ]
      test_3a <- df_3a[cv_folds[[fi]][[2]], ]
      test_4a <- df_4a[cv_folds[[fi]][[2]], ]
      
      pm25_conv <- conv_cv[, fi]
      df_1b <- cbind(df_1a, pm25_conv)
      df_2b <- cbind(df_2a, pm25_conv)
      df_3b <- cbind(df_3a, pm25_conv)
      df_4b <- cbind(df_4a, pm25_conv)
      
      train_1b <- df_1b[cv_folds[[fi]][[1]], ]
      train_2b <- df_2b[cv_folds[[fi]][[1]], ]
      train_3b <- df_3b[cv_folds[[fi]][[1]], ]
      train_4b <- df_4b[cv_folds[[fi]][[1]], ]
      
      test_1b <- df_1b[cv_folds[[fi]][[2]], ]
      test_2b <- df_2b[cv_folds[[fi]][[2]], ]
      test_3b <- df_3b[cv_folds[[fi]][[2]], ]
      test_4b <- df_4b[cv_folds[[fi]][[2]], ]
      
      set.seed(day + fi*1)
      fit_1a <- ranger(pm25_value ~ ., data = train_1a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*2)
      fit_2a <- ranger(pm25_value ~ ., data = train_2a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*3)
      fit_3a <- ranger(pm25_value ~ ., data = train_3a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*4)
      fit_4a <- ranger(pm25_value ~ ., data = train_4a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
  
      set.seed(day + fi*1)
      fit_1b <- ranger(pm25_value ~ ., data = train_1b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*2)
      fit_2b <- ranger(pm25_value ~ ., data = train_2b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*3)
      fit_3b <- ranger(pm25_value ~ ., data = train_3b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*4)
      fit_4b <- ranger(pm25_value ~ ., data = train_4b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      
      ## Get predictions for each parameter combination of interest
      df_save_cv[cv_folds[[fi]][[2]], paste0("p1a_", i)] <- 
        predict(fit_1a, test_1a)$predictions
      df_save_cv[cv_folds[[fi]][[2]], paste0("p2a_", i)] <- 
        predict(fit_2a, test_2a)$predictions
      df_save_cv[cv_folds[[fi]][[2]], paste0("p3a_", i)] <- 
        predict(fit_3a, test_3a)$predictions
      df_save_cv[cv_folds[[fi]][[2]], paste0("p4a_", i)] <- 
        predict(fit_4a, test_4a)$predictions
      
      df_save_cv[cv_folds[[fi]][[2]], paste0("p1b_", i)] <- 
        predict(fit_1b, test_1b)$predictions
      df_save_cv[cv_folds[[fi]][[2]], paste0("p2b_", i)] <- 
        predict(fit_2b, test_2b)$predictions
      df_save_cv[cv_folds[[fi]][[2]], paste0("p3b_", i)] <- 
        predict(fit_3b, test_3b)$predictions
      df_save_cv[cv_folds[[fi]][[2]], paste0("p4b_", i)] <- 
        predict(fit_4b, test_4b)$predictions
    }
  }
  rm(i, fi)
  
  for (i in seq(nrow(param))) {
    print(i)
    df_save_b1[, paste0("p", seq(4), "a_", i)] <- NA
    df_save_b1[, paste0("p", seq(4), "b_", i)] <- NA
    
    for (fi in seq(length(b1_folds))) {
      print(paste("fold number", fi)) 
      
      # Create training subsets for models
      train_1a <- df_1a[b1_folds[[fi]][[1]], ]
      train_2a <- df_2a[b1_folds[[fi]][[1]], ]
      train_3a <- df_3a[b1_folds[[fi]][[1]], ]
      train_4a <- df_4a[b1_folds[[fi]][[1]], ]
      
      test_1a <- df_1a[b1_folds[[fi]][[2]], ]
      test_2a <- df_2a[b1_folds[[fi]][[2]], ]
      test_3a <- df_3a[b1_folds[[fi]][[2]], ]
      test_4a <- df_4a[b1_folds[[fi]][[2]], ]
      
      pm25_conv <- conv_b1[, fi]
      df_1b <- cbind(df_1a, pm25_conv)
      df_2b <- cbind(df_2a, pm25_conv)
      df_3b <- cbind(df_3a, pm25_conv)
      df_4b <- cbind(df_4a, pm25_conv)
      
      train_1b <- df_1b[b1_folds[[fi]][[1]], ]
      train_2b <- df_2b[b1_folds[[fi]][[1]], ]
      train_3b <- df_3b[b1_folds[[fi]][[1]], ]
      train_4b <- df_4b[b1_folds[[fi]][[1]], ]
      
      test_1b <- df_1b[b1_folds[[fi]][[2]], ]
      test_2b <- df_2b[b1_folds[[fi]][[2]], ]
      test_3b <- df_3b[b1_folds[[fi]][[2]], ]
      test_4b <- df_4b[b1_folds[[fi]][[2]], ]
      
      set.seed(day + fi*1)
      fit_1a <- ranger(pm25_value ~ ., data = train_1a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*2)
      fit_2a <- ranger(pm25_value ~ ., data = train_2a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*3)
      fit_3a <- ranger(pm25_value ~ ., data = train_3a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*4)
      fit_4a <- ranger(pm25_value ~ ., data = train_4a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      
      set.seed(day + fi*1)
      fit_1b <- ranger(pm25_value ~ ., data = train_1b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*2)
      fit_2b <- ranger(pm25_value ~ ., data = train_2b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*3)
      fit_3b <- ranger(pm25_value ~ ., data = train_3b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*4)
      fit_4b <- ranger(pm25_value ~ ., data = train_4b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      
      ## Get predictions for each parameter combination of interest
      df_save_b1[b1_folds[[fi]][[2]], paste0("p1a_", i)] <- 
        predict(fit_1a, test_1a)$predictions
      df_save_b1[b1_folds[[fi]][[2]], paste0("p2a_", i)] <- 
        predict(fit_2a, test_2a)$predictions
      df_save_b1[b1_folds[[fi]][[2]], paste0("p3a_", i)] <- 
        predict(fit_3a, test_3a)$predictions
      df_save_b1[b1_folds[[fi]][[2]], paste0("p4a_", i)] <- 
        predict(fit_4a, test_4a)$predictions
      
      df_save_b1[b1_folds[[fi]][[2]], paste0("p1b_", i)] <- 
        predict(fit_1b, test_1b)$predictions
      df_save_b1[b1_folds[[fi]][[2]], paste0("p2b_", i)] <- 
        predict(fit_2b, test_2b)$predictions
      df_save_b1[b1_folds[[fi]][[2]], paste0("p3b_", i)] <- 
        predict(fit_3b, test_3b)$predictions
      df_save_b1[b1_folds[[fi]][[2]], paste0("p4b_", i)] <- 
        predict(fit_4b, test_4b)$predictions
    }
  }
  rm(i, fi)
  
  for (i in seq(nrow(param))) {
    print(i)
    df_save_b2[, paste0("p", seq(4), "a_", i)] <- NA
    df_save_b2[, paste0("p", seq(4), "b_", i)] <- NA
    
    for (fi in seq(length(b2_folds))) {
      print(paste("fold number", fi)) 
      
      # Create training subsets for models
      train_1a <- df_1a[b2_folds[[fi]][[1]], ]
      train_2a <- df_2a[b2_folds[[fi]][[1]], ]
      train_3a <- df_3a[b2_folds[[fi]][[1]], ]
      train_4a <- df_4a[b2_folds[[fi]][[1]], ]
      
      test_1a <- df_1a[b2_folds[[fi]][[2]], ]
      test_2a <- df_2a[b2_folds[[fi]][[2]], ]
      test_3a <- df_3a[b2_folds[[fi]][[2]], ]
      test_4a <- df_4a[b2_folds[[fi]][[2]], ]
      
      pm25_conv <- conv_b2[, fi]
      df_1b <- cbind(df_1a, pm25_conv)
      df_2b <- cbind(df_2a, pm25_conv)
      df_3b <- cbind(df_3a, pm25_conv)
      df_4b <- cbind(df_4a, pm25_conv)
      
      train_1b <- df_1b[b2_folds[[fi]][[1]], ]
      train_2b <- df_2b[b2_folds[[fi]][[1]], ]
      train_3b <- df_3b[b2_folds[[fi]][[1]], ]
      train_4b <- df_4b[b2_folds[[fi]][[1]], ]
      
      test_1b <- df_1b[b2_folds[[fi]][[2]], ]
      test_2b <- df_2b[b2_folds[[fi]][[2]], ]
      test_3b <- df_3b[b2_folds[[fi]][[2]], ]
      test_4b <- df_4b[b2_folds[[fi]][[2]], ]
      
      set.seed(day + fi*1)
      fit_1a <- ranger(pm25_value ~ ., data = train_1a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*2)
      fit_2a <- ranger(pm25_value ~ ., data = train_2a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*3)
      fit_3a <- ranger(pm25_value ~ ., data = train_3a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*4)
      fit_4a <- ranger(pm25_value ~ ., data = train_4a, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      
      set.seed(day + fi*1)
      fit_1b <- ranger(pm25_value ~ ., data = train_1b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*2)
      fit_2b <- ranger(pm25_value ~ ., data = train_2b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*3)
      fit_3b <- ranger(pm25_value ~ ., data = train_3b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      set.seed(day + fi*4)
      fit_4b <- ranger(pm25_value ~ ., data = train_4b, 
                       num.trees = NUM.TREES, mtry = param$mtry[i], 
                       min.node.size = param$min.node.size[i])
      
      ## Get predictions for each parameter combination of interest
      df_save_b2[b2_folds[[fi]][[2]], paste0("p1a_", i)] <- 
        predict(fit_1a, test_1a)$predictions
      df_save_b2[b2_folds[[fi]][[2]], paste0("p2a_", i)] <- 
        predict(fit_2a, test_2a)$predictions
      df_save_b2[b2_folds[[fi]][[2]], paste0("p3a_", i)] <- 
        predict(fit_3a, test_3a)$predictions
      df_save_b2[b2_folds[[fi]][[2]], paste0("p4a_", i)] <- 
        predict(fit_4a, test_4a)$predictions
      
      df_save_b2[b2_folds[[fi]][[2]], paste0("p1b_", i)] <- 
        predict(fit_1b, test_1b)$predictions
      df_save_b2[b2_folds[[fi]][[2]], paste0("p2b_", i)] <- 
        predict(fit_2b, test_2b)$predictions
      df_save_b2[b2_folds[[fi]][[2]], paste0("p3b_", i)] <- 
        predict(fit_3b, test_3b)$predictions
      df_save_b2[b2_folds[[fi]][[2]], paste0("p4b_", i)] <- 
        predict(fit_4b, test_4b)$predictions
    }
  }
  rm(i, fi)
  
  cvdf_list[[day_name]] <- df_save_cv
  b1df_list[[day_name]] <- df_save_b1
  b2df_list[[day_name]] <- df_save_b2
  
  saveRDS(list(df_save_cv = df_save_cv, 
               df_save_b1 = df_save_b1,
               df_save_b2 = df_save_b2),
    paste0("./pm25/cvresults/cv_rf_", day, ".rds"))
  
  rm(cv_folds, b1_folds, b2_folds, df, 
     fit_1a, fit_2a, fit_3a, fit_4a, fit_1b, fit_2b, fit_3b, fit_4b,
     df_save_cv, df_save_b1, df_save_b2,
     df_1a, df_2a, df_3a, df_4a, df_1b, df_2b, df_3b, df_4b,
     train_1a, train_2a, train_3a, train_4a, train_1b, train_2b, train_3b, train_4b,
     test_1a, test_2a, test_3a, test_4a, test_1b, test_2b, test_3b, test_4b,
     conv_cv, conv_b1, conv_b2
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

saveRDS(comb_list, "./pm25/cvresults/cv_rf_results.rds")
