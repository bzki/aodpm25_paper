# Fit overall (i.e. spatiotemporal) ranger model to entire data 
#  using 10-fold cross-validation methods.

# 1: Random CV (1)
# 2: Spatial cluster (1)
# 3: Spatial cluster (2)

# Just et al 2018 include days since 1970
# Hu et al. 2017 include dummy day variables. 
# Staffogia et al. 2019 include day of the year (integer) and day of the week

# My approach: 
# -- day (as integer), 
# -- day of week (as integer) -- can also include indicator for day of week. 

# - Try just 4 parameters: mtry = 4, 8, 12, 16 for nodesize = 5. 
# Models to run
# - 4 models with convolution PM2.5, using RF
# Not run but previously tried
# - 4 models with convolution PM2.5, using extratrees
# - 4 Random forest forcing (x,y)
# - 4 Random forest forcing (day, day of week)

# Read in the convolution PM2.5
library(dplyr)
library(ranger)
sessionInfo()

# Read data.
df_list <- readRDS("./pm25/daily_cv.rds")$df_list
df_stack <- do.call("rbind", df_list)
pred_list <- readRDS("./aod_holdouts/finalpred/allpreds.rds")
pred_stack <- do.call("rbind", pred_list) %>%
  select(cmaq_id, day, sl3_pred)
df_stack2 <- left_join(df_stack, pred_stack,
                       by = c("cmaq_id" = "cmaq_id", "day" = "day")) %>%
  select(-aod_missing, -year) %>%
  rename(aod_imputed = sl3_pred) %>%
  # Create a missing indicator and an AOD/GC Combined variable
  mutate(
    aod_missing = ifelse(is.na(aod_value), 1, 0),
    aod_gc_combine = ifelse(is.na(aod_value), gc_aod, aod_value)
  ) %>%
  # Create day of the week variable
  mutate(
    # Friday is 5; Saturday is 6; Sunday is 0; etc.
    # Day 182 is July 1, 2011
    day_of_week = (day - 177L) %% 7L
  )
rm(df_list, df_stack, pred_list, pred_stack)

# Reconstruct fold lists for cv/b1/b2 (and any others in the future)
cv_folds = b1_folds = b2_folds = vector("list", 10L)
for (fi in seq(length(cv_folds))) {
  cv_folds[[fi]][[1]] <- which(df_stack2$fold_cv != fi)
  cv_folds[[fi]][[2]] <- which(df_stack2$fold_cv == fi)
  
  b1_folds[[fi]][[1]] <- which(df_stack2$fold_b1 != fi)
  b1_folds[[fi]][[2]] <- which(df_stack2$fold_b1 == fi)
  
  b2_folds[[fi]][[1]] <- which(df_stack2$fold_b2 != fi)
  b2_folds[[fi]][[2]] <- which(df_stack2$fold_b2 == fi)
}
rm(fi)

# Prepare and stack the convolution PM2.5 created in previous script
# (see cv_rf_pm25.R)
comb_list <- readRDS("./pm25/cvresults/cv_rf_results.rds")
conv_cv <- comb_list$convcv_list
conv_b1 <- comb_list$convb1_list
conv_b2 <- comb_list$convb2_list
rm(comb_list)
days = 182:212
conv_cv_list <- vector("list", length(days))
names(conv_cv_list) <- paste0("day", days)
conv_b1_list = conv_b2_list = conv_cv_list

for (day in days) {
  day_name <- paste0("day", day)
  cv_df <- as.data.frame(conv_cv[[day_name]])
  cv_df[, "day"] <- day
  conv_cv_list[[day_name]] <- cv_df

  b1_df <- as.data.frame(conv_b1[[day_name]])
  b1_df[, "day"] <- day
  conv_b1_list[[day_name]] <- b1_df

  b2_df <- as.data.frame(conv_b2[[day_name]])
  b2_df[, "day"] <- day
  conv_b2_list[[day_name]] <- b2_df
}
rm(cv_df, b1_df, b2_df)
rm(day_name, day)
conv_cv_stack <- do.call("rbind", conv_cv_list)
conv_b1_stack <- do.call("rbind", conv_b1_list)
conv_b2_stack <- do.call("rbind", conv_b2_list)
rm(conv_cv_list, conv_b1_list, conv_b2_list,
   conv_cv, conv_b1, conv_b2)

# Parameter information
mtry_options <- c(4, 8, 12, 16)
node_options <- c(5)
param <- expand.grid(mtry = mtry_options, min.node.size = node_options) %>%
  arrange(mtry, min.node.size) %>%
  mutate(pnum = seq(nrow(.)))
rm(mtry_options, node_options)
NUM.TREES = 2000

# Create a dataset that will have all of the predictions and be saved
df_save_cv <- df_stack2 %>% select(day, cmaq_id, cmaq_x, cmaq_y, pm25_value, 
                            aod_value, aod_missing, aod_imputed, fold_cv)
df_save_b1 <- df_stack2 %>% select(day, cmaq_id, cmaq_x, cmaq_y, pm25_value, 
                            aod_value, aod_missing, aod_imputed, fold_b1)
df_save_b2 <- df_stack2 %>% select(day, cmaq_id, cmaq_x, cmaq_y, pm25_value, 
                            aod_value, aod_missing, aod_imputed, fold_b2)

# Create templates for the 4 principal models 
df_s <- df_stack2 %>% select(-cmaq_id, -fold_cv, -fold_b1, -fold_b2, -aod_value)

# No AOD, NO GC
df_1 <- df_s %>%
  select(-aod_imputed, -gc_aod, -aod_gc_combine)

# Without AOD but with GC,
df_2 <- df_s %>%
  select(-aod_imputed, -aod_gc_combine)

# Imputed AOD (leaving in GC)
df_3 <- df_s %>%
  select(-aod_gc_combine)

# AOD observed + GC where AOD missing
df_4 <- df_s %>%
  select(-aod_imputed, -gc_aod)

# Loop through each CV, each mtry, each fold, and each model
for (i in seq(nrow(param))) {
  print(paste("param number", i))
  # a => random forest
  # b => extra trees
  # c => rf + keep (x,y) at every split
  # d => rf + keep day, day of week at every split
  df_save_cv[, paste0("p", seq(4), "a_", i)] <- NA
  df_save_cv[, paste0("p", seq(4), "b_", i)] <- NA
  df_save_cv[, paste0("p", seq(4), "c_", i)] <- NA
  df_save_cv[, paste0("p", seq(4), "d_", i)] <- NA
  
  for (fi in seq(length(cv_folds))) {
    
    print(paste("fold number", fi)) 

    # Assign the appropriate convolution PM2.5
    pm25_conv <- conv_cv_stack[, paste0("pm25_f", fi)]
    df_1f <- cbind(df_1, pm25_conv)
    df_2f <- cbind(df_2, pm25_conv)
    df_3f <- cbind(df_3, pm25_conv)
    df_4f <- cbind(df_4, pm25_conv)
    rm(pm25_conv)
    
    train_1 <- df_1f[cv_folds[[fi]][[1]], ]
    train_2 <- df_2f[cv_folds[[fi]][[1]], ]
    train_3 <- df_3f[cv_folds[[fi]][[1]], ]
    train_4 <- df_4f[cv_folds[[fi]][[1]], ]
    
    test_1 <- df_1f[cv_folds[[fi]][[2]], ]
    test_2 <- df_2f[cv_folds[[fi]][[2]], ]
    test_3 <- df_3f[cv_folds[[fi]][[2]], ]
    test_4 <- df_4f[cv_folds[[fi]][[2]], ]
    
    set.seed(fi*1)
    fit_1a <- ranger(pm25_value ~ ., data = train_1, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*2)
    fit_2a <- ranger(pm25_value ~ ., data = train_2, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*3)
    fit_3a <- ranger(pm25_value ~ ., data = train_3, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*4)
    fit_4a <- ranger(pm25_value ~ ., data = train_4, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    
    ## Get predictions for each parameter combination of interest
    df_save_cv[cv_folds[[fi]][[2]], paste0("p1a_", i)] <- 
      predict(fit_1a, test_1)$predictions
    df_save_cv[cv_folds[[fi]][[2]], paste0("p2a_", i)] <- 
      predict(fit_2a, test_2)$predictions
    df_save_cv[cv_folds[[fi]][[2]], paste0("p3a_", i)] <- 
      predict(fit_3a, test_3)$predictions
    df_save_cv[cv_folds[[fi]][[2]], paste0("p4a_", i)] <- 
      predict(fit_4a, test_4)$predictions

    rm(fit_1a, fit_2a, fit_3a, fit_4a)
    gc()
    
    # set.seed(fi*1)
    # fit_1c <- ranger(pm25_value ~ ., data = train_1, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*2)
    # fit_2c <- ranger(pm25_value ~ ., data = train_2, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*3)
    # fit_3c <- ranger(pm25_value ~ ., data = train_3, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*4)
    # fit_4c <- ranger(pm25_value ~ ., data = train_4, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # 
    # df_save_cv[cv_folds[[fi]][[2]], paste0("p1c_", i)] <- 
    #   predict(fit_1c, test_1)$predictions
    # df_save_cv[cv_folds[[fi]][[2]], paste0("p2c_", i)] <- 
    #   predict(fit_2c, test_2)$predictions
    # df_save_cv[cv_folds[[fi]][[2]], paste0("p3c_", i)] <- 
    #   predict(fit_3c, test_3)$predictions
    # df_save_cv[cv_folds[[fi]][[2]], paste0("p4c_", i)] <- 
    #   predict(fit_4c, test_4)$predictions
    # rm(fit_1c, fit_2c, fit_3c, fit_4c)
    # gc()

    rm(train_1, train_2, train_3, train_4,
       test_1, test_2, test_3, test_4,
       df_1f, df_2f, df_3f, df_4f)
    gc()
  }
}
rm(i, fi)
saveRDS(df_save_cv, "./pm25/cvresults/rftemp_cv.rds")

# Loop through each CV, each mtry, each fold, and each model
for (i in seq(nrow(param))) {
  print(paste("param number", i))
  # a => random forest
  # b => extra trees
  # c => rf + keep (x,y) at every split
  # d => rf + keep day, day of week at every split

  df_save_b1[, paste0("p", seq(4), "a_", i)] <- NA
  df_save_b1[, paste0("p", seq(4), "b_", i)] <- NA
  df_save_b1[, paste0("p", seq(4), "c_", i)] <- NA
  df_save_b1[, paste0("p", seq(4), "d_", i)] <- NA
  
  for (fi in seq(length(b1_folds))) {
    
    print(paste("fold number", fi)) 
    
    # Assign the appropriate convolution PM2.5
    pm25_conv <- conv_b1_stack[, paste0("pm25_f", fi)]
    df_1f <- cbind(df_1, pm25_conv)
    df_2f <- cbind(df_2, pm25_conv)
    df_3f <- cbind(df_3, pm25_conv)
    df_4f <- cbind(df_4, pm25_conv)
    rm(pm25_conv)
    
    train_1 <- df_1f[b1_folds[[fi]][[1]], ]
    train_2 <- df_2f[b1_folds[[fi]][[1]], ]
    train_3 <- df_3f[b1_folds[[fi]][[1]], ]
    train_4 <- df_4f[b1_folds[[fi]][[1]], ]
    
    test_1 <- df_1f[b1_folds[[fi]][[2]], ]
    test_2 <- df_2f[b1_folds[[fi]][[2]], ]
    test_3 <- df_3f[b1_folds[[fi]][[2]], ]
    test_4 <- df_4f[b1_folds[[fi]][[2]], ]
    
    set.seed(fi*1)
    fit_1a <- ranger(pm25_value ~ ., data = train_1, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*2)
    fit_2a <- ranger(pm25_value ~ ., data = train_2, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*3)
    fit_3a <- ranger(pm25_value ~ ., data = train_3, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*4)
    fit_4a <- ranger(pm25_value ~ ., data = train_4, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    
    ## Get predictions for each parameter combination of interest
    df_save_b1[b1_folds[[fi]][[2]], paste0("p1a_", i)] <- 
      predict(fit_1a, test_1)$predictions
    df_save_b1[b1_folds[[fi]][[2]], paste0("p2a_", i)] <- 
      predict(fit_2a, test_2)$predictions
    df_save_b1[b1_folds[[fi]][[2]], paste0("p3a_", i)] <- 
      predict(fit_3a, test_3)$predictions
    df_save_b1[b1_folds[[fi]][[2]], paste0("p4a_", i)] <- 
      predict(fit_4a, test_4)$predictions
    
    rm(fit_1a, fit_2a, fit_3a, fit_4a)
    gc()
    
    # set.seed(fi*1)
    # fit_1c <- ranger(pm25_value ~ ., data = train_1, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*2)
    # fit_2c <- ranger(pm25_value ~ ., data = train_2, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*3)
    # fit_3c <- ranger(pm25_value ~ ., data = train_3, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*4)
    # fit_4c <- ranger(pm25_value ~ ., data = train_4, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # 
    # df_save_b1[b1_folds[[fi]][[2]], paste0("p1c_", i)] <- 
    #   predict(fit_1c, test_1)$predictions
    # df_save_b1[b1_folds[[fi]][[2]], paste0("p2c_", i)] <- 
    #   predict(fit_2c, test_2)$predictions
    # df_save_b1[b1_folds[[fi]][[2]], paste0("p3c_", i)] <- 
    #   predict(fit_3c, test_3)$predictions
    # df_save_b1[b1_folds[[fi]][[2]], paste0("p4c_", i)] <- 
    #   predict(fit_4c, test_4)$predictions
    # rm(fit_1c, fit_2c, fit_3c, fit_4c)
    # gc()
    
    rm(train_1, train_2, train_3, train_4,
       test_1, test_2, test_3, test_4,
       df_1f, df_2f, df_3f, df_4f)
    gc()
  }
}
rm(i, fi)
saveRDS(df_save_b1, "./pm25/cvresults/rftemp_b1.rds")


# Loop through each CV, each mtry, each fold, and each model
for (i in seq(nrow(param))) {
  print(paste("param number", i))
  # a => random forest
  # b => extra trees
  # c => rf + keep (x,y) at every split
  # d => rf + keep day, day of week at every split
  df_save_b2[, paste0("p", seq(4), "a_", i)] <- NA
  df_save_b2[, paste0("p", seq(4), "b_", i)] <- NA
  df_save_b2[, paste0("p", seq(4), "c_", i)] <- NA
  df_save_b2[, paste0("p", seq(4), "d_", i)] <- NA
  
  for (fi in seq(length(b2_folds))) {
    
    print(paste("fold number", fi)) 
    
    # Assign the appropriate convolution PM2.5
    pm25_conv <- conv_b2_stack[, paste0("pm25_f", fi)]
    df_1f <- cbind(df_1, pm25_conv)
    df_2f <- cbind(df_2, pm25_conv)
    df_3f <- cbind(df_3, pm25_conv)
    df_4f <- cbind(df_4, pm25_conv)
    rm(pm25_conv)
    
    train_1 <- df_1f[b2_folds[[fi]][[1]], ]
    train_2 <- df_2f[b2_folds[[fi]][[1]], ]
    train_3 <- df_3f[b2_folds[[fi]][[1]], ]
    train_4 <- df_4f[b2_folds[[fi]][[1]], ]
    
    test_1 <- df_1f[b2_folds[[fi]][[2]], ]
    test_2 <- df_2f[b2_folds[[fi]][[2]], ]
    test_3 <- df_3f[b2_folds[[fi]][[2]], ]
    test_4 <- df_4f[b2_folds[[fi]][[2]], ]
    
    set.seed(fi*1)
    fit_1a <- ranger(pm25_value ~ ., data = train_1, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*2)
    fit_2a <- ranger(pm25_value ~ ., data = train_2, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*3)
    fit_3a <- ranger(pm25_value ~ ., data = train_3, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    set.seed(fi*4)
    fit_4a <- ranger(pm25_value ~ ., data = train_4, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    
    ## Get predictions for each parameter combination of interest
    df_save_b2[b2_folds[[fi]][[2]], paste0("p1a_", i)] <- 
      predict(fit_1a, test_1)$predictions
    df_save_b2[b2_folds[[fi]][[2]], paste0("p2a_", i)] <- 
      predict(fit_2a, test_2)$predictions
    df_save_b2[b2_folds[[fi]][[2]], paste0("p3a_", i)] <- 
      predict(fit_3a, test_3)$predictions
    df_save_b2[b2_folds[[fi]][[2]], paste0("p4a_", i)] <- 
      predict(fit_4a, test_4)$predictions
    
    rm(fit_1a, fit_2a, fit_3a, fit_4a)
    gc()
    
    # set.seed(fi*1)
    # fit_1c <- ranger(pm25_value ~ ., data = train_1, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*2)
    # fit_2c <- ranger(pm25_value ~ ., data = train_2, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*3)
    # fit_3c <- ranger(pm25_value ~ ., data = train_3, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # set.seed(fi*4)
    # fit_4c <- ranger(pm25_value ~ ., data = train_4, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # 
    # df_save_b2[b2_folds[[fi]][[2]], paste0("p1c_", i)] <- 
    #   predict(fit_1c, test_1)$predictions
    # df_save_b2[b2_folds[[fi]][[2]], paste0("p2c_", i)] <- 
    #   predict(fit_2c, test_2)$predictions
    # df_save_b2[b2_folds[[fi]][[2]], paste0("p3c_", i)] <- 
    #   predict(fit_3c, test_3)$predictions
    # df_save_b2[b2_folds[[fi]][[2]], paste0("p4c_", i)] <- 
    #   predict(fit_4c, test_4)$predictions
    # rm(fit_1c, fit_2c, fit_3c, fit_4c)
    # gc()
    
    rm(train_1, train_2, train_3, train_4,
       test_1, test_2, test_3, test_4,
       df_1f, df_2f, df_3f, df_4f)
    gc()
  }
}
rm(i, fi)
saveRDS(df_save_b2, "./pm25/cvresults/rftemp_b2.rds")

comb_save <- list(df_save_cv = df_save_cv,
                  df_save_b1 = df_save_b1,
                  df_save_b2 = df_save_b2,
                  param = param)
saveRDS(comb_save, "./pm25/cvresults/rftemp_comb.rds")
