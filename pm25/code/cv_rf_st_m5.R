# ST models -- additional model (M5) -- see cv_rf_st.R for further details
# -- train solely on areas where AOD is actually observed
# -- predict everywhere by combining observed + imputed AOD

# - Be careful about calculating convolution here
# - A single set of features considered only
# - Restrict to nodesize = 5 and 4 mtry values
# - convolution, and force x/y split (just those)
# - Name these 5a/5c, consistent with before. 

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
attributes(param)$out.attrs <- NULL
rm(mtry_options, node_options)
NUM.TREES = 2000

# Create a dataset that will have all of the predictions and be saved
df_save_cv <- df_stack2 %>% select(day, cmaq_id, cmaq_x, cmaq_y, pm25_value, 
                                   aod_value, aod_missing, aod_imputed, fold_cv)
df_save_b1 <- df_stack2 %>% select(day, cmaq_id, cmaq_x, cmaq_y, pm25_value, 
                                   aod_value, aod_missing, aod_imputed, fold_b1)
df_save_b2 <- df_stack2 %>% select(day, cmaq_id, cmaq_x, cmaq_y, pm25_value, 
                                   aod_value, aod_missing, aod_imputed, fold_b2)

# Create template data.frame
df_5 <- df_stack2 %>% select(-cmaq_id, -fold_cv, -fold_b1, -fold_b2, 
                              -aod_gc_combine, -aod_missing)

# set.seed(51)
# Loop through each CV, each mtry, each fold, and each model
for (i in seq(nrow(param))) {
  print(paste("param number", i))
  # a => random forest
  # c => rf + keep (x,y) at every split
  df_save_cv[, paste0("p5a_", i)] <- NA
  df_save_cv[, paste0("p5c_", i)] <- NA

  for (fi in seq(length(cv_folds))) {
    
    print(paste("fold number", fi)) 
    
    # Assign the appropriate convolution PM2.5
    pm25_conv <- conv_cv_stack[, paste0("pm25_f", fi)]
    df_5f <- cbind(df_5, pm25_conv)
    rm(pm25_conv)
    
    train_5 <- df_5f[cv_folds[[fi]][[1]], ] %>%
      filter(!is.na(aod_value)) %>%
      select(-aod_imputed)
    
    test_5 <- df_5f[cv_folds[[fi]][[2]], ] %>%
      mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
      select(-aod_imputed)
    
    testthat::expect_identical(names(train_5), names(test_5))

    set.seed(fi*5)
    fit_5a <- ranger(pm25_value ~ ., data = train_5, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    
    ## Get predictions for each parameter combination of interest
    df_save_cv[cv_folds[[fi]][[2]], paste0("p5a_", i)] <- 
      predict(fit_5a, test_5)$predictions

    rm(fit_5a)
    gc()
    
    # set.seed(fi*5)
    # fit_5c <- ranger(pm25_value ~ ., data = train_5, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # 
    # df_save_cv[cv_folds[[fi]][[2]], paste0("p5c_", i)] <- 
    #   predict(fit_5c, test_5)$predictions
    # rm(fit_5c)
    # gc()

    rm(train_5, test_5, df_5f)
    gc()
  }
}
rm(i, fi)
saveRDS(df_save_cv, "./pm25/cvresults/rf_st_cv_results2.rds")


set.seed(52)
for (i in seq(nrow(param))) {
  print(paste("param number", i))
  # a => random forest
  # c => rf + keep (x,y) at every split
  df_save_b1[, paste0("p5a_", i)] <- NA
  df_save_b1[, paste0("p5c_", i)] <- NA
  
  for (fi in seq(length(b1_folds))) {
    
    print(paste("fold number", fi)) 
    
    # Assign the appropriate convolution PM2.5
    pm25_conv <- conv_b1_stack[, paste0("pm25_f", fi)]
    df_5f <- cbind(df_5, pm25_conv)
    rm(pm25_conv)
    
    train_5 <- df_5f[b1_folds[[fi]][[1]], ] %>%
      filter(!is.na(aod_value)) %>%
      select(-aod_imputed)
    
    test_5 <- df_5f[b1_folds[[fi]][[2]], ] %>%
      mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
      select(-aod_imputed)
    
    testthat::expect_identical(names(train_5), names(test_5))
    
    set.seed(fi*5)
    fit_5a <- ranger(pm25_value ~ ., data = train_5, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    
    ## Get predictions for each parameter combination of interest
    df_save_b1[b1_folds[[fi]][[2]], paste0("p5a_", i)] <- 
      predict(fit_5a, test_5)$predictions
    
    rm(fit_5a)
    gc()
    
    # set.seed(fi*5)
    # fit_5c <- ranger(pm25_value ~ ., data = train_5, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # 
    # df_save_b1[b1_folds[[fi]][[2]], paste0("p5c_", i)] <- 
    #   predict(fit_5c, test_5)$predictions
    # rm(fit_5c)
    # gc()
    
    rm(train_5, test_5, df_5f)
    gc()
  }
}
rm(i, fi)

saveRDS(df_save_b1, "./pm25/cvresults/rf_st_b1_results2.rds")


set.seed(53)
for (i in seq(nrow(param))) {
  print(paste("param number", i))
  # a => random forest
  # c => rf + keep (x,y) at every split
  df_save_b2[, paste0("p5a_", i)] <- NA
  df_save_b2[, paste0("p5c_", i)] <- NA
  
  for (fi in seq(length(b2_folds))) {
    
    print(paste("fold number", fi)) 
    
    # Assign the appropriate convolution PM2.5
    pm25_conv <- conv_b2_stack[, paste0("pm25_f", fi)]
    df_5f <- cbind(df_5, pm25_conv)
    rm(pm25_conv)
    
    train_5 <- df_5f[b2_folds[[fi]][[1]], ] %>%
      filter(!is.na(aod_value)) %>%
      select(-aod_imputed)
    
    test_5 <- df_5f[b2_folds[[fi]][[2]], ] %>%
      mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
      select(-aod_imputed)
    
    testthat::expect_identical(names(train_5), names(test_5))
    
    set.seed(fi*5)
    fit_5a <- ranger(pm25_value ~ ., data = train_5, 
                     num.trees = NUM.TREES, mtry = param$mtry[i], 
                     min.node.size = param$min.node.size[i])
    
    ## Get predictions for each parameter combination of interest
    df_save_b2[b2_folds[[fi]][[2]], paste0("p5a_", i)] <- 
      predict(fit_5a, test_5)$predictions
    
    rm(fit_5a)
    gc()
    
    # set.seed(fi*5)
    # fit_5c <- ranger(pm25_value ~ ., data = train_5, 
    #                  always.split.variables	= c("cmaq_x", "cmaq_y"),
    #                  num.trees = NUM.TREES, mtry = param$mtry[i], 
    #                  min.node.size = param$min.node.size[i])
    # 
    # df_save_b2[b2_folds[[fi]][[2]], paste0("p5c_", i)] <- 
    #   predict(fit_5c, test_5)$predictions
    # rm(fit_5c)
    # gc()
    
    rm(train_5, test_5, df_5f)
    gc()
  }
}
rm(i, fi)

saveRDS(df_save_b2, "./pm25/cvresults/rf_st_b2_results2.rds")

comb_save <- list(df_save_cv = df_save_cv,
                  df_save_b1 = df_save_b1,
                  df_save_b2 = df_save_b2,
                  param = param,
              descriptions = c("a = standard convolution layer random forest",
                               "c = include (x,y) split at every split"))
saveRDS(comb_save, "./pm25/cvresults/rf_st_results2.rds")
rm(list=ls())
