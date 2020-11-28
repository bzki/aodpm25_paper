# Addendum to full_rf.R

# Fit overall (i.e. spatiotemporal) ranger model to entire data 
#  Generate full predictions

# NOTE: This file contains an additional model to what's in full_rf.R
# The additional model trains only on areas where AOD is observed
# Prediction is done everywhere, using using imputed AOD where AOD not observed

# The additional models again fall into the three broad classes
# - (a) without convolution
# - (b) with convolution [MAIN MODEL]
# - (c) with convolution, and forcing (x,y)

# A note on including time/day in random forest models
# - Just et al 2018 include days since 1970
# - Hu et al. 2017 include dummy day variables. 
# - Staffogia et al. 2019 include day of the year (integer) and day of the week

# My current approach: 
# -- day (as integer), 
# -- day of week (as integer) -- can also include indicator for day of week. 

# Try just 4 parameters combinations: mtry = 4, 8, 12, 16 for nodesize = 5. 
# NOTE TO SELF:
# -- When you force a variable in the split, it is *in addition* to mtry
# See the ranger documentation for this explanation: help(ranger)

# Load libraries
library(dplyr)
library(ranger)
library(fields)
source("./functions/make_convolution.R", echo = T, max.deparse.length = Inf)
sessionInfo()

# Set parameter information
NUM.TREES = 2000 # We love trees, folks
mtry_options <- c(4, 8, 12, 16)
node_options <- c(5)
param <- expand.grid(mtry = mtry_options, min.node.size = node_options) %>%
  arrange(mtry, min.node.size) %>%
  mutate(pnum = seq(nrow(.)))
rm(mtry_options, node_options)
attributes(param)$out.attrs <- NULL

# Loop through *full* daily data, save into list
days = 182:212
raw_list <- vector("list", length(days))
names(raw_list) <- paste0("day", days)
for (day in days) {
    day_name <- paste0("day", day)
    raw_list[[day_name]] <- readRDS(paste0("./rawrds/n4_day", day, ".rds"))
    rm(day_name)
}
rm(day)

# Construct convolution PM for each day
df_list <- vector("list", length(days))
names(df_list) <- paste0("day", days)
for (day in days) {
    day_name <- paste0("day", day)
    print(day)
    df <- raw_list[[day_name]] %>%
      mutate(x_km = cmaq_x/1000, y_km = cmaq_y/1000)
    
    # Where is PM2.5 observed ("training") and no observed ("testing")?
    holdout_list <- list(vector("list", 2))
    holdout_list[[1]][[1]] <- which(!is.na(df$pm25_value))
    holdout_list[[1]][[2]] <- which(is.na(df$pm25_value))

    # Convolution PM2.5
    conv_pm25 <- 
      make_convolution(df_input = df, x_name = "x_km", y_name = "y_km",
                       var_name = "pm25_value", fold_list = holdout_list,
                       scale_num = 100000)

    # Save
    df$conv_pm25 <- as.vector(conv_pm25)
    df_list[[day_name]] <- df
    rm(day_name, conv_pm25, df, holdout_list)
}
rm(day)
rm(raw_list)

# Create a stacked version of the above -- with all days together
df_stack <- do.call("rbind", df_list)
rm(df_list)

# Read in full prediction AOD
pred_list <- readRDS("./aod_holdouts/finalpred/allpreds.rds")
pred_stack <- do.call("rbind", pred_list) %>%
  select(cmaq_id, day, sl3_pred)
rm(pred_list)

# Link together and save the full dataset
df_stack2 <- left_join(df_stack, pred_stack,
                       by = c("cmaq_id" = "cmaq_id", "day" = "day")) %>%
  # Remove useless variables here
  select(-year, -starts_with("narr"), -lat, -lon, -x_km, -y_km) %>%
  relocate(cmaq_id, day, cmaq_x, cmaq_y) %>%
  # Simpler name for the imputed AOD I am including
  rename(aod_imputed = sl3_pred) %>%
  # Create a missing indicator and an AOD/GC combined variable
  mutate(
    aod_missing = ifelse(is.na(aod_value), 1, 0),
    aod_gc_combine = ifelse(is.na(aod_value), gc_aod, aod_value)
  ) %>%
  # Create day of the week variable
  mutate(
    # Friday is 5; Saturday is 6; Sunday is 0; etc.
    # Day 182 is July 1, 2011, Friday
    day_of_week = (day - 177L) %% 7L
  )
rm(df_stack, pred_stack)

# Basic data.frame where we will save all predictions
df_save <- df_stack2 %>%
  select(cmaq_id, day, cmaq_x, cmaq_y, 
         pm25_value, aod_value, aod_imputed, gc_aod, conv_pm25)

## Create data.frames to be used for ranger
# 5: Train on areas where AOD observed; predict on combined AOD (obs + imp)
# (b) has convolution layer
df_5b <- df_stack2 %>%
  select(-cmaq_id, -aod_gc_combine, -aod_missing)

# (a) versions lack convolution PM2.5
df_5a <- df_5b %>% select(-conv_pm25)

# Train solely on areas where AOD is observed
train_5a <- df_5a %>% filter(!is.na(pm25_value)) %>%
      filter(!is.na(aod_value)) %>%
      select(-aod_imputed)
train_5b <- df_5b %>% filter(!is.na(pm25_value)) %>%
      filter(!is.na(aod_value)) %>%
      select(-aod_imputed)

# Use full data for predictions
# Replace missing AOD values with the imputed for the prediction model
pred_5a <- df_5a %>%
      mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
      select(-aod_imputed)
pred_5b <- df_5b %>%
      mutate(aod_value = ifelse(is.na(aod_value), aod_imputed, aod_value)) %>%
      select(-aod_imputed)

testthat::expect_identical(names(train_5a), names(pred_5a))
testthat::expect_identical(names(train_5b), names(pred_5b))

import_list_a <- vector("list", 4)
import_list_b <- vector("list", 4)
import_list_c <- vector("list", 4)

# Define a prediction matrix
pred_names <- c(
  paste0("p5a_", seq(nrow(param))), 
  paste0("p5b_", seq(nrow(param))),
  paste0("p5c_", seq(nrow(param)))
)
pred_matrix <- matrix(NA, nrow = nrow(df_5a), ncol = length(pred_names))
class(pred_matrix) <- "numeric"
colnames(pred_matrix) <- pred_names

source("./functions/rpred_chunked.R", echo = T, max.deparse.length = Inf)
# Predictions are not very memory-efficient in ranger.
# In order to get predictions in a manner that does not crash everything,
# I use a short helper function that performs predictions in chunks. 

for (i in seq(nrow(param))) {
  print(i)

  set.seed(5)
  print("Model 5a")
  fit_5a <- ranger(pm25_value ~ ., data = train_5a, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_a[[i]][[1]] <- fit_5a$variable.importance
  pred_matrix[, paste0("p5a_", i)] <- 
    rpred_chunked(fit_5a, pred_5a, chunk_size = 50000, verbose = F)
  rm(fit_5a)
  gc()
  
  set.seed(5)
  print("Model 5b")
  fit_5b <- ranger(pm25_value ~ ., data = train_5b, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_b[[i]][[1]] <- fit_5b$variable.importance
  pred_matrix[, paste0("p5b_", i)] <- 
    rpred_chunked(fit_5b, pred_5b, chunk_size = 50000, verbose = F)
  rm(fit_5b)
  gc()
  
  set.seed(5)  
  print("Model 5c")
  # Note: same training and prediction datasets as (b) -- doesnt change
  fit_5c <- ranger(pm25_value ~ ., data = train_5b, importance = "permutation",
                   always.split.variables	= c("cmaq_x", "cmaq_y"),
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_c[[i]][[1]] <- fit_5c$variable.importance
  pred_matrix[, paste0("p5c_", i)] <- 
    rpred_chunked(fit_5c, pred_5b, chunk_size = 50000, verbose = F)
  rm(fit_5c)
  gc()
  
}

df_save <- cbind(df_save, pred_matrix)

# Save everything
all_list <- list(
  df_save = df_save,
  import_list_a = import_list_a,
  import_list_b = import_list_b,
  import_list_c = import_list_c,
  param = param
)
saveRDS(all_list, "./pm25/rf_fullpred_2.rds")

