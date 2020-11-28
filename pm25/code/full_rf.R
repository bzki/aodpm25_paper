# Fit overall (i.e. spatiotemporal) ranger model to entire data 
#  Generate full predictions

# A note on including time/day in random forest models
# - Just et al 2018 include days since 1970
# - Hu et al. 2017 include dummy day variables. 
# - Staffogia et al. 2019 include day of the year (integer) and day of the week

# My current approach: 
# -- day (as integer), 
# -- day of week (as integer) -- can also include indicator for day of week. 

# I run the following models, all including convolution PM2.5
# (a)- 4 standard without convolution 
# (b)- 4 standard with convolution [MAIN MODEL]
# (c)- 4 with convolution, and forcing (x,y)
# Try just 4 parameters combinations: mtry = 4, 8, 12, 16 for nodesize = 5. 
# NOTE
# -- When you force a variable in the split, it is 
#  *in addition* to mtry other variables.

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
  # FIXME: Remove useless variables here
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
# 1: No AOD, No GC
df_1b <- df_stack2 %>%
  select(-cmaq_id,  -aod_value, 
         -aod_imputed, -aod_gc_combine, -gc_aod)
# 2: GC
df_2b <- df_stack2 %>%
  select(-cmaq_id,  -aod_value, 
         -aod_imputed, -aod_gc_combine)
# 3: AOD Imputed
df_3b <- df_stack2 %>%
  select(-cmaq_id,  -aod_value, 
         -aod_gc_combine)
# 4: AOD + GC combined variable
df_4b <- df_stack2 %>%
  select(-cmaq_id,  -aod_value, 
         -aod_imputed, -gc_aod)

# (a) versiosn lack convolution PM2.5
df_1a <- df_1b %>% select(-conv_pm25)
df_2a <- df_2b %>% select(-conv_pm25)
df_3a <- df_3b %>% select(-conv_pm25)
df_4a <- df_4b %>% select(-conv_pm25)

train_1a <- df_1a %>% filter(!is.na(pm25_value))
train_2a <- df_2a %>% filter(!is.na(pm25_value))
train_3a <- df_3a %>% filter(!is.na(pm25_value))
train_4a <- df_4a %>% filter(!is.na(pm25_value))

train_1b <- df_1b %>% filter(!is.na(pm25_value))
train_2b <- df_2b %>% filter(!is.na(pm25_value))
train_3b <- df_3b %>% filter(!is.na(pm25_value))
train_4b <- df_4b %>% filter(!is.na(pm25_value))

import_list_a <- vector("list", 4)
import_list_b <- vector("list", 4)
import_list_c <- vector("list", 4)

# Define a prediction matrix
pred_names <- c(
  paste0("p", seq(4), "a_1"),
  paste0("p", seq(4), "a_2"),
  paste0("p", seq(4), "a_3"),
  paste0("p", seq(4), "a_4"),
  paste0("p", seq(4), "b_1"),
  paste0("p", seq(4), "b_2"),
  paste0("p", seq(4), "b_3"),
  paste0("p", seq(4), "b_4"),
  paste0("p", seq(4), "c_1"),
  paste0("p", seq(4), "c_2"),
  paste0("p", seq(4), "c_3"),
  paste0("p", seq(4), "c_4")
)
pred_matrix <- matrix(NA, nrow = nrow(df_1a), ncol = length(pred_names))
class(pred_matrix) <- "numeric"
colnames(pred_matrix) <- pred_names

source("./functions/rpred_chunked.R", echo = T, max.deparse.length = Inf)
# Get predictions in a manner that does not crash everything.
# -- do it in smaller chunks, ideally.

for (i in seq(nrow(param))) {
  print(i)

  set.seed(1)
  print("Model 1a")
  fit_1a <- ranger(pm25_value ~ ., data = train_1a, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_a[[i]][[1]] <- fit_1a$variable.importance
  pred_matrix[, paste0("p1a_", i)] <- 
    rpred_chunked(fit_1a, df_1a, chunk_size = 50000, verbose = F)
  rm(fit_1a)
  gc()
  
  set.seed(2)
  print("Model 2a")
  fit_2a <- ranger(pm25_value ~ ., data = train_2a, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_a[[i]][[2]] <- fit_2a$variable.importance
  pred_matrix[, paste0("p2a_", i)] <- 
    rpred_chunked(fit_2a, df_2a, chunk_size = 50000, verbose = F)
  rm(fit_2a)
  gc()
  
  set.seed(3)
  print("Model 3a")
  fit_3a <- ranger(pm25_value ~ ., data = train_3a, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_a[[i]][[3]] <- fit_3a$variable.importance
  pred_matrix[, paste0("p3a_", i)] <- 
    rpred_chunked(fit_3a, df_3a, chunk_size = 50000, verbose = F)
  rm(fit_3a)
  gc()
  
  set.seed(4)
  print("Model 4a")
  fit_4a <- ranger(pm25_value ~ ., data = train_4a, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_a[[i]][[4]] <- fit_4a$variable.importance
  pred_matrix[, paste0("p4a_", i)] <- 
    rpred_chunked(fit_4a, df_4a, chunk_size = 50000, verbose = F)
  rm(fit_4a)
  gc()

  set.seed(1)
  print("Model 1b")
  fit_1b <- ranger(pm25_value ~ ., data = train_1b, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_b[[i]][[1]] <- fit_1b$variable.importance
  pred_matrix[, paste0("p1b_", i)] <- 
    rpred_chunked(fit_1b, df_1b, chunk_size = 50000, verbose = F)
  rm(fit_1b)
  gc()
  
  set.seed(2)
  print("Model 2b")
  fit_2b <- ranger(pm25_value ~ ., data = train_2b, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_b[[i]][[2]] <- fit_2b$variable.importance
  pred_matrix[, paste0("p2b_", i)] <- 
    rpred_chunked(fit_2b, df_2b, chunk_size = 50000, verbose = F)
  rm(fit_2b)
  gc()

  set.seed(3)
  print("Model 3b")
  fit_3b <- ranger(pm25_value ~ ., data = train_3b, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_b[[i]][[3]] <- fit_3b$variable.importance
  pred_matrix[, paste0("p3b_", i)] <- 
    rpred_chunked(fit_3b, df_3b, chunk_size = 50000, verbose = F)
  rm(fit_3b)
  gc()
  
  set.seed(4)
  print("Model 4b")
  fit_4b <- ranger(pm25_value ~ ., data = train_4b, importance = "permutation",
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_b[[i]][[4]] <- fit_4b$variable.importance
  pred_matrix[, paste0("p4b_", i)] <- 
    rpred_chunked(fit_4b, df_4b, chunk_size = 50000, verbose = F)
  rm(fit_4b)
  gc()
  
  set.seed(1)
  print("Model 1c")
  fit_1c <- ranger(pm25_value ~ ., data = train_1b, importance = "permutation",
                   always.split.variables	= c("cmaq_x", "cmaq_y"),
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_c[[i]][[1]] <- fit_1c$variable.importance
  pred_matrix[, paste0("p1c_", i)] <- 
    rpred_chunked(fit_1c, df_1b, chunk_size = 50000, verbose = F)
  rm(fit_1c)
  gc()
  
  set.seed(2)
  print("Model 2c")
  fit_2c <- ranger(pm25_value ~ ., data = train_2b, importance = "permutation",
                   always.split.variables	= c("cmaq_x", "cmaq_y"),
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_c[[i]][[2]] <- fit_2c$variable.importance
  pred_matrix[, paste0("p2c_", i)] <- 
    rpred_chunked(fit_2c, df_2b, chunk_size = 50000, verbose = F)
  rm(fit_2c)
  gc()
  
  set.seed(3)
  print("Model 3c")
  fit_3c <- ranger(pm25_value ~ ., data = train_3b, importance = "permutation",
                   always.split.variables	= c("cmaq_x", "cmaq_y"),
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_c[[i]][[3]] <- fit_3c$variable.importance
  pred_matrix[, paste0("p3c_", i)] <- 
    rpred_chunked(fit_3c, df_3b, chunk_size = 50000, verbose = F)
  rm(fit_3c)
  gc()
  
  set.seed(4)
  print("Model 4c")
  fit_4c <- ranger(pm25_value ~ ., data = train_4b, importance = "permutation",
                   always.split.variables	= c("cmaq_x", "cmaq_y"),
                   num.trees = NUM.TREES, mtry = param$mtry[i], 
                   min.node.size = param$min.node.size[i])
  import_list_c[[i]][[4]] <- fit_4c$variable.importance
  pred_matrix[, paste0("p4c_", i)] <- 
    rpred_chunked(fit_4c, df_4b, chunk_size = 50000, verbose = F)
  rm(fit_4c)
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
saveRDS(all_list, "./pm25/rf_fullpred.rds")

