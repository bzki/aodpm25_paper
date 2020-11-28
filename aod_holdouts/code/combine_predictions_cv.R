# Combine cross-validated data.frames from LatticeKrig and RandomForest
# Save in a data.frame list
# Examine cross-validated MSE in several ways

library(testthat) # Some tests
library(dplyr) # For easier manipulation and joining of data.frames
library(sessioninfo)
sessionInfo()
sessioninfo::session_info()

# Functions to calculate R2, MSE, RMSE, Int/Slope
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
  return(as.numeric(coef( lm( true ~ pred) )))
}

#### Create data.frame list of combined lk/rf results ####

days = 182:212
df_list <- vector("list", length(days)) 
names(df_list) <- paste0("day", days)

for (day in days) {
  print(day)
  
  lk_df <- readRDS(paste0("./aod_holdouts/cvdata/lk_cv_4_", day, ".rds")) %>%
    arrange(cmaq_id)
  rf_df <- readRDS(paste0("./aod_holdouts/cvdata/rf_cv_4_", day, ".rds")) %>%
    arrange(cmaq_id)
  
  testthat::expect_identical(lk_df$cmaq_id, rf_df$cmaq_id)
  testthat::expect_identical(lk_df$foldID, rf_df$foldID)
  testthat::expect_identical(lk_df$aod_value, rf_df$aod_value)
  
  comb_df <- left_join(lk_df, rf_df %>% dplyr::select(-aod_value, -foldID),
                       by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(day = day)
  
  day_name <- paste0("day", day)
  df_list[[day_name]] <- comb_df
  
}
rm(comb_df, lk_df, rf_df, day, day_name)
# Save
saveRDS(df_list, "./aod_holdouts/cvdata/cv_list.rds")

stacked_df <- do.call("rbind", df_list)

#### Assess MSE in several ways ####
# 1: For each fold separately, then averaged.
# 2: For each fold separately, then weighted averaged based on size of fold
# 3: Overall -- should be equivalent to (2)

cv_results <- stacked_df %>%
  group_by(day, foldID) %>%
  summarize(
    fold_size = n(),
    lk_mse = getMSE(true = aod_value, pred = lk_cv),
    rf_mse = getMSE(true = aod_value, pred = rf_cv)
  ) %>%
  ungroup(.) %>%
  group_by(day) %>%
  summarize(
    lk_mse_ave = mean(lk_mse),
    rf_mse_ave = mean(rf_mse),
    lk_mse_wtd = weighted.mean(lk_mse, fold_size),
    rf_mse_wtd = weighted.mean(rf_mse, fold_size)
  ) %>%
  ungroup(.)

cv_overall <- stacked_df %>%
  group_by(day) %>%
  summarize(
    lk_mse_wtd = getMSE(aod_value, lk_cv),
    rf_mse_wtd = getMSE(aod_value, rf_cv)
  ) %>%
  ungroup(.)

testthat::expect_equal(cv_overall$lk_mse_wtd, cv_results$lk_mse_wtd)
testthat::expect_equal(cv_overall$rf_mse_wtd, cv_results$rf_mse_wtd)

length(which(cv_results$lk_mse_ave < cv_results$rf_mse_ave))
length(which(cv_results$rf_mse_ave < cv_results$lk_mse_ave))

length(which(cv_results$lk_mse_wtd < cv_results$rf_mse_wtd))
length(which(cv_results$rf_mse_wtd < cv_results$lk_mse_wtd))

# LatticeKrig has a 2 to 1 advantage, roughly. Still -- variability in results

# Construct MSE based on quantiles
dist_q <-
  cut(stacked_df$cv_nndist_km,
      breaks = quantile(stacked_df$cv_nndist_km, 
                        seq(0, 1, 0.1)),
      include.lowest = T)

y_cv_split <- split(stacked_df$aod_value, dist_q)
lk_cv_split <- split(stacked_df$lk_cv, dist_q)
rf_cv_split <- split(stacked_df$rf_cv, dist_q)

quantile_df <- data.frame(dist_label = names(y_cv_split),
                          LK_MSE = NA, RF_MSE = NA)

for (sp_i in seq(length(y_cv_split))) {
  
  quantile_df[sp_i, "LK_MSE"] <- 
    getMSE(true = y_cv_split[[sp_i]],
           pred = lk_cv_split[[sp_i]])

  quantile_df[sp_i, "RF_MSE"] <- 
    getMSE(true = y_cv_split[[sp_i]],
           pred = rf_cv_split[[sp_i]])
  
}

plot(quantile_df$LK_MSE) 
points(quantile_df$RF_MSE, col = "blue", pch = 2)


#### Loop through and estimate daily weights ####
tmp_df <- df_list$day190
coef(lm(aod_value ~ -1 + lk_cv + rf_cv, data = tmp_df))
coef(lm(aod_value ~ -1 + lk_cv + rf_cv, data = stacked_df))

# - daily super learner
library(nnls)
coef_estimates <- matrix(NA, nrow = length(days), ncol = 3)
colnames(coef_estimates) <- c("day", "lk", "rf")
coef_estimates <- data.frame(coef_estimates)
nnls_estimates <- coef_estimates

for (day in days) {
  day_name <- paste0("day", day)
  day_df <- df_list[[day_name]]
  
  day_lm <- lm(aod_value ~ -1 + lk_cv + rf_cv, data = day_df)
  day_index <- which(days == day)
  coef_estimates[day_index, "day"] <- day
  coef_estimates[day_index, c("lk", "rf")] <- as.numeric(coef(day_lm))
  
  day_nnls <- nnls(A = as.matrix(day_df[, c("lk_cv", "rf_cv")]),
                   b = day_df$aod_value)
  
  nnls_estimates[day_index, "day"] <- day
  nnls_estimates[day_index, c("lk", "rf")] <- coef(day_nnls)
}

## Normalize and then adjust estimates so that max weight is 0.9 for method
nnls_normalize <- nnls_estimates %>%
  mutate(
    lk_norm = (lk / (lk + rf)),
    rf_norm = (rf / (lk + rf)),
    lk_norm_adj = pmin(lk_norm, 0.9),
    lk_norm_adj = pmax(lk_norm_adj, 0.1),
    rf_norm_adj = pmin(rf_norm, 0.9),
    rf_norm_adj = pmax(rf_norm_adj, 0.1)
  )

daily_sl <- nnls_normalize
# saveRDS(dail_sl, "./aod_holdouts/cvdata/sl_daily.rds")


# Function for doing normalized weights
do_normalize <- function(nnls_est, min.val = 0.1, max.val = 0.9) {
  nnls_norm <- nnls_est / sum(nnls_est)
  nnls_norm2 <- pmin(nnls_norm, 0.9)
  nnls_norm3 <- pmax(nnls_norm2, 0.1)
  return(nnls_norm3)
}

# - daily distance-based super learner. 
# -- how to do this?
# FIXME: Come back to this when you have a better sense of things. 
# for (day in days) {
#   print(day)
#   day_name <- paste0("day", day)
#   day_df <- df_list[[day_name]]
#   day_index <- which(days == day)
# 
#   # Put data into 5 to 10 evenly-spaced distances
#   dist_q <-
#     cut(day_df$cv_nndist_km,
#         c(seq(0, 200, 40), Inf),
#         include.lowest = T)
# 
#   y_cv_split <- split(day_df$aod_value, dist_q)
#   lk_cv_split <- split(day_df$lk_cv, dist_q)
#   rf_cv_split <- split(day_df$rf_cv, dist_q)
# 
#   # Estimate weights for each grouping
#   dweight_est <- NULL
#   for (di in seq(length(y_cv_split))) {
#     testthat::expect_equal(names(lk_cv_split), names(rf_cv_split))
#     di_nnls <- nnls(A = cbind(lk_cv_split[[di]], rf_cv_split[[di]]),
#                     b = y_cv_split[[di]])
#     di_nnls_norm <- do_normalize(coef(di_nnls))
#     testthat::expect_equal(sum(di_nnls_norm), 1)
#     dweight_est <- rbind(dweight_est, di_nnls_norm)
#   }
# 
#   dist_mid <- (seq(0, 225, by = 25) + seq(25, 250, by = 25)) / 2
#   plot(loess(dweight_est[, 1] ~ dist_mid))
# 
#   # Fit a smooth curve
# 
#   # Save the curve itself
# 
# }

#### - overall super learner ####
stacked_nnls <- nnls(A = as.matrix(stacked_df[, c("lk_cv", "rf_cv")]),
                     b = stacked_df$aod_value)

overall_sl <- do_normalize(coef(stacked_nnls))
names(overall_sl) <- c("lk", "rf")

# - overall distance-based super learner.
dist_q <-
  cut(stacked_df$cv_nndist_km,
      c(seq(0, 300, 25), Inf),
      include.lowest = T)

y_cv_split <- split(stacked_df$aod_value, dist_q)
lk_cv_split <- split(stacked_df$lk_cv, dist_q)
rf_cv_split <- split(stacked_df$rf_cv, dist_q)

# Estimate weights for each grouping
dweight_est <- NULL
for (di in seq(length(y_cv_split))) {
  testthat::expect_equal(names(lk_cv_split), names(rf_cv_split))
  di_nnls <- nnls(A = cbind(lk_cv_split[[di]], rf_cv_split[[di]]),
                  b = y_cv_split[[di]])
  di_nnls_norm <- do_normalize(coef(di_nnls))
  testthat::expect_equal(sum(di_nnls_norm), 1)
  dweight_est <- rbind(dweight_est, di_nnls_norm)
}

dist_mid <- (seq(0, 300, by = 25) + seq(25, 325, by = 25)) / 2
plot(dist_mid, dweight_est[, 1])
points(loess.smooth(dist_mid, dweight_est[, 1]), type = "l")
overall_dist <- data.frame(cbind(dist_mid, dweight_est))
colnames(overall_dist) <- c("dist", "lk", "rf")
rownames(overall_dist) <- NULL

weight_results <- list(
  overall_sl = overall_sl,
  daily_sl = daily_sl,
  overall_dist = overall_dist
)
saveRDS(weight_results, "./aod_holdouts/cvdata/sl_weights.rds")
