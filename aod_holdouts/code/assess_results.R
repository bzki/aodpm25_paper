# This file will assess test results using various methods
# LK alone, RF alone, and various methods of combining.

# library(testthat) # Some tests
library(fields)
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

# Function for doing normalized weights
do_normalize <- function(nnls_est, min.val = 0.1, max.val = 0.9) {
  nnls_norm <- nnls_est / sum(nnls_est)
  nnls_norm2 <- pmin(nnls_norm, 0.9)
  nnls_norm3 <- pmax(nnls_norm2, 0.1)
  return(nnls_norm3)
}

#### Read in super learner weights ####
# -- See combine_predictions_cv.R for the creation of these weights.
sl_weights <- readRDS("./aod_holdouts/cvdata/sl_weights.rds")
overall_sl <- sl_weights$overall_sl
daily_sl <- sl_weights$daily_sl
overall_dist <- sl_weights$overall_dist
# Fit loess with defaults to overall_dist
loess_fit <- loess(lk ~ dist, data = overall_dist, span = 2/3, degree = 1,
                   control = loess.control(surface = "direct"))
plot(overall_dist$dist, overall_dist$lk)
points(overall_dist$dist, predict(loess_fit, overall_dist), type = "l") 

#### Loop through data ####
days = 182:212
df_list <- vector("list", length(days))
names(df_list) <- paste0("day", days)
for (day in days) {
  print(day)
  # Read in predicted data and join
  lk_df <- readRDS(paste0("./aod_holdouts/testdata/lk_test_4_", day, ".rds"))
  rf_df <- readRDS(paste0("./aod_holdouts/testdata/rf_test_4_", day, ".rds"))
  testthat::expect_identical(lk_df$cmaq_id, rf_df$cmaq_id)
  testthat::expect_identical(lk_df$aod_value, rf_df$aod_value)
  comb_df <- left_join(lk_df, rf_df %>% dplyr::select(-aod_value),
                       by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(day = day)
  test_xy <- as.matrix(comb_df[, c("cmaq_x", "cmaq_y")])
  
  # Read in full training data
  train_xy <- as.matrix(
    readRDS(paste0("./rawrds/split_4_", day, ".rds")) %>%
    filter(split != "Test") %>%
    dplyr::select(cmaq_x, cmaq_y)
  )
  
  # Calculate distance between test and training data
  nn_meters <- 
    apply(rdist(test_xy, train_xy), 1, min)
  comb_df$nn_m <- nn_meters
  comb_df$nn_km <- nn_meters/1000
  
  day_weight <- daily_sl[which(daily_sl$day == day), ]
  
  dist_df <- data.frame(dist = comb_df$nn_km)
  dist_df$weight1 <- predict(loess_fit, dist_df)
  # plot(dist_df$dist, dist_df$weight1)
  dist_df$weight2 <- pmax(dist_df$weight1, 0.1)
  dist_df$weight3 <- pmin(dist_df$weight2, 0.9)
  
  comb_df$dist_weight <- dist_df$weight3
  rm(dist_df)
  
  # Create various combinations of interest using the weights I calculated
  comb_df <- comb_df %>%
    mutate(
      ave_comb = (lk_predict + rf_predict)/2,
      sl_1 = overall_sl[1]*lk_predict + overall_sl[2]*rf_predict,
      sl_2 = day_weight$lk_norm_adj*lk_predict + day_weight$rf_norm_adj*rf_predict,
      sl_3 = dist_weight*lk_predict + (1 - dist_weight)*rf_predict
    )
  day_name <- paste0("day", day)
  df_list[[day_name]] <- comb_df
  rm(day_weight, comb_df, nn_meters, test_xy, train_xy, lk_df, rf_df)
}

#### Summarize results by day and overall ####
res_names <- c("day", 
               "LK_R2", "RF_R2", "Ave_R2", "SL1_R2", "SL2_R2", "SL3_R2",
               "LK_RMSE", "RF_RMSE", "Ave_RMSE", "SL1_RMSE", "SL2_RMSE", "SL3_RMSE",
               "LK_MSE", "RF_MSE", "Ave_MSE", "SL1_MSE", "SL2_MSE", "SL3_MSE",
               "LK_Int", "LK_Slope", 
               "RF_Int", "RF_Slope",
               "Ave_Int", "Ave_Slope",
               "SL1_Int", "SL1_Slope",
               "SL2_Int", "SL2_Slope",
               "SL3_Int", "SL3_Slope")

res_df <- matrix(NA, nrow = length(days), ncol = length(res_names))
colnames(res_df) <- res_names
res_df <- as.data.frame(res_df)

q_names <- c("day", 
             paste0("Q_", seq(10)),
             paste0("LK_MSE_", seq(10)),
             paste0("RF_MSE_", seq(10)))
q_df <- matrix(NA, nrow = length(days), ncol = length(q_names))
colnames(q_df) <- q_names
q_df <- as.data.frame(q_df)

for (day in days) {
  print(day)
  day_index <- which(days == day)
  day_name <- paste0("day", day)
  
  day_df <- df_list[[day_name]]
  res_df$LK_R2[day_index] <- getR2(day_df$aod_value, day_df$lk_predict)
  res_df$RF_R2[day_index] <- getR2(day_df$aod_value, day_df$rf_predict)
  res_df$Ave_R2[day_index] <- getR2(day_df$aod_value, day_df$ave_comb)
  res_df$SL1_R2[day_index] <- getR2(day_df$aod_value, day_df$sl_1)
  res_df$SL2_R2[day_index] <- getR2(day_df$aod_value, day_df$sl_2)
  res_df$SL3_R2[day_index] <- getR2(day_df$aod_value, day_df$sl_3)
  
  res_df$LK_RMSE[day_index] <- getRMSE(day_df$aod_value, day_df$lk_predict)
  res_df$RF_RMSE[day_index] <- getRMSE(day_df$aod_value, day_df$rf_predict)
  res_df$Ave_RMSE[day_index] <- getRMSE(day_df$aod_value, day_df$ave_comb)
  res_df$SL1_RMSE[day_index] <- getRMSE(day_df$aod_value, day_df$sl_1)
  res_df$SL2_RMSE[day_index] <- getRMSE(day_df$aod_value, day_df$sl_2)
  res_df$SL3_RMSE[day_index] <- getRMSE(day_df$aod_value, day_df$sl_3)
  
  res_df$LK_MSE[day_index] <- getMSE(day_df$aod_value, day_df$lk_predict)
  res_df$RF_MSE[day_index] <- getMSE(day_df$aod_value, day_df$rf_predict)
  res_df$Ave_MSE[day_index] <- getMSE(day_df$aod_value, day_df$ave_comb)
  res_df$SL1_MSE[day_index] <- getMSE(day_df$aod_value, day_df$sl_1)
  res_df$SL2_MSE[day_index] <- getMSE(day_df$aod_value, day_df$sl_2)
  res_df$SL3_MSE[day_index] <- getMSE(day_df$aod_value, day_df$sl_3)
  
  res_df[day_index, c("LK_Int", "LK_Slope")] <- 
    getIntSlope(day_df$aod_value, day_df$lk_predict)
  res_df[day_index, c("RF_Int", "RF_Slope")] <- 
    getIntSlope(day_df$aod_value, day_df$rf_predict)
  res_df[day_index, c("Ave_Int", "Ave_Slope")] <- 
    getIntSlope(day_df$aod_value, day_df$ave_comb)
  res_df[day_index, c("SL1_Int", "SL1_Slope")] <- 
    getIntSlope(day_df$aod_value, day_df$sl_1)
  res_df[day_index, c("SL2_Int", "SL2_Slope")] <- 
    getIntSlope(day_df$aod_value, day_df$sl_2)
  res_df[day_index, c("SL3_Int", "SL3_Slope")] <- 
    getIntSlope(day_df$aod_value, day_df$sl_3)
  
  res_df[day_index, "day"] <- day
  
  
  ## Daily quantile
  xy_test_q10 <-
    cut(day_df$nn_km,
        breaks = quantile(day_df$nn_km, 
                          seq(0, 1, by = 0.1)),
        include.lowest = T)
  
  # Split the AOD values by the groupings
  y_test_split <- split(day_df$aod_value, xy_test_q10)
  lk_test_split <- split(day_df$lk_predict, xy_test_q10)
  rf_test_split <- split(day_df$rf_predict, xy_test_q10)
  
  q_df[day_index, paste0("Q_", seq(10))] <- names(y_test_split)
  
  for (sp_i in seq(length(y_test_split))) {
    
    q_df[day_index, paste0("LK_MSE_", sp_i)] <- 
      getMSE(true = y_test_split[[sp_i]],
             pred = lk_test_split[[sp_i]])
    q_df[day_index, paste0("RF_MSE_", sp_i)] <- 
      getMSE(true = y_test_split[[sp_i]],
             pred = rf_test_split[[sp_i]])
  
  }
  
}

length(which(res_df$LK_R2 > res_df$RF_R2))
length(which(res_df$LK_MSE < res_df$RF_MSE))

colMeans(res_df) * 100

## Final -- summarize all data combined
stacked_df <- do.call("rbind", df_list)
lk_rmse <- getRMSE(stacked_df$aod_value, stacked_df$lk_predict)
rf_rmse <- getRMSE(stacked_df$aod_value, stacked_df$rf_predict)
ave_rmse <- getRMSE(stacked_df$aod_value, stacked_df$ave_comb)
sl1_rmse <- getRMSE(stacked_df$aod_value, stacked_df$sl_1)
sl2_rmse <- getRMSE(stacked_df$aod_value, stacked_df$sl_2)
sl3_rmse <- getRMSE(stacked_df$aod_value, stacked_df$sl_3)

lk_mse <- getMSE(stacked_df$aod_value, stacked_df$lk_predict)
rf_mse <- getMSE(stacked_df$aod_value, stacked_df$rf_predict)
ave_mse <- getMSE(stacked_df$aod_value, stacked_df$ave_comb)
sl1_mse <- getMSE(stacked_df$aod_value, stacked_df$sl_1)
sl2_mse <- getMSE(stacked_df$aod_value, stacked_df$sl_2)
sl3_mse <- getMSE(stacked_df$aod_value, stacked_df$sl_3)

lk_r2 <- getR2(stacked_df$aod_value, stacked_df$lk_predict)
rf_r2 <- getR2(stacked_df$aod_value, stacked_df$rf_predict)
ave_r2 <- getR2(stacked_df$aod_value, stacked_df$ave_comb)
sl1_r2 <- getR2(stacked_df$aod_value, stacked_df$sl_1)
sl2_r2 <- getR2(stacked_df$aod_value, stacked_df$sl_2)
sl3_r2 <- getR2(stacked_df$aod_value, stacked_df$sl_3)

lk_is <- getIntSlope(stacked_df$aod_value, stacked_df$lk_predict)
rf_is <- getIntSlope(stacked_df$aod_value, stacked_df$rf_predict)
ave_is <- getIntSlope(stacked_df$aod_value, stacked_df$ave_comb)
sl1_is <- getIntSlope(stacked_df$aod_value, stacked_df$sl_1)
sl2_is <- getIntSlope(stacked_df$aod_value, stacked_df$sl_2)
sl3_is <- getIntSlope(stacked_df$aod_value, stacked_df$sl_3)


c(lk_rmse, rf_rmse, ave_rmse, sl1_rmse, sl2_rmse, sl3_rmse) * 100
c(lk_r2, rf_r2, ave_r2, sl1_r2, sl2_r2, sl3_r2) * 100
round(c(lk_mse, rf_mse, ave_mse, sl1_mse, sl2_mse, sl3_mse) * 100, 3)

(sl1_mse/lk_mse - 1) * 100
(sl3_mse/lk_mse - 1) * 100

(sl1_rmse/lk_rmse - 1) * 100 # -2.30 %
(sl3_rmse/lk_rmse - 1) * 100 # -2.34 %

#### Make overall table ####
library(xtable)
final_names <- c("Method", "R2", "MSE (x100)", 
                 "Intercept", "Slope")
final_df <- matrix(NA, nrow = 6, ncol = length(final_names))
colnames(final_df) <- final_names
final_df <- as.data.frame(final_df)
final_df$Method <- 
  c("LatticeKrig", "Random Forest", "LK-RF Average",
    "SuperLearner: Overall Weights", "SuperLearner - Daily Weights",
    "SuperLearner: Distance-based weights")
final_df$R2 <-
  c(lk_r2, rf_r2, ave_r2, sl1_r2, sl2_r2, sl3_r2)
final_df[, "MSE (x100)"] <- 
  c(lk_mse, rf_mse, ave_mse, sl1_mse, sl2_mse, sl3_mse) * 100
final_df[, c("Intercept", "Slope")] <-
  rbind(lk_is, rf_is, ave_is, sl1_is, sl2_is, sl3_is)

knitr::kable(final_df, digits = 3)
print(xtable(final_df))
print(xtable(final_df, digits = 3), include.rownames = F,
      file = "./paper_plots/sl_plots/result_table.txt")

#### Alternate tables with and RMSE ####
altf_names <- c("Method", "R2", "RMSE", "Intercept", "Slope")
altf_df <- matrix(NA, nrow = 6, ncol = length(altf_names))
colnames(altf_df) <- altf_names
altf_df <- as.data.frame(altf_df)
altf_df$Method <- 
  c("LatticeKrig", "Random Forest", "LK-RF Average",
    "SuperLearner: Overall Weights", "SuperLearner - Daily Weights",
    "SuperLearner: Distance-based weights")
altf_df$R2 <-
  c(lk_r2, rf_r2, ave_r2, sl1_r2, sl2_r2, sl3_r2)
altf_df[, "RMSE"] <- 
  c(lk_rmse, rf_rmse, ave_rmse, sl1_rmse, sl2_rmse, sl3_rmse)*100
altf_df[, c("Intercept", "Slope")] <-
  rbind(lk_is, rf_is, ave_is, sl1_is, sl2_is, sl3_is)

knitr::kable(altf_df, digits = c(0, 3, 2, 2, 2))
print(xtable(altf_df, digits = c(0, 0, 3, 2, 2, 2)))
print(xtable(altf_df, digits = c(0, 0, 3, 2, 2, 2)), include.rownames = F,
      file = "./paper_plots/sl_plots/result_table_alt.txt")

#### Figures of MSE/R2/RMSE  ####
mse_range <- range(res_df %>% dplyr::select(contains("_MSE")))
r2_range <- range(res_df %>% dplyr::select(contains("_R2")))

res_df$day_1 <- res_df$day - 181
library(RColorBrewer)
display.brewer.pal(n = 6, name = "Set1")
brew_colors <- RColorBrewer::brewer.pal(n = 6, name = "Set1")
mse_set <- paste0(c("RF", "LK", "Ave", "SL1", "SL2", "SL3"), "_MSE")
r2_set <- paste0(c("RF", "LK", "Ave", "SL1", "SL2", "SL3"), "_R2")
pch_set <- c("R", "L", "A", "1", "2", "3")
pch_set2 <- c(0, 1, 3, 7, 8, 9)

# par(mfrow = c(2, 1))
png("./paper_plots/sl_plots/r2mse_4.png",
    height = 900, width = 700)
par(mar = c(4, 4.5, 2.5, 0.5))
m_layout <- matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE)
layout(mat = m_layout, heights = c(0.45,0.45,0.1))

for (mi in seq(mse_set)) {
  if (mi == 1) {
    plot(res_df$day_1, res_df[, mse_set[mi]] * 100, ylim = mse_range * 100,
         pch = pch_set2[mi], col = brew_colors[mi],
         ylab = "MSE (x100)", xlab = "", cex.lab = 1.75, cex = 1.5, cex.axis = 1.75,
         cex.main = 1.75,
         main = "Daily MSE")
  }
  else {
    points(res_df$day_1, res_df[, mse_set[mi]] * 100, cex = 1.5, 
           pch = pch_set2[mi], col = brew_colors[mi])
  }
}
for (mi in seq(r2_set)) {
  if (mi == 1) {
    plot(res_df$day_1, res_df[, r2_set[mi]] * 100, ylim = r2_range * 100,
         pch = pch_set2[mi], col = brew_colors[mi],
         ylab = "R2 (x100)", xlab = "July 2011 Day", cex.lab = 1.75, cex = 1.5, cex.axis = 1.75,
         cex.main = 1.75,
         main = "Daily R2")
  }
  else {
    points(res_df$day_1, res_df[, r2_set[mi]] * 100, cex = 1.5,
           pch = pch_set2[mi], col = brew_colors[mi])
  }
}
par(mar = c(0, 0, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top", inset = 0,
       legend = c("RF", "LK", "Ave", "SL1", "SL2", "SL3"), 
       col=brew_colors, pch = pch_set2, 
       horiz = TRUE, cex = 1.75)
dev.off()

#### Re-do in ggplot2 ####
library(ggplot2)
library(tidyr)
library(forcats)
res_df2 <- res_df %>%
  mutate(MSE_RF = RF_MSE * 100,
         MSE_LK = LK_MSE * 100,
         MSE_Ave = Ave_MSE * 100,
         MSE_SL1 = SL1_MSE * 100,
         MSE_SL2 = SL2_MSE * 100,
         MSE_SL3 = SL3_MSE * 100,
         R2_RF = RF_R2,
         R2_LK = LK_R2,
         R2_Ave = Ave_R2,
         R2_SL1 = SL1_R2,
         R2_SL2 = SL2_R2,
         R2_SL3 = SL3_R2,
         RMSE_RF = RF_RMSE,
         RMSE_LK = LK_RMSE,
         RMSE_Ave = Ave_RMSE,
         RMSE_SL1 = SL1_RMSE,
         RMSE_SL2 = SL2_RMSE,
         RMSE_SL3 = SL3_RMSE,
         RMSE_RF_x100 = RF_RMSE*100,
         RMSE_LK_x100 = LK_RMSE*100,
         RMSE_Ave_x100 = Ave_RMSE*100,
         RMSE_SL1_x100 = SL1_RMSE*100,
         RMSE_SL2_x100 = SL2_RMSE*100,
         RMSE_SL3_x100 = SL3_RMSE*100,
)

best_df <- res_df2 %>%
  select(MSE_RF:MSE_SL3)
# Add rank for each method, and show which is best on a particular day
best_names <- c("RF", "LK", "Ave", "SL1", "SL2", "SL3")
rank_days <- as.data.frame(
  t(apply(best_df %>% select(MSE_RF:MSE_SL3),
                     1, function(x) { rank(x) }))
)
colnames(rank_days) <- paste0(best_names, "_rank")
best_df$best_col <- apply(best_df %>% select(MSE_RF:MSE_SL3), 1, function(x) {
                            best_names[which.min(x)]
                          })
best_df$worst_col <- apply(best_df %>% select(MSE_RF:MSE_SL3), 1, function(x) {
  best_names[which.max(x)]
})


table(best_df$best_col)
rank_sum <- 
  rank_days %>% select(LK_rank, RF_rank, Ave_rank:SL3_rank) %>%
  summarize(
    across(.fns = list(median = median, mean = mean))
  )

rank_mean <- rank_sum %>% select(ends_with("_mean"))
rank_median <- rank_sum %>% select(ends_with("_median"))
testthat::expect_equivalent(rank_median,
                            as.integer(rank_median))
best_days = table(best_df$best_col)[c("LK", "RF", "Ave", "SL1", "SL2", "SL3")]
worst_days = table(best_df$worst_col)[c("LK", "RF", "Ave", "SL1", "SL2", "SL3")]

# Data.frame for averages across days rather than pooling all test data together
# Add to the final_df output
final_df2 <- cbind(final_df, 
                   Median_rank = as.integer(as.vector(t(rank_median))),
                   Mean_rank = as.vector(t(rank_mean)),
                   Best_days = as.integer(as.vector(best_days)),
                   Worst_days = as.integer(as.vector(worst_days))) %>%
  mutate(Worst_days = replace(Worst_days, is.na(Worst_days), 0L))


knitr::kable(final_df2, digits = 3)
print(xtable(final_df2,  digits = c(0, 0, 3, 3, 2, 2, 0, 2, 0, 0)))
print(xtable(final_df2, digits = c(0, 0, 3, 3, 2, 2, 0, 2, 0, 0)), include.rownames = F,
      file = "./paper_plots/sl_plots/result_table2.txt")

final_df3 <- final_df2 %>% select(-Median_rank, -Mean_rank)
knitr::kable(final_df3, digits = 3)
print(xtable(final_df3,  digits = c(0, 0, 3, 3, 2, 2, 0, 0)))
print(xtable(final_df3, digits = c(0, 0, 3, 3, 2, 2, 0, 0)), include.rownames = F,
      file = "./paper_plots/sl_plots/result_table3.txt")

#### Construct alternative version with RMSE x100 ####
altf_df2 <- cbind(altf_df, 
                   Median_rank = as.integer(as.vector(t(rank_median))),
                   Mean_rank = as.vector(t(rank_mean)),
                   Best_days = as.integer(as.vector(best_days)),
                   Worst_days = as.integer(as.vector(worst_days))) %>%
  mutate(Worst_days = replace(Worst_days, is.na(Worst_days), 0L))


knitr::kable(altf_df2, digits = 3)
print(xtable(altf_df2,  digits = c(0, 0,  3, 2, 2, 2, 0, 2, 0, 0)))
print(xtable(altf_df2, digits = c(0, 0,  3, 2, 2, 2, 0, 2, 0, 0)), include.rownames = F,
      file = "./paper_plots/sl_plots/result_table_alt2.txt")

altf_df3 <- altf_df2 %>% select(-Median_rank, -Mean_rank)
knitr::kable(altf_df3, digits = c(0, 3, 2, 2, 2))
print(xtable(altf_df3,  digits = c(0, 0,  3, 2, 2, 2, 0, 0)))
print(xtable(altf_df3, digits = c(0, 0,  3, 2, 2, 2, 0, 0)), include.rownames = F,
      file = "./paper_plots/sl_plots/result_table_alt3.txt")

# Construct a "long" data.frame
long_mse <- 
  pivot_longer(res_df2 %>% dplyr::select(day_1, MSE_RF:MSE_SL3), 
               cols = MSE_RF:MSE_SL3,
               names_to = "Method",
               names_prefix = "MSE_") %>%
  mutate(Method = as_factor(Method))

long_rmse <- 
  pivot_longer(res_df2 %>% dplyr::select(day_1, RMSE_RF:RMSE_SL3) %>%
                 mutate(RMSE_RF = RMSE_RF*100,
                        RMSE_LK = RMSE_LK*100,
                        RMSE_Ave = RMSE_Ave*100,
                        RMSE_SL1 = RMSE_SL1*100,
                        RMSE_SL2 = RMSE_SL2*100,
                        RMSE_SL3 = RMSE_SL3*100), 
               cols = RMSE_RF:RMSE_SL3,
               names_to = "Method",
               names_prefix = "RMSE_") %>%
  mutate(Method = as_factor(Method))

long_r2 <- 
  pivot_longer(res_df2 %>% dplyr::select(day_1, R2_RF:R2_SL3), 
               cols = R2_RF:R2_SL3,
               names_to = "Method",
               names_prefix = "R2_") %>%
  mutate(Method = as_factor(Method))

mse_plot <- 
  ggplot() + 
  geom_point(data = long_mse, 
             aes(x = day_1, y = value, color = Method, shape = Method), 
             alpha = 0.5, size = 2.5) +
  scale_color_brewer(palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  xlab("Day (July 2011)") +
  ylab("MSE (x100)") +
  ggtitle("MSE") 

r2_plot <- 
ggplot() + 
  geom_point(data = long_r2, 
             aes(x = day_1, y = value, color = Method, shape = Method), 
             alpha = 0.5, size = 2.5) +
  scale_color_brewer(palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  xlab("Day (July 2011)") +
  ylab("R2") +
  ggtitle("R2")

rmse_plot <- 
  ggplot() + 
  geom_point(data = long_rmse, 
             aes(x = day_1, y = value, color = Method, shape = Method), 
             alpha = 0.5, size = 2.5) +
  scale_color_brewer(palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  xlab("Day (July 2011)") +
  ylab("RMSE (x100)") +
  ggtitle("RMSE") 


library(ggpubr)
combined_mser2 <- 
  ggarrange(mse_plot, r2_plot, nrow = 2, ncol = 1,
          common.legend = T, legend = "bottom")
ggsave("./paper_plots/sl_plots/combined_mser2.png", combined_mser2)
combined_rmser2 <- 
  ggarrange(rmse_plot, r2_plot, nrow = 2, ncol = 1,
            common.legend = T, legend = "bottom")
ggsave("./paper_plots/sl_plots/combined_rmser2.png", combined_rmser2,
       width = 8.5, height = 8.5)

# ggplot2 is not correctly ascertaining my colors if color is inside aes()
# p = ggplot() + 
#   geom_point(data = res_df2, aes(x = day_1, y = RF_MSE), col = brew_colors[1], alpha = 0.5) +
#   geom_point(data = res_df2, aes(x = day_1, y = LK_MSE), col = brew_colors[2], alpha = 0.5) + 
#   geom_point(data = res_df2, aes(x = day_1, y = Ave_MSE), col = brew_colors[3], alpha = 0.5) + 
#   geom_point(data = res_df2, aes(x = day_1, y = SL1_MSE), col = brew_colors[4], alpha = 0.5) + 
#   geom_point(data = res_df2, aes(x = day_1, y = SL2_MSE), col = brew_colors[5], alpha = 0.5) + 
#   geom_point(data = res_df2, aes(x = day_1, y = SL3_MSE), col = brew_colors[6], alpha = 0.5) + 
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
#   xlab("Day (July 2011)") +
#   ylab("MSE (x100)") +
#   ggtitle("MSE for each day")

#### Quantile Plot ####
# Summarize by quantile
xy_test_q10 <-
  cut(stacked_df$nn_km,
      breaks = quantile(stacked_df$nn_km, seq(0, 1, by = 0.1)),
      include.lowest = T)

# Split the AOD values by the groupings
y_test_split <- split(stacked_df$aod_value, xy_test_q10)
lk_test_split <- split(stacked_df$lk_predict, xy_test_q10)
rf_test_split <- split(stacked_df$rf_predict, xy_test_q10)

quant_df <- data.frame(q_label = as_factor(names(y_test_split)))
quant_df$LK_MSE <- NA
quant_df$RF_MSE <- NA
quant_df$LK_RMSE <- NA
quant_df$RF_RMSE <- NA
quant_df$LK_R2 <- NA
quant_df$RF_R2 <- NA

for (sp_i in seq(length(y_test_split))) {
  
  quant_df[sp_i, "LK_MSE"] <- 
    getMSE(true = y_test_split[[sp_i]],
           pred = lk_test_split[[sp_i]])
  quant_df[sp_i, "RF_MSE"] <- 
    getMSE(true = y_test_split[[sp_i]],
           pred = rf_test_split[[sp_i]])

  quant_df[sp_i, "LK_RMSE"] <- 
    getRMSE(true = y_test_split[[sp_i]],
           pred = lk_test_split[[sp_i]])
  quant_df[sp_i, "RF_RMSE"] <- 
    getRMSE(true = y_test_split[[sp_i]],
           pred = rf_test_split[[sp_i]])
  
  quant_df[sp_i, "LK_R2"] <- 
    getR2(true = y_test_split[[sp_i]],
           pred = lk_test_split[[sp_i]])
  quant_df[sp_i, "RF_R2"] <- 
    getR2(true = y_test_split[[sp_i]],
           pred = rf_test_split[[sp_i]])
  
}

quant_range <- range(quant_df[, -1]) * 100
png("./paper_plots/sl_plots/quant_summary.png",
    width = 800, height = 600)
par(mar = c(4.2, 4.2, 4, 2))
plot(seq(10), quant_df$LK_MSE*100, ylim = quant_range,
     xlab = "Distance grouping (km)", xaxt = "n", ylab = "MSE (x100)",
     main = "MSE by distance to nearest training point",
     pch = pch_set2[2], col = brew_colors[2], cex = 1.75,
     cex.lab = 1.5, cex.main = 1.5)
points(seq(10), quant_df$RF_MSE*100,
       pch = pch_set2[1], col = brew_colors[1], cex = 1.75)
axis(side = 1, at = seq(10), labels = quant_df$q_label)
legend(x = "topleft", 
       legend = c("RF", "LK"), 
       col=brew_colors[1:2], pch = pch_set2[1:2], 
       horiz = TRUE, cex = 1.75)
dev.off()

lkqmean <- colMeans(q_df[, paste0("LK_MSE_", seq(10))])
rfqmean <- colMeans(q_df[, paste0("RF_MSE_", seq(10))])
plot(lkqmean)
points(rfqmean, pch = 2)

#### quantile redo ggplot2 ####
# -- use Set1 colorbrewer
library(tidyr)
library(forcats)
long_quantmse <- 
  pivot_longer(quant_df %>% 
                 mutate(quant_id = seq(nrow(.)),
                        MSE_LK = LK_MSE * 100,
                        MSE_RF = RF_MSE * 100) %>%
                 select(-LK_MSE, -RF_MSE), 
               cols = c(MSE_RF, MSE_LK),
               names_to = "Method",
               names_prefix = "MSE_") %>%
  mutate(Method = as_factor(Method))

long_quantrmse <- 
  pivot_longer(quant_df %>% 
                 mutate(quant_id = seq(nrow(.)),
                        RMSE_LK = LK_RMSE * 100,
                        RMSE_RF = RF_RMSE * 100) %>%
                 select(-LK_RMSE, -RF_RMSE), 
               cols = c(RMSE_RF, RMSE_LK),
               names_to = "Method",
               names_prefix = "RMSE_") %>%
  mutate(Method = as_factor(Method))

long_quantr2 <- 
  pivot_longer(quant_df %>% 
                 mutate(quant_id = seq(nrow(.)),
                        R2_LK = LK_R2,
                        R2_RF = RF_R2) %>%
                 select(-LK_R2, -RF_R2), 
               cols = c(R2_RF, R2_LK),
               names_to = "Method",
               names_prefix = "R2_") %>%
  mutate(Method = as_factor(Method))

qmse_plot <- 
  ggplot() + 
  geom_point(data = long_quantmse, 
             aes(x = q_label, y = value, color = Method, shape = Method), 
             alpha = 0.5, size = 2.5) +
  scale_color_brewer(palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  xlab("Distance quantile") +
  ylab("MSE (x100)") +
  ggtitle("MSE by distance to nearest observed point") 

qrmse_plot <- 
  ggplot() + 
  geom_point(data = long_quantrmse, 
             aes(x = q_label, y = value, color = Method, shape = Method), 
             alpha = 0.5, size = 4) +
  scale_color_brewer(palette = "Set1") + 
  xlab("Distance quantile (km)") +
  ylab("RMSE (x100)") +
  ggtitle("RMSE by distance to nearest observed point") +
  theme_grey(
    base_size = 15,
    base_family = ""
  ) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

qr2_plot <- 
  ggplot() + 
  geom_point(data = long_quantr2, 
             aes(x = q_label, y = value, color = Method, shape = Method), 
             alpha = 0.5, size = 2.5) +
  scale_color_brewer(palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  xlab("Distance quantile") +
  ylab("R2") +
  ggtitle("R2 by distance to nearest observed point") 

comb_quant <- ggarrange(qmse_plot, qr2_plot, nrow = 2, ncol = 1,
                        common.legend = T, legend = "bottom")
ggsave("./paper_plots/sl_plots/quant_mser2.png", comb_quant,
       width = 7, height = 7)

ggsave("./paper_plots/sl_plots/quant_msesum.png", qmse_plot)

ggsave("./paper_plots/sl_plots/quant_rmsesum.png", qrmse_plot,
       width = 12, height = 8)

