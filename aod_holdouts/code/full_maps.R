# Make full maps
# Includes the following predictors
# LatticeKrig (now adding in elevation in fixed part of model)
# RandomForest
# 4 Average/SuperLearner methods
# -- See "combine_predictions_cv_full.R" for superLearner weight construction

# NOTE (potential future)
# Methods of improving final predictions for the combined methods, and 
#   perhaps in the cross-validation stage as well.
# -> If LatticeKrig is outside bounds of training data, set pmin/pmax to 
#   observed training data. Or -- use random forest directly in these cases.
# -> Possibly look for edge issues with LatticeKrig and 
#   switch to Random Forest in such cases
# -> Use these kinds of procedures in the final maps.


library(dplyr) # For easier manipulation and joining of data.frames
library(fields) # Fast distance calculations
library(sessioninfo)
sessionInfo()
sessioninfo::session_info()

# Function for doing normalized weights
# do_normalize <- function(nnls_est, min.val = 0.1, max.val = 0.9) {
#   nnls_norm <- nnls_est / sum(nnls_est)
#   nnls_norm2 <- pmin(nnls_norm, max.val)
#   nnls_norm3 <- pmax(nnls_norm2, min.val)
#   return(nnls_norm3)
# }

# Chunked nearest neighbor distance calculation
nn_chunked <- function(x1, x2, chunk_size = 5000) {
  start_row <- 1
  end_row <- min(chunk_size, nrow(x1))
  nn_dist <- NULL
  while (end_row <= nrow(x1)) {
    print(paste0("Distance for ", start_row, " to ", end_row))
    nn_dist <- 
      c(nn_dist, 
        apply(
          rdist(x1[start_row:end_row, , drop = F], x2), 
        1, min)
      )
    if (end_row == nrow(x1)) {
      break
    } else {
      start_row <- start_row + chunk_size
      end_row <- min(end_row + chunk_size, nrow(x1))
    }
  }
  return(nn_dist)
}

# test_x1 <- cbind(seq(100), seq(100))
# test_x2 <- cbind(seq(10), seq(10))
# test_x3 <- cbind(rnorm(50), rnorm(50))
# testthat::expect_identical(
#   nn_chunked(x1 = test_x1, x2 = test_x2, chunk_size = 3),
#   apply(rdist(test_x1, test_x2), 1, min))
# testthat::expect_identical(
#   nn_chunked(x1 = test_x1, x2 = test_x2, chunk_size = 6),
#   apply(rdist(test_x1, test_x2), 1, min))
# testthat::expect_identical(
#   nn_chunked(x1 = test_x1, x2 = test_x2, chunk_size = 1),
#   apply(rdist(test_x1, test_x2), 1, min))
# testthat::expect_identical(
#   nn_chunked(x1 = test_x1, x2 = test_x2, chunk_size = 2),
#   apply(rdist(test_x1, test_x2), 1, min))
# testthat::expect_identical(
#   nn_chunked(x1 = test_x1, x2 = test_x2, chunk_size = 5),
#   apply(rdist(test_x1, test_x2), 1, min))
# testthat::expect_identical(
#   nn_chunked(x1 = test_x2, x2 = test_x3, chunk_size = 3),
#   apply(rdist(test_x2, test_x3), 1, min))
# rm(test_x1, test_x2, test_x3)

#### Read in super learner weights ####
# -- See combined_predictions_cv_full.R for creation of these. 
sl_weights <- readRDS("./aod_holdouts/cvdata_full/slfull_weights.rds")
overall_sl <- sl_weights$overall_sl
daily_sl <- sl_weights$daily_sl
overall_dist <- sl_weights$overall_dist
# Fit loess with defaults to overall_dist
loess_fit <- loess(lk ~ dist, data = overall_dist, span = 2/3, degree = 1,
                   control = loess.control(surface = "direct"))
pdf(file = "./aod_holdouts/paper_plots/sl_plots/slfull_weights.pdf",
    width = 7, height = 5)
plot(overall_dist$dist, overall_dist$lk, xlab = "Distance (km)",
     ylab = "LK Weight",
     cex = 1.5)
points(overall_dist$dist, predict(loess_fit, overall_dist), type = "l") 
dev.off()
#### Construct daily data.frames ####
days = 182:212
df_list <- vector("list", length(days))
names(df_list) <- paste0("day", days)
lk_full <- readRDS("./aod_holdouts/finalpred/lk_final_elev.rds")$fpredict_list
rf_full <- readRDS("./aod_holdouts/finalpred/rf_final.rds")$fpredict_list

# Contruct combined estimators
# -- Will focus on using AOD where observed, predictions where not
# But we also provide full predictions.

for (day in days) {

  print(day)
  day_name <- paste0("day", day)
  
  # Get distance from non-observed to observed
  pred_df <- lk_full[[day_name]] %>% 
    dplyr::filter(holdout == 1) %>%
    dplyr::select(cmaq_id, x_km, y_km)
  
  pred_xy <- as.matrix(pred_df[, c("x_km", "y_km")])
  obs_xy <- as.matrix(lk_full[[day_name]] %>%
    dplyr::filter(holdout == 0) %>%
    dplyr::select(x_km, y_km))
  
  # Need to do a chunked version of distance calculation to prevent crashing
  # -- break into 5 to 10k chunks
  nn_dist <- nn_chunked(x1=pred_xy, x2=obs_xy, chunk_size=5000)
  testthat::expect_equal(length(nn_dist), nrow(pred_xy))
  pred_df$nn_dist_km <- nn_dist
  rm(nn_dist)
  
  day_weight <- daily_sl[which(daily_sl$day == day), ]
  dist_df <- data.frame(dist = pred_df$nn_dist_km)
  dist_df$weight1 <- predict(loess_fit, dist_df)
  # plot(dist_df$dist, dist_df$weight1)
  dist_df$weight2 <- pmax(dist_df$weight1, 0.1)
  dist_df$weight3 <- pmin(dist_df$weight2, 0.9)
  pred_df$dist_weight <- dist_df$weight3
  rm(dist_df)
  
  # Approach -- focus on test dataset, combine predictions, then join
  #  into full data.frames
  # Read in predicted data and join
  lk_df <- lk_full[[day_name]] %>%
    #dplyr::filter(holdout == 1) %>%
    arrange(cmaq_id)
  rf_df <- rf_full[[day_name]] %>%
    #dplyr::filter(holdout == 1) %>%
    arrange(cmaq_id) %>%
    rename(rf_predict = ranger_predict)
  
  testthat::expect_identical(lk_df$cmaq_id, rf_df$cmaq_id)
  testthat::expect_identical(lk_df$aod_value, rf_df$aod_value)
  testthat::expect_identical(lk_df$holdout, rf_df$holdout)
  
  comb_df <- left_join(lk_df, rf_df %>% dplyr::select(-aod_value, -holdout),
                       by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(day = day) %>%
    left_join(pred_df %>% dplyr::select(-x_km, -y_km),
              by = c("cmaq_id" = "cmaq_id"))
  
  ## Plot code
  # plot(obs_xy, cex = 0.2)
  # points(pred_xy, col = "blue", cex = 0.2)
  # dist_points <- which(comb_df$nn_dist_km > quantile(comb_df$nn_dist_km, 0.99))
  # max_pt <- which.max(comb_df$nn_dist_km)
  # points(comb_df[max_pt, c("x_km", "y_km")], cex = 1, col = "red", pch = 16)
  # points(comb_df[dist_points, c("x_km", "y_km")], cex = 0.5, col = "purple")
  # tmp_dist <- as.vector(rdist(as.matrix(comb_df[max_pt, c("x_km", "y_km")]), obs_xy))
  # points(obs_xy[which.min(tmp_dist), , drop = F], cex = 1, col = "orange", pch = 16)
  
  # Create various combinations of interest using the weights I calculated
  full_df <- comb_df %>%
    mutate(
      ave_fpred = (lk_predict + rf_predict)/2,
      sl1_fpred = overall_sl[1]*lk_predict + overall_sl[2]*rf_predict,
      sl2_fpred = day_weight$lk_norm_adj*lk_predict + day_weight$rf_norm_adj*rf_predict,
      sl3_fpred = dist_weight*lk_predict + (1 - dist_weight)*rf_predict
    ) %>%
    rename(lk_fpred = lk_predict,
           rf_fpred = rf_predict) %>%
    mutate(
      rf_pred = if_else(!is.na(aod_value), aod_value, rf_fpred),
      lk_pred = if_else(!is.na(aod_value), aod_value, lk_fpred),
      ave_pred = if_else(!is.na(aod_value), aod_value, ave_fpred),
      sl1_pred = if_else(!is.na(aod_value), aod_value, sl1_fpred),
      sl2_pred = if_else(!is.na(aod_value), aod_value, sl2_fpred),
      sl3_pred = if_else(!is.na(aod_value), aod_value, sl3_fpred),
      lkrf_diff = lk_pred - rf_pred,
      lkrf_fdiff = lk_fpred - rf_fpred
    )

  df_list[[day_name]] <- full_df
  rm(day_weight, comb_df, full_df, lk_df, rf_df, obs_xy, pred_xy, pred_df)
}

saveRDS(df_list, "./aod_holdouts/finalpred/allpreds.rds")

#### Create maps ####
# -- Aiming to avoid using spplot 
df_list <- readRDS("./aod_holdouts/finalpred/allpreds.rds")
library(sf)
library(sp)
library(ggplot2)
library(ggpubr)

#### 1: Average plots daily ####
xy_matrix <- as.matrix(df_list$day182[, c("cmaq_x", "cmaq_y")])
full_id <- df_list$day182$cmaq_id
aod_matrix = matrix(NA, nrow = nrow(xy_matrix), ncol = length(days))
class(aod_matrix) <- "numeric"
colnames(aod_matrix) <- paste0("day", days)
lk_matrix = rf_matrix = ave_matrix = sl1_matrix = sl2_matrix = sl3_matrix = 
  lkrf_matrix = aod_matrix

for (day in days) {
  print(day)
  day_name <- paste0("day", day)
  # Check that order is same for all datasets
  testthat::expect_identical(full_id, df_list[[day_name]]$cmaq_id)
  # Record daily observed and predicted values into respective matrices
  aod_matrix[, day_name] <- df_list[[day_name]]$aod_value
  lk_matrix[, day_name] <- df_list[[day_name]]$lk_pred
  rf_matrix[, day_name] <- df_list[[day_name]]$rf_pred
  ave_matrix[, day_name] <- df_list[[day_name]]$ave_pred
  sl1_matrix[, day_name] <- df_list[[day_name]]$sl1_pred
  sl2_matrix[, day_name] <- df_list[[day_name]]$sl2_pred
  sl3_matrix[, day_name] <- df_list[[day_name]]$sl3_pred
  lkrf_matrix[, day_name] <- df_list[[day_name]]$lkrf_diff
}

aod_ave <- rowMeans(aod_matrix, na.rm = T)
lk_ave <- rowMeans(lk_matrix)
rf_ave <- rowMeans(rf_matrix)
ave_ave <- rowMeans(ave_matrix)
sl1_ave <- rowMeans(sl1_matrix)
sl2_ave <- rowMeans(sl2_matrix)
sl3_ave <- rowMeans(sl3_matrix)

overall_df <- data.frame(xy_matrix, cmaq_id = full_id,
                         aod_ave, lk_ave, rf_ave, ave_ave,
                         sl1_ave, sl2_ave, sl3_ave) %>%
  mutate(lk_rf_diff = lk_ave - rf_ave,
         lk_aod_diff = lk_ave - aod_ave,
         rf_aod_diff = rf_ave - aod_ave,
         ave_aod_diff = ave_ave - aod_ave,
         sl1_aod_diff = sl1_ave - aod_ave,
         sl2_aod_diff = sl2_ave - aod_ave,
         sl3_aod_diff = sl3_ave - aod_ave)
  
all_range <- range(overall_df %>% 
                     dplyr::select(aod_ave:sl3_ave),
                   na.rm = T)

p_aod <- ggplot(data=overall_df) +
  geom_point(aes(x = cmaq_x, y = cmaq_y, color = aod_ave),
             shape = 15, size = 0.6) + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        legend.text = element_text(size = 12)) + 
  ggtitle(paste0("Observed")) +
  scale_color_viridis_c(limits = all_range, option = "plasma") +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 20)) +
  coord_fixed()

make_plot1 <- function(df, var_name, main, val_range,
                       title_size = 18) {
  ggp <- ggplot(data=df) +
    geom_point(aes(x = cmaq_x, y = cmaq_y, color = var_name),
               shape = 15, size = 0.6) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = title_size, face = "bold"),
          panel.background = element_blank()) + 
    ggtitle(main) +
    scale_color_viridis_c(limits = val_range, option = "plasma") + 
    coord_fixed(ratio = 1) # Important for spatial plots
  return(ggp)
}

make_plot2 <- function(df, var_name, main, val_range,
                       title_size = 18) {
  ggp <- ggplot(data=df) +
    geom_point(aes(x = cmaq_x, y = cmaq_y, color = var_name),
               shape = 15, size = 0.6) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = title_size, face = "bold"),
          panel.background = element_blank()) + 
    ggtitle(main) +
    scale_color_viridis_c(limits = val_range, option = "plasma") + 
    guides(fill = guide_colourbar(barwidth = 2, barheight = 20)) +
    coord_fixed(ratio = 1) # Important for spatial plots
  return(ggp)
}

p_aod <- make_plot1(overall_df, aod_ave, main = "Observed", val_range = all_range)
p_aod2 <- make_plot2(overall_df, aod_ave, main = "Observed", val_range = all_range)

p_lk <- make_plot1(overall_df, lk_ave, main = "LatticeKrig", val_range = all_range)
p_rf <- make_plot1(overall_df, rf_ave, main = "Random Forest", val_range = all_range)
p_ave <- make_plot1(overall_df, ave_ave, main = "Average", val_range = all_range)
p_sl1 <- make_plot1(overall_df, sl1_ave, 
                    main = "Super Learner: Overall weights", val_range = all_range)
p_sl2 <- make_plot1(overall_df, sl2_ave, 
                    main = "Super Learner: Daily weights", val_range = all_range)
p_sl3 <- make_plot1(overall_df, sl3_ave, 
                    main = "Super Learner: Distance-based weights", val_range = all_range)

comb_plot <-
  ggarrange(p_aod, p_lk, p_rf, p_ave, p_sl1, p_sl2, p_sl3,
            ncol=2, nrow=4, common.legend = TRUE, legend="bottom")

comb_plot

sub_plot <- 
  ggarrange(p_aod, p_lk, p_rf, 
            font.label = list(size = 14),
            ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
sub_plot

ggsave(file = "./aod_holdouts/paper_plots/sl_plots/average.png",
       sub_plot, width = 8, height = 11.5)

#### Modify: Add state borders, larger legend ####
# Read in previously created sf object in USLCC with cmaq_id correctly formatted
x <- readRDS("cmaq/region_assignment.rds")
str(x)
st_crs(x)
# Attach data to the observed points
overall_sf <- left_join(x, overall_df, by = c("cmaq_id" = "cmaq_id"))

# Read in a U.S states shapefile and convert to USLCC
states_sf <- sf::read_sf("./cmaq/shapefile_tl2010_us_state_2010/US_state_2010.shp",
                         stringsAsFactors = F) %>%
  filter(!(NAME10 %in% c("Alaska", "Hawaii", "Puerto Rico")))
states_lcc_sf <- st_transform(states_sf,
                              crs = st_crs(x))

make_sfplot <- function(main_df, state_df, barwidth = 1, barheight = 30,
                        bar_direction = "vertical",
                        var_name, main, val_range, title_size = 18,
                        point_size = 0.4, border_size = 0.3) {
  sf_plot <- ggplot() +
    geom_sf(data = main_df, 
            aes_string(color = var_name), 
            shape = 15, size = point_size) + 
    geom_sf(data = st_geometry(state_df), color = "black", 
            size = border_size, alpha = 0) +
    # To plot projected rather than lat/lon
    coord_sf(datum = NULL) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = title_size, face = "bold"),
          panel.background = element_blank()) + 
    ggtitle(main) +
    scale_color_viridis_c(limits = val_range, option = "plasma") +
    guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight,
                                   direction = bar_direction))
  return(sf_plot)
}

p_aod_sf <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "aod_ave", main = "Observed",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)
p_lk_sf <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "lk_ave", main = "LatticeKrig",
                       point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                       barwidth = 15, barheight = 1)
p_rf_sf <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "rf_ave", main = "Random Forest",
                       point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                       barwidth = 15, barheight = 1)

p_ave_sf <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                       var_name = "ave_ave", main = "LK and RF average",
                       point_size = 0.3, title_size = 16,
                       val_range = all_range, bar_direction = "horizontal",
                       barwidth = 15, barheight = 1)

p_sl1_sf <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "sl1_ave", main = "SL: Overall",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

p_sl2_sf <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "sl2_ave", main = "SL: Daily",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

p_sl3_sf <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "sl3_ave", main = "SL: Distance-based",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

sf_plot <- 
  ggarrange(p_aod_sf, p_lk_sf, p_rf_sf, 
            font.label = list(size = 14),
            ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggsave(file = "./aod_holdouts/paper_plots/sl_plots/average_withstates_8x11p5.png",
       sf_plot, width = 8, height = 11.5)
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/average_withstates_5x10.png",
       sf_plot, width = 5, height = 10)



# Get full panel with all methods
sf_plot_all <- 
  ggarrange(p_aod_sf, p_lk_sf, p_rf_sf, p_ave_sf, p_sl1_sf, p_sl2_sf, p_sl3_sf,
            font.label = list(size = 14),
            ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom")
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/allaverage_withstates_8x11p5.png",
       sf_plot_all, width = 8, height = 11.5)
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/allaverage_withstates_5x10.png",
       sf_plot_all, width = 5, height = 10)

# Repeat with letters as titles -- can try to use the ggarrange labeling
p_aod_sf2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "aod_ave", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)
p_lk_sf2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                       var_name = "lk_ave", main = "",
                       point_size = 0.3, title_size = 16,
                       val_range = all_range, bar_direction = "horizontal",
                       barwidth = 15, barheight = 1)
p_rf_sf2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                       var_name = "rf_ave", main = "",
                       point_size = 0.3, title_size = 16,
                       val_range = all_range, bar_direction = "horizontal",
                       barwidth = 15, barheight = 1)

p_ave_sf2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "ave_ave", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

p_sl1_sf2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "sl1_ave", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

p_sl2_sf2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "sl2_ave", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

p_sl3_sf2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "sl3_ave", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = all_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)
sf_plot_all2 <- 
  ggarrange(p_aod_sf2, p_lk_sf2, p_rf_sf2, p_ave_sf2, p_sl1_sf2, p_sl2_sf2, p_sl3_sf2,
            font.label = list(size = 14), labels = list("(a)", "(b)", "(c)", "(d)", "(e)",
                                                        "(f)", "(g)"),
            ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom")
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/allaverage2_withstates_8x11p5.png",
       sf_plot_all2, width = 8, height = 11.5)

#### Final adjustment -- exclude values > 1 from the plot for clarity ####
overall_sf2 <- overall_sf %>%
  mutate(aod_ave_m = replace(aod_ave, which(aod_ave > 1), NA))
m_aod_range <- 
  range(c(overall_sf2$aod_ave_m, overall_sf2$ave_ave, 
        overall_sf2$lk_ave, overall_sf2$rf_ave, overall_sf2$sl1_ave,
        overall_sf2$sl2_ave, overall_sf2$sl3_ave), na.rm = T)

p_aod_sf3 <- make_sfplot(main_df = overall_sf2, state_df = states_lcc_sf,
                         var_name = "aod_ave_m", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = m_aod_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)
p_lk_sf3 <- make_sfplot(main_df = overall_sf2, state_df = states_lcc_sf,
                        var_name = "lk_ave", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = m_aod_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)
p_rf_sf3 <- make_sfplot(main_df = overall_sf2, state_df = states_lcc_sf,
                        var_name = "rf_ave", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = m_aod_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

p_ave_sf3 <- make_sfplot(main_df = overall_sf2, state_df = states_lcc_sf,
                         var_name = "ave_ave", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = m_aod_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)

p_sl1_sf3 <- make_sfplot(main_df = overall_sf2, state_df = states_lcc_sf,
                         var_name = "sl1_ave", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = m_aod_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)

p_sl2_sf3 <- make_sfplot(main_df = overall_sf2, state_df = states_lcc_sf,
                         var_name = "sl2_ave", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = m_aod_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)

p_sl3_sf3 <- make_sfplot(main_df = overall_sf2, state_df = states_lcc_sf,
                         var_name = "sl3_ave", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = m_aod_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)
sf_plot_all3 <- 
  ggarrange(p_aod_sf3, p_lk_sf3, p_rf_sf3, p_ave_sf3, 
            p_sl1_sf3, p_sl2_sf3, p_sl3_sf3,
            font.label = list(size = 14), 
            labels = list("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"),
            ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom")
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/allaverage3_withstates_8x11p5.png",
       sf_plot_all3, width = 8, height = 11.5)


# Try wide format
p_aod_sfw <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "aod_ave", main = "Observed",
                        point_size = 0.3, title_size = 12,
                        val_range = all_range, bar_direction = "vertical",
                        barwidth = 1, barheight = 10)
p_lk_sfw <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                       var_name = "lk_ave", main = "LatticeKrig",
                       point_size = 0.3, title_size = 12,
                       val_range = all_range, bar_direction = "vertical",
                       barwidth = 1, barheight = 10)
p_rf_sfw <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                       var_name = "rf_ave", main = "Random Forest",
                       point_size = 0.3, title_size = 12,
                       val_range = all_range, bar_direction = "vertical",
                       barwidth = 1, barheight = 10)
sf_plot_wide <- 
  ggarrange(p_aod_sfw, p_lk_sfw, p_rf_sfw, 
            font.label = list(size = 12),
            ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/average_wide.png",
       sf_plot_wide, width = 8.5, height = 3)

# Plot with just observed and SL: distance-based
p_aod_sfw2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                         var_name = "aod_ave", main = "Observed",
                         point_size = 0.3, title_size = 12,
                         val_range = m_aod_range, bar_direction = "vertical",
                         barwidth = 1, barheight = 10)

p_sl3_sfw2 <- make_sfplot(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "sl3_ave", main = "SL: Distance-based",
                        point_size = 0.3, title_size = 12,
                        val_range = m_aod_range, bar_direction = "vertical",
                        barwidth = 1, barheight = 10)
sf_plot_2 <- 
  ggarrange(p_aod_sfw2, p_sl3_sfw2, 
            font.label = list(size = 14),
            ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/ave_obs_sl3_7x3.png",
       sf_plot_2, width = 7, height = 3)
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/ave_obs_sl3_9x4.png",
       sf_plot_2, width = 9, height = 4)

#### Differences from observed ####
library(scales)
diff_obs_range <- range(overall_df %>% 
                          dplyr::select(lk_aod_diff:sl3_aod_diff),
                        na.rm = T)

make_sfplot_diff <- function(main_df, state_df, barwidth = 1, barheight = 30,
                             bar_direction = "vertical",
                             var_name, main, val_range, title_size = 18,
                             point_size = 0.4, border_size = 0.3) {
  sf_plot <- ggplot() +
    geom_sf(data = main_df, 
            aes_string(color = var_name), 
            shape = 15, size = point_size) + 
    geom_sf(data = st_geometry(state_df), color = "black", 
            size = border_size, alpha = 0) +
    # To plot projected rather than lat/lon
    coord_sf(datum = NULL) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = title_size, face = "bold"),
          panel.background = element_blank()) + 
    ggtitle(main) +
    scale_color_gradient2(limits = val_range, oob = squish) + 
    guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight,
                                   direction = bar_direction))
  return(sf_plot)
}

squish_range <- c(-0.3, 0.3)

p_lk_sfd <- make_sfplot_diff(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "lk_aod_diff", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = squish_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)
p_rf_sfd <- make_sfplot_diff(main_df = overall_sf, state_df = states_lcc_sf,
                        var_name = "rf_aod_diff", main = "",
                        point_size = 0.3, title_size = 16,
                        val_range = squish_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1)

p_ave_sfd <- make_sfplot_diff(main_df = overall_sf, state_df = states_lcc_sf,
                         var_name = "ave_aod_diff", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = squish_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)

p_sl1_sfd <- make_sfplot_diff(main_df = overall_sf, state_df = states_lcc_sf,
                         var_name = "sl1_aod_diff", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = squish_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)

p_sl2_sfd <- make_sfplot_diff(main_df = overall_sf, state_df = states_lcc_sf,
                         var_name = "sl2_aod_diff", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = squish_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)

p_sl3_sfd <- make_sfplot_diff(main_df = overall_sf, state_df = states_lcc_sf,
                         var_name = "sl3_aod_diff", main = "",
                         point_size = 0.3, title_size = 16,
                         val_range = squish_range, bar_direction = "horizontal",
                         barwidth = 15, barheight = 1)

sf_plot_diff <- 
  ggarrange(p_lk_sfd, p_rf_sfd, p_ave_sfd, p_sl1_sfd, p_sl2_sfd, p_sl3_sfd,
            font.label = list(size = 14), labels = list("(a)", "(b)", "(c)", 
                                                        "(d)", "(e)", "(f)"),
            ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/alldiff_withstates.png",
       sf_plot_diff, width = 8*0.8, height = 11.5*0.8)
ggsave(file = "./aod_holdouts/paper_plots/sl_plots/alldiff_withstates_8x11.png",
       sf_plot_diff, width = 8, height = 11.5)

#### 2: Difference of LK/RF averages ####
(plot_avediff <- 
  ggplot(data=overall_df) +
  geom_point(aes(x = cmaq_x, y = cmaq_y, color = lk_rf_diff),
             shape = 15, size = 0.6) + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  ggtitle(label = "") +
  scale_color_gradient2() + 
  coord_fixed())

ggsave(paste0("./aod_holdouts/paper_plots/sl_plots/diff_lkrf.png"), 
       plot_avediff)

(plot_avediff2 <- 
  ggplot(data=overall_df) +
  geom_point(aes(x = cmaq_x, y = cmaq_y, color = lk_rf_diff),
             shape = 15, size = 0.6) + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.height = unit(2, "cm")) + 
  ggtitle(label = "") +
  scale_color_gradient2() + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 30)))
ggsave(paste0("./aod_holdouts/paper_plots/sl_plots/diff_lkrf2.png"), 
       plot_avediff2)

## Try sf version ##

(plot_avediff2_sf <- ggplot() +
  geom_sf(data = overall_sf, 
          aes(color = lk_rf_diff), 
          shape = 15, size = 0.5) + 
  geom_sf(data = st_geometry(states_lcc_sf), color = "grey", 
          size = 0.3, alpha = 0) +
  # To plot projected rather than lat/lon
  coord_sf(datum = NULL) + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.background = element_blank()) + 
  scale_color_gradient2() + 
  guides(color = guide_colourbar(barwidth = 1, barheight = 20)))
ggsave(paste0("./aod_holdouts/paper_plots/sl_plots/diff_lkrf2_states.png"), 
       plot_avediff2_sf, width = 8, height = 5)
# Stopping point

#### 3: Daily (LK - RF) differences and daily AOD/LK/RF Plots ####
diff_range <- range(lkrf_matrix)
aod_range <- range(
  c(range(aod_matrix, na.rm = T),
    range(lk_matrix), range(rf_matrix)))
for (pi in seq(ncol(lk_matrix))) {
  print(range(lk_matrix[, pi]))
}

library(scales)
for (day in days) {
  print(day)
  day_name <- paste0("day", day)
  cur_df <- df_list[[day_name]]
  
  # AOD, LK, RF plots -- common range
  # -- Create daily panel of 3 plots
  aod_value <- cur_df$aod_value
  lk_pred <- cur_df$lk_pred
  rf_pred <- cur_df$rf_pred
  d_aod <- make_plot1(cur_df, aod_value, main = paste0("Day ", day, ": Observed"), 
                      val_range = aod_range, title_size = 8)
  d_lk <- make_plot1(cur_df, lk_pred, main = "LatticeKrig", 
                     val_range = aod_range, title_size = 8)
  d_rf <- make_plot1(cur_df, rf_pred, main = "Random Forest", 
                     val_range = aod_range, title_size = 8)
  
  daily_plot <- 
    ggarrange(d_aod, d_lk, d_rf, 
              ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
  ggsave(file = paste0("./aod_holdouts/paper_plots/sl_plots/daily/day_", day, ".png"),
         daily_plot, width = 2.5, height = 5.5, units = "in")
  rm(aod_value, lk_pred, rf_pred)
  
  # Difference plots
  ## 1: No truncation of values
  plot_diff1 <- 
    ggplot(data=cur_df) +
    geom_point(aes(x = cmaq_x, y = cmaq_y, color = lkrf_diff),
               shape = 15, size = 0.6) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    ggtitle(paste0("Day ", day)) +
    scale_color_gradient2(limits = diff_range) +
    coord_fixed()
  ggsave(paste0("./aod_holdouts/paper_plots/sl_plots/daily/diff_maxrange_", day, ".png"), 
         plot_diff1, width = 5.3, height = 3.95)
  
  ## 2: Keep range between -0.4 and 0.4 -- add note about this. 
  plot_diff2 <- 
    ggplot(data=cur_df) +
    geom_point(data=cur_df, 
               aes(x = cmaq_x, y = cmaq_y, color = lkrf_diff),
               shape = 15, size = 0.6) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    labs(title = paste0("Day ", day), 
         caption = "Values outside range are truncated") +
    scale_color_gradient2(limits = c(-0.4, 0.4),
                          oob = squish) +
    coord_fixed()
  ggsave(paste0("./aod_holdouts/paper_plots/sl_plots/daily/diff_limrange_", day, ".png"), 
         plot_diff2, width = 5.3, height = 3.95)
  
}

#### Construct daily plots with all methods + state borders ####
# -- Collect range of AOD values for all methods on all days
daily_aod_range <- 
  range(c(aod_matrix, lk_matrix, rf_matrix, sl1_matrix, sl2_matrix, sl3_matrix),
        na.rm = T)

# -- Loop through days
for (day in days) {

  print(day)
  day_name <- paste0("day", day)
  
  # Read in the daily data.frame and attach to the sf data.frame
  daily_df <- left_join(x, df_list[[day_name]],
                        by = c("cmaq_id" = "cmaq_id"))
  
  # Create sf plots using your make_sfplot function
  pd_aod <- make_sfplot(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "aod_value", 
                        main = paste0("Observed: July ", day - 182 + 1, ", 2011"),
                        point_size = 0.4, title_size = 10,
                        val_range = daily_aod_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  pd_lk <- make_sfplot(main_df = daily_df, state_df = states_lcc_sf,
                       var_name = "lk_pred", main = "LK",
                       point_size = 0.4, title_size = 10,
                       val_range = daily_aod_range, bar_direction = "horizontal",
                       barwidth = 15, barheight = 1, border_size = 0.2)
  pd_rf <- make_sfplot(main_df = daily_df, state_df = states_lcc_sf,
                       var_name = "rf_pred", main = "RF",
                       point_size = 0.4, title_size = 10,
                       val_range = daily_aod_range, bar_direction = "horizontal",
                       barwidth = 15, barheight = 1, border_size = 0.2)
  pd_ave <- make_sfplot(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "ave_pred", main = "LK and RF average",
                        point_size = 0.4, title_size = 10,
                        val_range = daily_aod_range, bar_direction = "horizontal",
                          barwidth = 15, barheight = 1, border_size = 0.2)
  pd_sl1 <- make_sfplot(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "sl1_pred", main = "SL: Overall",
                        point_size = 0.4, title_size = 10,
                        val_range = daily_aod_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  pd_sl2 <- make_sfplot(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "sl2_pred", main = "SL: Daily",
                        point_size = 0.4, title_size = 10,
                        val_range = daily_aod_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  pd_sl3 <- make_sfplot(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "sl3_pred", main = "SL: Distance-based",
                        point_size = 0.4, title_size = 10,
                        val_range = daily_aod_range, bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)

  daily_plot_all <- 
    ggarrange(pd_aod, pd_lk, pd_rf, pd_ave, pd_sl1, pd_sl2, pd_sl3,
              ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom")
  ggsave(
    file = paste0("./aod_holdouts/paper_plots/sl_plots/daily/day_all_", 
                  day, ".png"),
     daily_plot_all, width = 8, height = 11.5, scale = 0.6)
  ggsave(
    file = paste0("./aod_holdouts/paper_plots/sl_plots/daily/jpeg/day_all_", 
                  day, ".jpeg"),
    daily_plot_all, width = 8, height = 11.5, scale = 0.6)
  
}


#### Re-do, but truncate observations to (-0.05, 1) ####
library(scales)
make_sfplot_squish <- function(main_df, state_df, barwidth = 1, barheight = 30,
                        bar_direction = "vertical",
                        var_name, main, val_range = c(-0.05, 1), title_size = 18,
                        point_size = 0.4, border_size = 0.3) {
  sf_plot <- ggplot() +
    geom_sf(data = main_df, 
            aes_string(color = var_name), 
            shape = 15, size = point_size) + 
    geom_sf(data = st_geometry(state_df), color = "black", 
            size = border_size, alpha = 0) +
    # To plot projected rather than lat/lon
    coord_sf(datum = NULL) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = title_size, face = "bold"),
          panel.background = element_blank()) + 
    ggtitle(main) +
    scale_color_viridis_c(limits = val_range, option = "plasma",
                          oob = squish) +
    guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight,
                                   direction = bar_direction))
  return(sf_plot)
}

# -- Loop through days
for (day in days) {
  
  print(day)
  day_name <- paste0("day", day)
  
  # Read in the daily data.frame and attach to the sf data.frame
  daily_df <- left_join(x, df_list[[day_name]],
                        by = c("cmaq_id" = "cmaq_id"))
  
  # Create sf plots using your make_sfplot function
  pd_aod <- make_sfplot_squish(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "aod_value", 
                        main = paste0("Observed: July ", day - 182 + 1, ", 2011"),
                        point_size = 0.4, title_size = 10,
                        val_range = c(-0.05, 1), bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  pd_lk <- make_sfplot_squish(main_df = daily_df, state_df = states_lcc_sf,
                       var_name = "lk_pred", main = "LK",
                       point_size = 0.4, title_size = 10,
                       val_range = c(-0.05, 1), bar_direction = "horizontal",
                       barwidth = 15, barheight = 1, border_size = 0.2)
  pd_rf <- make_sfplot_squish(main_df = daily_df, state_df = states_lcc_sf,
                       var_name = "rf_pred", main = "RF",
                       point_size = 0.4, title_size = 10,
                       val_range = c(-0.05, 1), bar_direction = "horizontal",
                       barwidth = 15, barheight = 1, border_size = 0.2)
  pd_ave <- make_sfplot_squish(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "ave_pred", main = "LK and RF average",
                        point_size = 0.4, title_size = 10,
                        val_range = c(-0.05, 1), bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  pd_sl1 <- make_sfplot_squish(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "sl1_pred", main = "SL: Overall",
                        point_size = 0.4, title_size = 10,
                        val_range = c(-0.05, 1), bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  pd_sl2 <- make_sfplot_squish(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "sl2_pred", main = "SL: Daily",
                        point_size = 0.4, title_size = 10,
                        val_range = c(-0.05, 1), bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  pd_sl3 <- make_sfplot_squish(main_df = daily_df, state_df = states_lcc_sf,
                        var_name = "sl3_pred", main = "SL: Distance-based",
                        point_size = 0.4, title_size = 10,
                        val_range = c(-0.05, 1), bar_direction = "horizontal",
                        barwidth = 15, barheight = 1, border_size = 0.2)
  
  daily_plot_all <- 
    ggarrange(pd_aod, pd_lk, pd_rf, pd_ave, pd_sl1, pd_sl2, pd_sl3,
              ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom")

  ggsave(
    file = paste0("./aod_holdouts/paper_plots/sl_plots/daily/jpeg/day_all2_", 
                  day, ".jpeg"),
    daily_plot_all, width = 8, height = 11.5, scale = 0.6, dpi = 300)
  
}


#### Construct difference plots with state borders added for each day ####
daily_diff_range <- range(lkrf_matrix)
library(scales)

for (day in days) {
  
  print(day)
  day_name <- paste0("day", day)
  
  # Read in the daily data.frame and attach to the sf data.frame
  daily_df <- left_join(x, df_list[[day_name]],
                        by = c("cmaq_id" = "cmaq_id"))
  
  # Create sf plots using your make_sfplot function
  pd_diff <- make_sfplot_diff(main_df = daily_df, state_df = states_lcc_sf,
                              var_name = "lkrf_diff", main = "",
                              point_size = 0.4, title_size = 10,
                              val_range = c(-0.4, 0.4), 
                              bar_direction = "vertical",
                              barwidth = 1, barheight = 8,
                              border_size = 0.2) + 
    labs(title = paste0("July ", day - 182 + 1, ", 2011"), 
         caption = "Values outside range (-0.4, 0.4) are truncated")
  
  ggsave(
    file = paste0("./aod_holdouts/paper_plots/sl_plots/daily/lkrfdiff_", 
                  day, ".png"),
    pd_diff, width = 5.3, height = 3.95)
  
  ggsave(
    file = paste0("./aod_holdouts/paper_plots/sl_plots/daily/jpeg/lkrfdiff_", 
                  day, ".jpeg"),
    pd_diff, width = 5.3, height = 3.95)
  
}
