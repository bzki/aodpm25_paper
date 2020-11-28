# Create full maps on random forest trained on all observed data
# Also obtain tables of variable importance
# For now, focus is on model (a) -- straightforward random forest with
#  convolution layer.
# In future, other tables/maps can be created and put into Supplemental

library(dplyr)
library(ggplot2)
library(ggpubr)
sessionInfo()

# These are the cross-validation summary tables to help decide mtry
# Focus on model (a)
st_results <- readRDS("./pm25/cvresults/cvbest_st.rds")
print(matrix(st_results$descriptions, ncol = 1))
print(st_results$mgrid)
# Both mtry = 4 an 8 are preferred, depending on the setting and feature set.
# I will provide both: mtry = 4 and mtry = 8. One can go in the supplemental.
# Confusingly, the cross-validation model (a) correspond to full model (b)

# Read in the predictions
pred_list <- readRDS("./pm25/rf_fullpred.rds")
pred_list2 <- readRDS("./pm25/rf_fullpred_2.rds")
df_pred <- pred_list$df_save %>%
  left_join(pred_list2$df_save %>% select(cmaq_id, day, p5a_1:p5c_4),
            by = c("cmaq_id" = "cmaq_id", "day" = "day")) %>%
  select(cmaq_id:conv_pm25, contains(c("b_1", "b_2")))

# Create partial predictions -- use observed where we have it. 
df_ppred <- df_pred %>%
  mutate(
    p1b_1p = if_else(is.na(pm25_value), p1b_1, pm25_value),
    p2b_1p = if_else(is.na(pm25_value), p2b_1, pm25_value),
    p3b_1p = if_else(is.na(pm25_value), p3b_1, pm25_value),
    p4b_1p = if_else(is.na(pm25_value), p4b_1, pm25_value),
    p5b_1p = if_else(is.na(pm25_value), p5b_1, pm25_value),
    p1b_2p = if_else(is.na(pm25_value), p1b_2, pm25_value),
    p2b_2p = if_else(is.na(pm25_value), p2b_2, pm25_value),
    p3b_2p = if_else(is.na(pm25_value), p3b_2, pm25_value),
    p4b_2p = if_else(is.na(pm25_value), p4b_2, pm25_value),
    p5b_2p = if_else(is.na(pm25_value), p5b_2, pm25_value)
  )

sd(df_ppred$p3b_1p - df_ppred$p1b_1p)
sd(df_ppred$p3b_1p - df_ppred$p2b_1p)
sd(df_ppred$p3b_1p - df_ppred$p4b_1p)

#### Get average PM2.5 and prediction for each location across days ####
names(df_ppred)
df_mean <- df_pred %>% filter(day == 182L) %>% select(cmaq_id, cmaq_x, cmaq_y)
library(tidyr)
# Names to summarize
sum_names <- c("pm25_value", 
               paste0("p", seq(5), "b_1p"), 
               paste0("p", seq(5), "b_2p"))

for (nm in sum_names) {
  print(nm)
  wide_df <- pivot_wider(
    data = df_ppred[, c("cmaq_id", "day", nm)],
    id_cols = c(cmaq_id, day),
    names_from = day,
    values_from = all_of(nm),
    names_prefix = nm
  )
  # Create row means
  testthat::expect_identical(wide_df$cmaq_id, df_mean$cmaq_id)
  df_mean[, nm] <- rowMeans(wide_df %>% select(-cmaq_id), na.rm = T)
}

#### Show maps of averages ####
plot_range <- range(df_mean %>% select(pm25_value:p5b_2p), na.rm = T)

(p_pm25 <- ggplot(data = df_mean) +
  geom_point(aes(x = cmaq_x, y = cmaq_y, color = pm25_value),
             shape = 15, size = 2) + 
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
  scale_color_viridis_c(limits = plot_range, option = "plasma") +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 20)) +
  coord_fixed()
)

# loop through and create the plots
# https://stackoverflow.com/a/4856950/12265198
# Save the plots in a list named as variable name
# Save with default settings using ggsave for now
pm25_plotlist <- vector("list", length(sum_names))
names(pm25_plotlist) <- sum_names
temp_labels <- c(
  "(1) No AOD as predictor",
  "(2) GEOS-Chem as predictor",
  "(3) Imputed AOD as predictor",
  "(4) AOD/GEOS-Chem as predictor",
  "(5) Training with observed AOD")
descriptions <- 
  c("Observed", 
    paste0(temp_labels, " (mtry = 4)"),
    paste0(temp_labels, " (mtry = 8)"))
size_vector <- c(3, rep(0.6, length(sum_names) - 1))

for (var_name in sum_names) {
  vi <- which(sum_names == var_name)
  print(var_name)
  print(p <- ggplot(data = df_mean) +
     geom_point(aes_string(x = "cmaq_x", y = "cmaq_y", color = var_name),
                shape = 15, size = size_vector[vi]) + 
     theme(axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks = element_blank(),
           legend.title = element_blank(),
           plot.title = element_text(hjust = 0.5),
           panel.background = element_blank(),
           legend.text = element_text(size = 12)) +
     ggtitle(descriptions[vi]) +
     scale_color_viridis_c(limits = plot_range, option = "plasma") +
     guides(colour = guide_colourbar(barwidth = 1, barheight = 20)) +
     coord_fixed()
  )
  ggsave(paste0("./pm25/fullmaps/", var_name, ".png"), p)
  pm25_plotlist[[var_name]] <- p
  rm(p)
}

# Show maps of differences -- c vs. a, b, d
df_mean <- df_mean %>%
  mutate(
    diffb_3m1_mtry4 = p3b_1p - p1b_1p,
    diffb_3m2_mtry4 = p3b_1p - p2b_1p,
    diffb_3m4_mtry4 = p3b_1p - p4b_1p,
    diffb_3m5_mtry4 = p3b_1p - p5b_1p,
    diffb_3m1_mtry8 = p3b_2p - p1b_2p,
    diffb_3m2_mtry8 = p3b_2p - p2b_2p,
    diffb_3m4_mtry8 = p3b_2p - p4b_2p,
    diffb_3m5_mtry8 = p3b_2p - p5b_2p
  )

summary(df_mean$diffb_3m1_mtry4)
summary(df_mean$diffb_3m2_mtry4)
summary(df_mean$diffb_3m4_mtry4)
summary(df_mean$diffb_3m5_mtry4)
summary(df_mean$diffb_3m1_mtry8)
summary(df_mean$diffb_3m2_mtry8)
summary(df_mean$diffb_3m4_mtry8)
summary(df_mean$diffb_3m5_mtry8)

names_2 <- colnames(df_mean)[which(colnames(df_mean) == "diffb_3m1_mtry4"):ncol(df_mean)]
descriptions_2 <-
  rep(
  c("No AOD/GEOS-Chem",
    "GEOS-Chem",
    "AOD/GEOS-Chem combination",
    "Observed AOD"),
  2)
diff_range <- range(df_mean[, names_2])
diff_plotlist <- vector('list', length(names_2))
names(diff_plotlist) <- names_2

for (var_name in names_2) {
  vi <- which(names_2 == var_name)
  print(var_name)
  print(p <- ggplot(data = df_mean) +
          geom_point(aes_string(x = "cmaq_x", y = "cmaq_y", color = var_name),
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
          #ggtitle("") +
          ggtitle(descriptions_2[vi]) +
          scale_color_gradient2(limits = diff_range) +
          guides(colour = guide_colourbar(barwidth = 1, barheight = 20)) +
          coord_fixed()
  )
  ggsave(paste0("./pm25/fullmaps/", var_name, ".png"), p)
  diff_plotlist[[var_name]] <- p
  rm(p)
}

library(ggpubr)
print(p3_mtry4 <- ggplot(data = df_mean) +
        geom_point(aes_string(x = "cmaq_x", y = "cmaq_y", color = "p3b_1p"),
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
        ggtitle("") +
        scale_color_viridis_c(limits = plot_range, option = "plasma") +
        guides(colour = guide_colourbar(barwidth = 1, barheight = 20)) +
        coord_fixed()
)
print(p3_mtry8 <- ggplot(data = df_mean) +
        geom_point(aes_string(x = "cmaq_x", y = "cmaq_y", color = "p3b_2p"),
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
        ggtitle("") +
        scale_color_viridis_c(limits = plot_range, option = "plasma") +
        guides(colour = guide_colourbar(barwidth = 1, barheight = 20)) +
        coord_fixed()
)
ggsave("./pm25/fullmaps/mainmodel_mtry4.png",
       p3_mtry4)
ggsave("./pm25/fullmaps/mainmodel_mtry8.png",
       p3_mtry8)

# Add geom_sf version with state boundaries
library(sf)
x <- readRDS("cmaq/region_assignment.rds")
str(x)
st_crs(x)
# Attach the lat/lon to your dataset coordinates
df_mean_sf <- left_join(x, df_mean, by = c("cmaq_id" = "cmaq_id"))

# Read in a U.S states shapefile and convert to USLCC
us_states <- read_sf("./cmaq/shapefile_tl2010_us_state_2010/US_state_2010.shp",
                     stringsAsFactors = F) %>%
  filter(!(NAME10 %in% c("Alaska", "Hawaii", "Puerto Rico")))
us_states_uslcc <- st_transform(us_states, crs = st_crs(x))

make_sfplot <- function(main_df, state_df, barwidth = 1, barheight = 30,
                        var_name, main = "", val_range = NULL, title_size = 18,
                        point_size = 0.4, border_size = 0.3, border_color = "black") {
  if (is.null(val_range)) {
    val_range <- range(main_df[[var_name]])
  }
  sf_plot <- ggplot() +
    geom_sf(data = main_df, 
            aes_string(color = var_name), 
            shape = 15, size = point_size) + 
    geom_sf(data = st_geometry(state_df), color = border_color, 
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
    # Note change here from "plasma" to "viridis"
    scale_color_viridis_c(limits = val_range, option = "viridis") +
    guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight))
  
  return(sf_plot)
}

p3_mtry4_sf <- make_sfplot(main_df=df_mean_sf, state_df=us_states_uslcc,
                           var_name="p3b_1p", border_color = "black",
                           barheight = 20)
p3_mtry8_sf <- make_sfplot(main_df=df_mean_sf, state_df=us_states_uslcc,
                           var_name="p3b_2p", border_color = "black",
                           barheight = 20)
ggsave("./pm25/fullmaps/mainmodel_mtry4_states.png",
       p3_mtry4_sf, width = 8, height = 6)
ggsave("./pm25/fullmaps/mainmodel_mtry8_states.png",
       p3_mtry8_sf, width = 8, height = 6)

(diff_plot_mtry4 <- 
  ggarrange(diff_plotlist$diffb_3m1_mtry4, diff_plotlist$diffb_3m2_mtry4, 
            diff_plotlist$diffb_3m4_mtry4, diff_plotlist$diffb_3m5_mtry4,
            common.legend = T, legend = "right", ncol = 2, nrow = 2))
ggsave("./pm25/fullmaps/diffcompare_mtry4.png", diff_plot_mtry4)
(diff_plot_mtry8 <- 
    ggarrange(diff_plotlist$diffb_3m1_mtry8, diff_plotlist$diffb_3m2_mtry8, 
              diff_plotlist$diffb_3m4_mtry8, diff_plotlist$diffb_3m5_mtry8,
              common.legend = T, legend = "right", ncol = 2, nrow = 2))
ggsave("./pm25/fullmaps/diffcompare_mtry8.png", diff_plot_mtry8)

## Add sf versions of these plots
make_sfplot2 <- function(main_df, state_df, barwidth = 1, barheight = 30,
                        var_name, main = "", val_range = NULL, title_size = 18,
                        point_size = 0.4, border_size = 0.3, border_color = "black") {
  if (is.null(val_range)) {
    val_range <- range(main_df[[var_name]])
  }
  sf_plot <- ggplot() +
    geom_sf(data = main_df, 
            aes_string(color = var_name), 
            shape = 15, size = point_size) + 
    geom_sf(data = st_geometry(state_df), color = border_color, 
            size = border_size, alpha = 0) +
    # To plot projected rather than lat/lon
    coord_sf(datum = NULL) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = title_size),
          panel.background = element_blank(),
          legend.text = element_text(size = 12)) + 
    ggtitle(main) +
    scale_color_gradient2(limits = val_range) +
    guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight))
  
  return(sf_plot)
}

diff_plotlist_sf <- vector('list', length(names_2))
names(diff_plotlist_sf) <- names_2

for (var_name in names_2) {
  vi <- which(names_2 == var_name)
  print(var_name)
  p <- make_sfplot2(main_df=df_mean_sf, state_df=us_states_uslcc,
                    var_name=var_name, border_color = "grey",
                    barheight = 25, main = descriptions_2[vi],
                    title_size = 12, val_range = diff_range,
                    point_size = 0.5)
  diff_plotlist_sf[[var_name]] <- p
  rm(p); gc()
}

diff_plot_mtry4_sf <- 
    ggarrange(diff_plotlist_sf$diffb_3m1_mtry4, diff_plotlist_sf$diffb_3m2_mtry4, 
              diff_plotlist_sf$diffb_3m4_mtry4, diff_plotlist_sf$diffb_3m5_mtry4,
              common.legend = T, legend = "right", ncol = 2, nrow = 2)
ggsave("./pm25/fullmaps/diffcompare_mtry4_states.png", diff_plot_mtry4_sf,
       width = 8, height = 6)

diff_plot_mtry8_sf <- 
    ggarrange(diff_plotlist_sf$diffb_3m1_mtry8, diff_plotlist_sf$diffb_3m2_mtry8, 
              diff_plotlist_sf$diffb_3m4_mtry8, diff_plotlist_sf$diffb_3m5_mtry8,
              common.legend = T, legend = "right", ncol = 2, nrow = 2)
ggsave("./pm25/fullmaps/diffcompare_mtry8_states.png", diff_plot_mtry8_sf,
       width = 8, height = 6)

#### Re-do plots with (a) - (d) labels ####
diff_plotlist_sf2 <- vector('list', length(names_2))
names(diff_plotlist_sf2) <- names_2
use_letters <- rep(letters[1:4], 2)
for (var_name in names_2) {
  vi <- which(names_2 == var_name)
  print(var_name)
  p <- make_sfplot2(main_df=df_mean_sf, state_df=us_states_uslcc,
                    var_name=var_name, border_color = "grey",
                    barheight = 25, main = paste0("(", use_letters[vi], ")"),
                    title_size = 12, val_range = diff_range,
                    point_size = 0.5)
  diff_plotlist_sf2[[var_name]] <- p
  rm(p); gc()
}

diff_plot_mtry4_sf2 <- 
  ggarrange(diff_plotlist_sf2$diffb_3m1_mtry4, diff_plotlist_sf2$diffb_3m2_mtry4, 
            diff_plotlist_sf2$diffb_3m4_mtry4, diff_plotlist_sf2$diffb_3m5_mtry4,
            common.legend = T, legend = "right", ncol = 2, nrow = 2)
ggsave("./pm25/fullmaps/diffcompare_mtry4_states2.png", diff_plot_mtry4_sf2,
       width = 8, height = 6)

diff_plot_mtry8_sf2 <- 
  ggarrange(diff_plotlist_sf2$diffb_3m1_mtry8, diff_plotlist_sf2$diffb_3m2_mtry8, 
            diff_plotlist_sf2$diffb_3m4_mtry8, diff_plotlist_sf2$diffb_3m5_mtry8,
            common.legend = T, legend = "right", ncol = 2, nrow = 2)
ggsave("./pm25/fullmaps/diffcompare_mtry8_states2.png", diff_plot_mtry8_sf2,
       width = 8, height = 6)

#### ST Model: Daily PM2.5 Observed and Difference Plots ####
# -- save these as jpeg to conserve some space.
# st models: create both mtry = 4 and mtry = 8 panels
# panels of 5 figures
library(scales)
make_sfplot3 <- function(main_df, state_df, barwidth = 1, barheight = 30,
                         var_name, main = "", val_range = NULL, title_size = 18,
                         point_size = 0.4, border_size = 0.3, border_color = "black") {
  sf_plot <- ggplot() +
    geom_sf(data = main_df, 
            aes_string(color = var_name), 
            shape = 15, size = point_size) + 
    geom_sf(data = st_geometry(state_df), color = border_color, 
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
    scale_color_viridis_c(limits = val_range, option = "viridis",
                          oob = squish) +
    guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight))
  
  return(sf_plot)
}

make_sfplot3diff <- function(main_df, state_df, barwidth = 1, barheight = 30,
                         var_name, main = "", val_range = NULL, title_size = 18,
                         point_size = 0.4, border_size = 0.3, border_color = "black",
                         pm25_col = "green4", pm25_size = 0.3) {

  pm25_tmp <- main_df %>%
    filter(!is.na(pm25_value))
  
  sf_plot <- ggplot() +
    geom_sf(data = main_df, 
            aes_string(color = var_name), 
            shape = 15, size = point_size) + 
    geom_sf(data = st_geometry(state_df), color = border_color, 
            size = border_size, alpha = 0) +
    geom_sf(data = st_geometry(pm25_tmp), color = pm25_col,
            size = pm25_size, shape = 20) +
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
    guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight))
  
  return(sf_plot)
}

days = 182:212
for (day in days) {
  
  print(day)
  daily_sf <- left_join(x, df_ppred[df_ppred$day == day, ], 
                           by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(
      diffb_3m1_mtry4 = p3b_1p - p1b_1p,
      diffb_3m2_mtry4 = p3b_1p - p2b_1p,
      diffb_3m4_mtry4 = p3b_1p - p4b_1p,
      diffb_3m5_mtry4 = p3b_1p - p5b_1p,
      diffb_3m1_mtry8 = p3b_2p - p1b_2p,
      diffb_3m2_mtry8 = p3b_2p - p2b_2p,
      diffb_3m4_mtry8 = p3b_2p - p4b_2p,
      diffb_3m5_mtry8 = p3b_2p - p5b_2p
    )
  
  # day_range_p4 <- range(c(daily_sf$p1b_1p, daily_sf$p2b_1p, daily_sf$p3b_1p,
  #                         daily_sf$p4b_1p, daily_sf$p5b_1p))
  # 
  # day_range_p8 <- range(c(daily_sf$p1b_2p, daily_sf$p2b_2p, daily_sf$p3b_2p,
  #                         daily_sf$p4b_2p, daily_sf$p5b_2p))
  
  # day_range_p4[1] <- max(day_range_p4[1], 0)
  # day_range_p4[2] <- min(day_range_p4[2], 35)
  # day_range_p8[1] <- max(day_range_p8[1], 0)
  # day_range_p8[2] <- min(day_range_p8[2], 35)
  
  # Set range to (0, 30) for all days and truncate outlier values
  day_range_p4 <- c(0, 30)
  day_range_p8 <- c(0, 30)
  
  # Create 10 figures using make_sfplot3 function
  # -- for now, I will refrain from creating a constant scale day to day
  p_m1_p4 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                         var_name="p1b_1p", border_color = "black",
                         barheight = 30, val_range = day_range_p4,
                         title_size = 14,
                         main = paste0("July ", day - 182 + 1, ", 2011: M1"))
  
  p_m2_p4 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                         var_name="p2b_1p", border_color = "black",
                         barheight = 30, val_range = day_range_p4,
                         title_size = 14,
                         main = "M2")
  
  p_m3_p4 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                         var_name="p3b_1p", border_color = "black",
                         barheight = 30, val_range = day_range_p4,
                         title_size = 14,
                         main = "M3")
  
  p_m4_p4 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                         var_name="p4b_1p", border_color = "black",
                         barheight = 30, val_range = day_range_p4,
                         title_size = 14,
                         main = "M4")
  
  p_m5_p4 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                         var_name="p5b_1p", border_color = "black",
                         barheight = 30, val_range = day_range_p4,
                         title_size = 14,
                         main = "M5")
  
  p_p4 <- 
    ggarrange(p_m1_p4, p_m2_p4, p_m3_p4, p_m4_p4, p_m5_p4,
              common.legend = T, legend = "right", ncol = 2, nrow = 3)
  ggsave(paste0("./pm25/fullmaps/jpeg/daily_mtry4_", day, ".jpeg"), 
                p_p4, width = 8.5*0.8, height = 11.5*0.8)
  
  p_m3_p4_x <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                           var_name="p3b_1p", border_color = "black",
                           barheight = 20, val_range = day_range_p4,
                           main = paste0("July ", day - 182 + 1, ", 2011: M3"))
  ggsave(paste0("./pm25/fullmaps/jpeg/daily_m3_mtry4_", day, ".jpeg"), 
         p_m3_p4_x, width = 8, height = 6)
  
  
  # Repeat for mtry = 8
  p_m1_p8 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                          var_name="p1b_2p", border_color = "black",
                          barheight = 30, val_range = day_range_p8,
                          title_size = 14,
                          main = paste0("July ", day - 182 + 1, ", 2011: M1"))
  
  p_m2_p8 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                          var_name="p2b_2p", border_color = "black",
                          barheight = 30, val_range = day_range_p8,
                          title_size = 14,
                          main = "M2")
  
  p_m3_p8 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                          var_name="p3b_2p", border_color = "black",
                          barheight = 30, val_range = day_range_p8,
                          title_size = 14,
                          main = "M3")
  
  p_m4_p8 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                          var_name="p4b_2p", border_color = "black",
                          barheight = 30, val_range = day_range_p8,
                          title_size = 14,
                          main = "M4")
  
  p_m5_p8 <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                          var_name="p5b_2p", border_color = "black",
                          barheight = 30, val_range = day_range_p8,
                          title_size = 14,
                          main = "M5")
  
  p_p8 <- 
    ggarrange(p_m1_p8, p_m2_p8, p_m3_p8, p_m4_p8, p_m5_p8,
              common.legend = T, legend = "right", ncol = 2, nrow = 3)
  ggsave(paste0("./pm25/fullmaps/jpeg/daily_mtry8_", day, ".jpeg"), 
         p_p8, width = 8.5, height = 11.5, dpi = 200)
  
  p_m3_p8_x <- make_sfplot3(main_df=daily_sf, state_df=us_states_uslcc,
                            var_name="p3b_2p", border_color = "black",
                            barheight = 20, val_range = day_range_p8,
                            main = paste0("July ", day - 182 + 1, ", 2011: M3"))
  ggsave(paste0("./pm25/fullmaps/jpeg/daily_m3_mtry8_", day, ".jpeg"), 
         p_m3_p8_x, width = 8, height = 6)
  
  
  ## Differences
  diff_use_range = c(-2, 2)
  p_m3m1_p4 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                          var_name="diffb_3m1_mtry4", border_color = "black",
                          barheight = 20, val_range = diff_use_range,
                          title_size = 14,
                          main = paste0("July ", day - 182 + 1, ", 2011: M3 - M1"))
  p_m3m2_p4 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                                var_name="diffb_3m2_mtry4", border_color = "black",
                                barheight = 20, val_range = diff_use_range,
                                title_size = 14,
                                main = "M3 - M2")
  p_m3m4_p4 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                                var_name="diffb_3m4_mtry4", border_color = "black",
                                barheight = 20, val_range = diff_use_range,
                                title_size = 14,
                                main = "M3 - M4")
  p_m3m5_p4 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                                var_name="diffb_3m5_mtry4", border_color = "black",
                                barheight = 20, val_range = diff_use_range,
                                title_size = 14,
                                main = "M3 - M5")
  
  p_diff_p4 <- 
    ggarrange(p_m3m1_p4, p_m3m2_p4, p_m3m4_p4, p_m3m5_p4,
              common.legend = T, legend = "right", ncol = 2, nrow = 2)
  ggsave(paste0("./pm25/fullmaps/jpeg/diff_mtry4_", day, ".jpeg"), 
         p_diff_p4, width = 6.5, height = 6.5, dpi = 200)
  
  # repeat for mtry = 8
  p_m3m1_p8 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                                var_name="diffb_3m1_mtry8", border_color = "black",
                                barheight = 20, val_range = diff_use_range,
                                title_size = 14,
                                main = paste0("July ", day - 182 + 1, ", 2011: M3 - M1"))
  p_m3m2_p8 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                                var_name="diffb_3m2_mtry8", border_color = "black",
                                barheight = 20, val_range = diff_use_range,
                                title_size = 14,
                                main = "M3 - M2")
  p_m3m4_p8 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                                var_name="diffb_3m4_mtry8", border_color = "black",
                                barheight = 20, val_range = diff_use_range,
                                title_size = 14,
                                main = "M3 - M4")
  p_m3m5_p8 <- make_sfplot3diff(main_df=daily_sf, state_df=us_states_uslcc,
                                var_name="diffb_3m5_mtry8", border_color = "black",
                                barheight = 20, val_range = diff_use_range,
                                title_size = 14,
                                main = "M3 - M5")
  
  p_diff_p8 <- 
    ggarrange(p_m3m1_p8, p_m3m2_p8, p_m3m4_p8, p_m3m5_p8,
              common.legend = T, legend = "right", ncol = 2, nrow = 2)
  ggsave(paste0("./pm25/fullmaps/jpeg/diff_mtry8_", day, ".jpeg"), 
         p_diff_p8, width = 6.5, height = 6.5, dpi = 200)
  
}

################################################################
#### Daily Model: TODO ####


#### Summarize variable importance ####
# structure is: import_list_b[[mtry]][[1=a, 2=b, 3=c, 4=d]]
sort(pred_list$import_list_b[[1]][[1]], decreasing = T)
sort(pred_list$import_list_b[[1]][[2]], decreasing = T)
sort(pred_list$import_list_b[[1]][[3]], decreasing = T)
sort(pred_list$import_list_b[[1]][[4]], decreasing = T)
# Model 5 was run separately
sort(pred_list2$import_list_b[[1]][[1]], decreasing = T)

import_list_b_mtry4 <- pred_list$import_list_b[[1]]
import_list_b_mtry4_m5 <- pred_list2$import_list_b[[1]][[1]]
import_list_b_mtry8 <- pred_list$import_list_b[[2]]
import_list_b_mtry8_m5 <- pred_list2$import_list_b[[2]][[1]]

# Create a single data.frame for all of these
# - Get set of names
library(dplyr)
imp_names <- 
  union(
  union(
    union(names(import_list_b_mtry4[[1]]),
          names(import_list_b_mtry4[[2]])),
    union(names(import_list_b_mtry4[[3]]),
          names(import_list_b_mtry4[[4]]))
  ),
  names(import_list_b_mtry4_m5))

df_1 <- data.frame(M1a = import_list_b_mtry4[[1]])
df_1$Feature <- rownames(df_1)
rownames(df_1) <- NULL
df_2 <- data.frame(M2a = import_list_b_mtry4[[2]])
df_2$Feature <- rownames(df_2)
rownames(df_2) <- NULL
df_3 <- data.frame(M3a = import_list_b_mtry4[[3]])
df_3$Feature <- rownames(df_3)
rownames(df_3) <- NULL
df_4 <- data.frame(M4a = import_list_b_mtry4[[4]])
df_4$Feature <- rownames(df_4)
rownames(df_4) <- NULL
df_5 <- data.frame(M5a = import_list_b_mtry4_m5)
df_5$Feature <- rownames(df_5)
rownames(df_5) <- NULL

# Create a data.frame with these names as columns
imp_df <- data.frame(Feature = imp_names, stringsAsFactors = F) %>%
  left_join(df_1) %>% left_join(df_2) %>% left_join(df_3) %>% 
  left_join(df_4) %>% left_join(df_5)

matrix(imp_df$Feature, ncol = 1)
# [1,] "day"           
# [2,] "cmaq_x"        
# [3,] "cmaq_y"        
# [4,] "elev"          
# [5,] "emissi11"      
# [6,] "forest_cover"  
# [7,] "high"          
# [8,] "limi"          
# [9,] "local"         
# [10,] "is"            
# [11,] "pd"            
# [12,] "nldas_pevapsfc"
# [13,] "nldas_dlwrfsfc"
# [14,] "nldas_dswrfsfc"
# [15,] "nldas_cape"    
# [16,] "nldas_fpcsfc"  
# [17,] "nldas_pcpsfc"  
# [18,] "nldas_rh2m"    
# [19,] "nldas_tmp2m"   
# [20,] "nldas_vgrd10m" 
# [21,] "nldas_ugrd10m" 
# [22,] "nldas_pressfc" 
# [23,] "conv_pm25"     
# [24,] "aod_missing"   
# [25,] "day_of_week"   
# [26,] "gc_aod"        
# [27,] "aod_imputed"   
# [28,] "aod_gc_combine"
# [29,] "aod_value"
Description_Features <-
  c("Day",
    "CMAQ-X Coordinate", "CMAQ-Y Coordinate",
    "Elevation", "EPA 2011 emission inventory",
    "Percent forest cover", "Total length of highway", "Total length of limited-access road",
    "Total length of local road", "Impervious surface (%)", "Population density",
    "Potential evaporation", "Surface DW longwave radiation flux", 
    "Surface DW shortwave radiation flux", "Convective available potential energy",
    "Faction of total precipitation that is convective", 
    "Precipitation", "Relative humidity", "Temperature", "v-direction wind speed",
    "u-direction wind-speed", "Pressure at surface", "Convolution layer PM2.5",
    "AOD Missing Indicator", "Day of the Week", "GEOS-Chem", 
    "Imputed AOD", "AOD/GEOS-Chem combination", "Observed AOD")

# Check for correctness
cbind(imp_df$Feature, Description_Features)

imp_df2 <- cbind(imp_df, Description_Features) %>%
  relocate(Feature, Description_Features)

# Order it in a way that makes sense -- try Model 3a
sort_imp_df <- imp_df2 %>%
  arrange(desc(M3a))

# Write latex using xtable -- use just 2 decimal points
library(xtable)
print(xtable(sort_imp_df %>% select(-Feature), digits = 2), 
      include.rownames = F, 
      file = "./pm25/tables/st_importance_main_mtry4.txt")

# Re-do the above table for mtry8
rm(imp_df, imp_df2, sort_imp_df)
df_1m8 <- data.frame(M1a = import_list_b_mtry8[[1]])
df_1m8$Feature <- rownames(df_1m8)
rownames(df_1m8) <- NULL
df_2m8 <- data.frame(M2a = import_list_b_mtry8[[2]])
df_2m8$Feature <- rownames(df_2m8)
rownames(df_2m8) <- NULL
df_3m8 <- data.frame(M3a = import_list_b_mtry8[[3]])
df_3m8$Feature <- rownames(df_3m8)
rownames(df_3m8) <- NULL
df_4m8 <- data.frame(M4a = import_list_b_mtry8[[4]])
df_4m8$Feature <- rownames(df_4m8)
rownames(df_4m8) <- NULL
df_5m8 <- data.frame(M5a = import_list_b_mtry8_m5)
df_5m8$Feature <- rownames(df_5m8)
rownames(df_5m8) <- NULL

# Create a data.frame with these names as columns
imp_m8_df <- data.frame(Feature = imp_names, stringsAsFactors = F) %>%
  left_join(df_1m8) %>% left_join(df_2m8) %>% left_join(df_3m8) %>% 
  left_join(df_4m8) %>% left_join(df_5m8)

imp_m8_df2 <- cbind(imp_m8_df, Description_Features) %>%
  relocate(Feature, Description_Features)

# Order it in a way that makes sense -- try Model 3a
sort_m8_imp_df <- imp_m8_df2 %>%
  arrange(desc(M3a))

# Write latex using xtable -- use just 2 decimal points
library(xtable)
print(xtable(sort_m8_imp_df %>% select(-Feature), digits = 2), 
      include.rownames = F, 
      file = "./pm25/tables/st_importance_main_mtry8.txt")
