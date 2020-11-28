library(sp)
library(dplyr)
df_sum <- as.data.frame(readRDS("./rawrds/observed_summary.rds"))
(max_day <- df_sum$day[which.max(df_sum$pct_observed)])
(min_day <- df_sum$day[which.min(df_sum$pct_observed)])
# Day 182 and 193 are the max/min days

# 1: Load in the full datasets and the holdouts
full_df <- readRDS("./rawrds/n4_day182.rds") %>%
  select(cmaq_id, cmaq_x, cmaq_y)
day_182 <- readRDS("./rawrds/aod_182.rds") %>%
  select(cmaq_id, aod_value) %>%
  rename(obs_182 = aod_value)
day_193 <- readRDS("./rawrds/aod_193.rds") %>%
  select(cmaq_id, aod_value) %>%
  rename(obs_193 = aod_value)

# Focus is entirely on S4 at this point
holdout_4 <- readRDS("./aod_holdouts/holdout_4_circleB.rds")

day_182$s4_182 <- NA
day_182$s4_182[-holdout_4$day182$holdout_index] <- 
  day_182$obs_182[-holdout_4$day182$holdout_index]

day_193$s4_193 <- NA
day_193$s4_193[-holdout_4$day193$holdout_index] <- 
  day_193$obs_193[-holdout_4$day193$holdout_index]

comb_df <- full_df %>% left_join(day_182) %>% left_join(day_193)
data_df <- comb_df
coordinates(comb_df) <- ~cmaq_x + cmaq_y

library(latticeExtra)
png("./aod_holdouts/paper_plots/general/minmax3_aod.png",
    width = 850, height = 850)
par(cex = 0.1)
print(
  spplot(comb_df, c("obs_182", "obs_193",
                    "s4_182", "s4_193"),
         names.attr = c("Day 182: Full observed", "Day 193: Full observed",
                        "Day 182: Training", "Day 193: Training"),
         colorkey = TRUE, cex = 0.35,
         layout = c(2, 2)) + 
    layer_(sp.points(comb_df, col='grey', cex = 0.2))
)
dev.off()

##################################
#### Simplify and use ggplot2 ####

library(ggplot2)
library(ggpubr)
str(data_df)
make_plot <- function(df, var_name, main, val_range, title_size = 18) {
  
  ggp <- ggplot(data=df) +
    geom_point(aes_string(x = "cmaq_x", y = "cmaq_y", color = var_name),
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

aod_range <- range(c(data_df$obs_182, data_df$obs_193), na.rm = T)
# obs_182 <- data_df$obs_182
# obs_193 <- data_df$obs_193
# s4_182 <- data_df$s4_182
# s4_193 <- data_df$s4_193
o_1 <- make_plot(data_df, "obs_182", val_range = aod_range, title_size = 13,
          main = "July 1, 2011: Observed AOD")
o_2 <- make_plot(data_df, "obs_193", val_range = aod_range, title_size = 13,
          main = "July 12, 2011: Observed AOD")
s_1 <- make_plot(data_df, "s4_182", val_range = aod_range, title_size = 13,
          main = "July 1, 2011: Training")
s_2 <- make_plot(data_df, "s4_193", val_range = aod_range, title_size = 13,
          main = "July 12, 2011: Training")

gg_fig <- 
  ggarrange(o_1, o_2, s_1, s_2, 
            ncol = 2, nrow = 2, common.legend = TRUE,
            legend = "bottom")

ggsave("./aod_holdouts/paper_plots/general/twodays.png", gg_fig)

# Additional modifications
# (1) Make the legend color key bigger
make_plot2 <- function(df, var_name, main, val_range, title_size = 18,
                       point_size = 0.6) {
  
  ggp <- ggplot(data=df) +
    geom_point(aes_string(x = "cmaq_x", y = "cmaq_y", color = var_name),
               shape = 15, size = point_size) + 
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
    guides(color = guide_colourbar(barwidth = 1, barheight = 30)) +
    coord_fixed(ratio = 1) # Important for spatial plots if not using geom_sf
  
  return(ggp)
}

make_plot2(data_df, "obs_182", val_range = aod_range, title_size = 13,
          main = "July 1, 2011: Observed AOD")
o_1b <- make_plot2(data_df, "obs_182", val_range = aod_range, title_size = 13,
                 main = "July 1, 2011: Observed AOD")
o_2b <- make_plot2(data_df, "obs_193", val_range = aod_range, title_size = 13,
                 main = "July 12, 2011: Observed AOD")
s_1b <- make_plot2(data_df, "s4_182", val_range = aod_range, title_size = 13,
                 main = "July 1, 2011: Training")
s_2b <- make_plot2(data_df, "s4_193", val_range = aod_range, title_size = 13,
                 main = "July 12, 2011: Training")
gg_fig2 <- 
  ggarrange(o_1b, o_2b, s_1b, s_2b,
            ncol = 2, nrow = 2, common.legend = TRUE,
            legend = "right")
gg_fig2
ggsave("./aod_holdouts/paper_plots/general/twodays_biglegend_12x8.png", gg_fig2,
       width = 12, height = 8)
ggsave("./aod_holdouts/paper_plots/general/twodays_biglegend_sq.png", gg_fig2,
       width = 8.5, height = 8.5)

# (2) See about adding state outlines (i.e., geom_sf())
library(sf)
x <- readRDS("cmaq/region_assignment.rds")
str(x)
st_crs(x)
# Attach the lat/lon to your dataset coordinates
x2 <- left_join(x, data_df, by = c("cmaq_id" = "cmaq_id"))

# Read in a U.S states shapefile and convert to USLCC
# Now using better shapefile from IPUMS NHGIS for state boundaries
us_states <- read_sf("./cmaq/shapefile_tl2010_us_state_2010/US_state_2010.shp",
                     stringsAsFactors = F) %>%
  filter(!(NAME10 %in% c("Alaska", "Hawaii", "Puerto Rico")))
us_states_uslcc <- st_transform(us_states,
                                crs = st_crs(x))

# Old code using spData -- the borders here are less precise
# library(spData) # See: https://nowosad.github.io/spData/
# us_states_longlat <- spData::us_states # Note this is an sf object
# st_crs(us_states_longlat) # "+proj=longlat +datum=NAD83 +no_defs"
# us_states_uslcc <- st_transform(us_states_longlat,
#                                 crs = st_crs(x))

# Try getting borders
make_sfplot <- function(main_df, state_df, barwidth = 1, barheight = 30,
                        var_name, main, val_range, title_size = 18,
                        point_size = 0.6, border_size = 0.3) {
  
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
  guides(color = guide_colourbar(barwidth = barwidth, barheight = barheight))
  
return(sf_plot)
}

o_1c <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "obs_182", val_range = aod_range, title_size = 13,
                    main = "July 1, 2011: Observed AOD",
                    point_size = 0.3)
o_2c <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "obs_193", val_range = aod_range, title_size = 13,
                    main = "July 12, 2011: Observed AOD",
                    point_size = 0.3)
s_1c <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "s4_182", val_range = aod_range, title_size = 13,
                    main = "July 1, 2011: Training",
                    point_size = 0.3)
s_2c <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "s4_193", val_range = aod_range, title_size = 13,
                    main = "July 12, 2011: Training",
                    point_size = 0.3)
gg_fig_sf <- 
  ggarrange(o_1c, o_2c, s_1c, s_2c,
            ncol = 2, nrow = 2, common.legend = TRUE,
            legend = "right")
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_12x8.png", gg_fig_sf,
       width = 12, height = 8)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_9x6.png", gg_fig_sf,
       width = 9, height = 6)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_sq.png", gg_fig_sf,
       width = 8.5, height = 8.5)

# Create modified version for Atmosphere with (a), (b), (c), (d) headings
gg_fig_sf_2 <- 
  ggarrange(o_1c, o_2c, s_1c, s_2c, labels = list("(a)", "(b)", "(c)", "(d)"),
            ncol = 2, nrow = 2, common.legend = TRUE,
            legend = "right")
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_12x8_2.png", gg_fig_sf_2,
       width = 12, height = 8)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_9x6_2.png", gg_fig_sf_2,
       width = 9, height = 6)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_sq_2.png", gg_fig_sf_2,
       width = 8.5, height = 8.5)

o_1d <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "obs_182", val_range = aod_range, title_size = 15,
                    main = "(a)",
                    point_size = 0.3)
o_2d <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "obs_193", val_range = aod_range, title_size = 15,
                    main = "(b)",
                    point_size = 0.3)
s_1d <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "s4_182", val_range = aod_range, title_size = 15,
                    main = "(c)",
                    point_size = 0.3)
s_2d <- make_sfplot(main_df = x2, state_df = us_states_uslcc, 
                    var_name = "s4_193", val_range = aod_range, title_size = 15,
                    main = "(d)",
                    point_size = 0.3)
gg_fig_sf_atmos <- 
  ggarrange(o_1d, o_2d, s_1d, s_2d,
            ncol = 2, nrow = 2, common.legend = TRUE,
            legend = "right")

ggsave("./aod_holdouts/paper_plots/general/twodays_sf_12x8_atmos.png", 
       gg_fig_sf_atmos,
       width = 12, height = 8)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_9x6_atmos.png", 
       gg_fig_sf_atmos,
       width = 9, height = 6)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_sq_atmos.png", 
       gg_fig_sf_atmos,
       width = 8.5, height = 8.5)

# Updated -- remove values over 1.0 from plot to aid in visual clarity
# -- add to caption in manuscript
x3 <- x2 %>%
  mutate(obs_182m = replace(obs_182, which(obs_182 > 1), NA),
         obs_193m = replace(obs_193, which(obs_193 > 1), NA),
         s4_182m = replace(s4_182, which(s4_182 > 1), NA),
         s4_193m = replace(s4_193, which(s4_193 > 1), NA))

o_1e <- make_sfplot(main_df = x3, state_df = us_states_uslcc, 
                    var_name = "obs_182m", val_range = c(-0.05, 1), title_size = 15,
                    main = "(a)",
                    point_size = 0.3)
o_2e <- make_sfplot(main_df = x3, state_df = us_states_uslcc, 
                    var_name = "obs_193m", val_range = c(-0.05, 1), title_size = 15,
                    main = "(b)",
                    point_size = 0.3)
s_1e <- make_sfplot(main_df = x3, state_df = us_states_uslcc, 
                    var_name = "s4_182m", val_range = c(-0.05, 1), title_size = 15,
                    main = "(c)",
                    point_size = 0.3)
s_2e <- make_sfplot(main_df = x3, state_df = us_states_uslcc, 
                    var_name = "s4_193m", val_range = c(-0.05, 1), title_size = 15,
                    main = "(d)",
                    point_size = 0.3)
gg_fig_sf_atmos2 <- 
  ggarrange(o_1e, o_2e, s_1e, s_2e,
            ncol = 2, nrow = 2, common.legend = TRUE,
            legend = "right")
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_12x8_atmos2.png", 
       gg_fig_sf_atmos2,
       width = 12, height = 8)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_9x6_atmos2.png", 
       gg_fig_sf_atmos2,
       width = 9, height = 6)
ggsave("./aod_holdouts/paper_plots/general/twodays_sf_sq_atmos2.png", 
       gg_fig_sf_atmos2,
       width = 8.5, height = 8.5)

#### Panel plot for every day ####
# Read in every day and record the AOD range
days = 182:212
all_aod_range <- NULL
for (day in days) {
  print(day)
  print(all_aod_range)
  tmp_df <- readRDS(paste0("./rawrds/split_4_", day, ".rds"))
  all_aod_range <- range(all_aod_range, tmp_df$aod_value)
  rm(tmp_df)
}


# Use split_4_XXX.rds files to help do this quickly
# Link to the spatial file

# Save: 
days = 182:212
for (day in days) {
  print(day)
  
  daily_df <- readRDS(paste0("./rawrds/split_4_", day, ".rds")) %>%
    mutate(train = ifelse(split != "Test", 1, 0))
  sf_df <- left_join(x, daily_df, by = c("cmaq_id" = "cmaq_id"))
  sf_train_df <- left_join(x, daily_df %>% filter(train == 1),
                           by = c("cmaq_id" = "cmaq_id"))
  
  daily_obs <- 
    make_sfplot(main_df = sf_df, state_df = us_states_uslcc, 
                var_name = "aod_value", val_range = all_aod_range, title_size = 13,
                main = paste0("July ", (day + 1 - 182), ", 2011: Observed"),
                point_size = 0.3)
  
  daily_train <- 
    make_sfplot(main_df = sf_train_df, state_df = us_states_uslcc, 
                var_name = "aod_value", val_range = all_aod_range, title_size = 13,
                main = paste0("July ", (day + 1 - 182), ", 2011: Training"),
                point_size = 0.3)
  
  gg_daily <- ggarrange(daily_obs, daily_train,
                        ncol = 1, nrow = 2, common.legend = TRUE,
                        legend = "right")
  ggsave(paste0(
    "./aod_holdouts/paper_plots/daily_observed/obstrain_", day, ".png"), 
    gg_daily, width = 8, height = 12)
  
  rm(daily_obs, daily_train, sf_df, sf_train_df); gc()
}

#### Plot (2 colors) for every day explicitly showing location of testing ####
make_sfplot_bin <- function(main_df, state_df,
                            var_name, main, title_size = 18,
                            point_size = 0.6, border_size = 0.3) {
  
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
          legend.text = element_text(size=title_size - 4),
          plot.title = element_text(hjust = 0.5, size = title_size, face = "bold"),
          panel.background = element_blank()) + 
    ggtitle(main) 
  
  return(sf_plot)
}

library(forcats)
days = 182:212
for (day in days) {
  print(day)
  
  daily_df <- readRDS(paste0("./rawrds/split_4_", day, ".rds")) %>%
    mutate(train = as_factor(ifelse(split != "Test", "Training", "Testing")))
  sf_df <- left_join(x, daily_df, by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(train = fct_relevel(train, "Testing", "Training"))

  daily_obs <- 
    make_sfplot_bin(main_df = sf_df, state_df = us_states_uslcc, 
                var_name = "train", title_size = 13,
                main = paste0("July ", (day + 1 - 182), ", 2011: Training and testing"),
                point_size = 0.6) + 
    guides(colour = guide_legend(override.aes = list(size=6)))
  
  ggsave(paste0(
    "./aod_holdouts/paper_plots/daily_observed/bin_obstrain_", day, ".png"), 
    daily_obs, width = 5, height = 4)
  
  ggsave(paste0(
    "./aod_holdouts/paper_plots/daily_observed/jpeg/bin_obstrain_", day, ".jpeg"), 
    daily_obs, width = 5, height = 4)
  
  rm(daily_obs, sf_df); gc()
}
