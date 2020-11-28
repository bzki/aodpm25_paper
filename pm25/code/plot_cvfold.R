# PM25 cross-validation fold visualization
# -- As outputted from pm25_cvfold.R
# Creating a nicer plot of the spatial block cross-validation folds

library(sf)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(forcats)
sessionInfo()

# Read in previously created CMAQ sf file 
x <- readRDS("cmaq/region_assignment.rds")
str(x)
st_crs(x)

# Read in a U.S states shapefile and convert to USLCC
us_states <- read_sf("./cmaq/shapefile_tl2010_us_state_2010/US_state_2010.shp",
                     stringsAsFactors = F) %>%
  filter(!(NAME10 %in% c("Alaska", "Hawaii", "Puerto Rico")))
us_states_uslcc <- st_transform(us_states, crs = st_crs(x))

cv_list <- readRDS("./pm25/blockcv_overall.rds")
df <- cv_list$df
sbcv <- cv_list$sbcv

df$Fold <- fct_inseq(as_factor(df$fold_block))
df_sf <- left_join(x, df, by = c("cmaq_id" = "cmaq_id"))

# Create colors for folds
my_brewer <- brewer.pal(n = 10L, name = "Set3")
names(my_brewer) <- levels(df_sf$Fold)

p_sf <- ggplot() +
  geom_sf(data = st_geometry(us_states_uslcc), color = "black", 
          size = 0.1, alpha = 0.2) +
  geom_sf(data = df_sf, 
          aes_string(color = "Fold"), 
          shape = 16, size = 1, alpha = 0.75) + 
  # To plot projected rather than lat/lon
  coord_sf(datum = NULL) + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        #legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        panel.background = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 5)) + 
  ggtitle("Constant spatial cross-validation folds") +
  scale_colour_manual(name = "Fold", values = my_brewer,
                      na.translate = F) +
  guides(colour = guide_legend(keywidth = 1, keyheight = 1, 
                               override.aes = list(size=2)))
  #guides(colour = guide_legend(override.aes = list(size=1.5)))

ggsave(paste0("./pm25/cvplots/blockcv_overall_states.jpeg"),
       p_sf, width = 5, height = 3.5)
ggsave(paste0("./pm25/cvplots/blockcv_overall_states.png"), 
       p_sf, width = 5, height = 3.5)

rm(sbcv, p_sf, df_sf, df, cv_list, my_brewer)

#### Daily plots ####
fold_list <- readRDS("./pm25/daily_cv.rds")
days = 182:212
## Constant spatial CV fold ##
for (day in days) {
  print(day)
  day_name <- paste0("day", day)
  df_sf <- left_join(x, 
                     fold_list$df_list[[day_name]] %>% 
                       select(cmaq_x, cmaq_y, cmaq_id, fold_b1, pm25_value),
                     by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(Fold = fct_inseq(as_factor(fold_b1)))
  
  # Create colors for folds
  my_brewer <- brewer.pal(n = 10L, name = "Set3")
  names(my_brewer) <- levels(df_sf$Fold)
  
  p_sf <- ggplot() +
    geom_sf(data = st_geometry(us_states_uslcc), color = "black", 
            size = 0.1, alpha = 0.2) +
    geom_sf(data = df_sf, 
            aes_string(color = "Fold"), 
            shape = 16, size = 1, alpha = 0.75) + 
    # To plot projected rather than lat/lon
    coord_sf(datum = NULL) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          #legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          panel.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 5)) + 
    ggtitle(paste0("Constant spatial CV: July ", day - 182 + 1, ", 2011")) +
    scale_colour_manual(name = "Fold", values = my_brewer,
                        na.translate = F) +
    guides(colour = guide_legend(keywidth = 1, keyheight = 1, 
                                 override.aes = list(size=2)))
  
  ggsave(paste0("./pm25/cvplots/b1_jpeg/b1_", day, ".jpeg"),
         p_sf, width = 5, height = 3.5, dpi = 200)
  rm(p_sf, df_sf, my_brewer)
}

## Varying spatial CV ##
days = 182:212
# Varying spatially clustered cross-validation: b2
for (day in days) {
  print(day)
  day_name <- paste0("day", day)
  df_sf <- left_join(x, 
                  fold_list$df_list[[day_name]] %>% 
                    select(cmaq_x, cmaq_y, cmaq_id, fold_b2, pm25_value),
                  by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(Fold = fct_inseq(as_factor(fold_b2)))

  # Create colors for folds
  my_brewer <- brewer.pal(n = 10L, name = "Set3")
  names(my_brewer) <- levels(df_sf$Fold)
  
  p_sf <- ggplot() +
    geom_sf(data = st_geometry(us_states_uslcc), color = "black", 
            size = 0.1, alpha = 0.2) +
    geom_sf(data = df_sf, 
            aes_string(color = "Fold"), 
            shape = 16, size = 1, alpha = 0.75) + 
    # To plot projected rather than lat/lon
    coord_sf(datum = NULL) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          #legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          panel.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 5)) + 
    ggtitle(paste0("Varying spatial CV: July ", day - 182 + 1, ", 2011")) +
    scale_colour_manual(name = "Fold", values = my_brewer,
                        na.translate = F) +
    guides(colour = guide_legend(keywidth = 1, keyheight = 1, 
                                 override.aes = list(size=2)))
  
  ggsave(paste0("./pm25/cvplots/b2_jpeg/b2_", day, ".jpeg"),
         p_sf, width = 5, height = 3.5, dpi = 200)
  rm(p_sf, df_sf, my_brewer)
}

## Repeat for random (non-spatial) CV: fold_cv (not b1 or b2) ##
days = 182:212
for (day in days) {
  print(day)
  day_name <- paste0("day", day)
  df_sf <- left_join(x, 
                     fold_list$df_list[[day_name]] %>% 
                       select(cmaq_x, cmaq_y, cmaq_id, fold_cv, pm25_value),
                     by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(Fold = fct_inseq(as_factor(fold_cv)))
  
  # Create colors for folds
  my_brewer <- brewer.pal(n = 10L, name = "Set3")
  names(my_brewer) <- levels(df_sf$Fold)
  
  p_sf <- ggplot() +
    geom_sf(data = st_geometry(us_states_uslcc), color = "black", 
            size = 0.1, alpha = 0.2) +
    geom_sf(data = df_sf, 
            aes_string(color = "Fold"), 
            shape = 16, size = 1, alpha = 0.75) + 
    # To plot projected rather than lat/lon
    coord_sf(datum = NULL) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          #legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          panel.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 5)) + 
    ggtitle(paste0("Random CV: July ", day - 182 + 1, ", 2011")) +
    scale_colour_manual(name = "Fold", values = my_brewer,
                        na.translate = F) +
    guides(colour = guide_legend(keywidth = 1, keyheight = 1, 
                                 override.aes = list(size=2)))
  
  ggsave(paste0("./pm25/cvplots/ran_jpeg/ran_", day, ".jpeg"),
         p_sf, width = 5, height = 3.5, dpi = 200)
  rm(p_sf, df_sf, my_brewer)
}
