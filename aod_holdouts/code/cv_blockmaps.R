# Create a nicer plot of the spatial block cross-validation folds
library(sf)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
sessionInfo()

# Create colors for folds
my_brewer <- brewer.pal(n = 10L, name = "Set3")

# Read in previously created CMAQ sf file 
x <- readRDS("cmaq/region_assignment.rds")
str(x)
st_crs(x)

# Read in a U.S states shapefile and convert to USLCC
us_states <- read_sf("./cmaq/shapefile_tl2010_us_state_2010/US_state_2010.shp",
                     stringsAsFactors = F) %>%
  filter(!(NAME10 %in% c("Alaska", "Hawaii", "Puerto Rico")))
us_states_uslcc <- st_transform(us_states, crs = st_crs(x))

days = 182:212L

for (day in days) {
  
  print(day)  
  daily_list <- readRDS(paste0("./rawrds/splitcv_4_", day, ".rds"))
  df_nontest <- left_join(x, 
                          daily_list$train$df_train %>% 
                            select(cmaq_id, foldID, cmaq_x, cmaq_y),
                          by = c("cmaq_id" = "cmaq_id")) %>%
    mutate(fold = as.factor(foldID))
  
  names(my_brewer) <- levels(df_nontest$fold)
  
  p_sf <- ggplot() +
      geom_sf(data = df_nontest, 
              aes_string(color = "fold"), 
              shape = 15, size = 0.2) + 
      geom_sf(data = st_geometry(us_states_uslcc), color = "black", 
              size = 0.2, alpha = 0) +
      # To plot projected rather than lat/lon
      coord_sf(datum = NULL) + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks = element_blank(),
            #legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
            panel.background = element_blank()) + 
      ggtitle(paste0("July ", day - 182 + 1, ", 2011: Spatial cross-validation folds")) +
      scale_colour_manual(name = "Fold", values = my_brewer,
                          na.translate = F) +
    guides(colour = guide_legend(override.aes = list(size=2)))
    
    ggsave(paste0("./aod_holdouts/paper_plots/cvplots/blockcv_", day, ".jpeg"), 
           p_sf, width = 5, height = 3.5, dpi = 200)
    rm(daily_list, p_sf, df_nontest)
    gc()
}


