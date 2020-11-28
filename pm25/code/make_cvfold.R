# 3 forms of cross-validation creation for daily PM2.5 data
library(dplyr)
library(ggplot2)
library(forcats)
library(sp)
library(blockCV)
sessionInfo()

if (!dir.exists("./pm25/")) {
  dir.create("./pm25/")
  dir.create("./pm25/cvplots/")
}

# 10-fold cross-validation
# -- Loop through each day
# -- Remove narr_ variables, lat/lon
# -- Limit to observed PM2.5 only 
# -- Use shuffle and n mod 10 technique

days = 182:212
df_list <- vector("list", length(days))
names(df_list) <- paste0("day", days)

for (day in days) {
  
  print(day)
  x <- readRDS(paste0("./rawrds/n4_day", day, ".rds"))
  testthat::expect_identical(day, x$day[1])
  x_range <- range(x$cmaq_x)
  y_range <- range(x$cmaq_y)
  
  df <- x %>%
    filter(!is.na(pm25_value)) %>%
    dplyr::select(-starts_with("narr"), -lon, -lat)

  # Key steps -- shuffle data, then use modulo 10 to create folds
  # https://stats.stackexchange.com/a/48161/208322
  set.seed(day)
  df2 <- df[sample.int(nrow(df), replace = F), ]
  testthat::expect_setequal(df2$cmaq_id, df$cmaq_id)
  mod_10 <- seq(nrow(df2)) %% 10L 
  mod_10 <- replace(mod_10, which(mod_10 == 0L), 10L)

  df2$fold <- as_factor(as.character(mod_10))
  df2$aod_missing <- fct_relevel(
    as_factor(
      ifelse(is.na(df2$aod_value), "Missing", "Observed")
    ), "Observed", "Missing")
  png(paste0("./pm25/cvplots/fold1_", day, ".png"),
      width = 800, height = 600)
  print(
    ggplot(data=df2) +
    geom_point(aes(x = cmaq_x, y = cmaq_y, color = fold),
               size = 3, alpha = 0.75) + 
    xlim(x_range) +
    ylim(y_range) + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          #legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          #panel.background = element_blank(),
          legend.text = element_text(size = 12)) + 
    ggtitle(paste0("PM2.5 CV folds on day ", day)) +
    scale_colour_brewer(palette = "Set3", name = "Fold") +
    coord_fixed()
  )
  dev.off()  
  png(paste0("./pm25/cvplots/fold1m_", day, ".png"),
                 width = 800, height = 600)
  print(
    ggplot(data=df2) +
      geom_point(aes(x = cmaq_x, y = cmaq_y, color = fold, shape = aod_missing),
                 size = 3, alpha = 0.75) + 
      xlim(x_range) +
      ylim(y_range) + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks = element_blank(),
            #legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            #panel.background = element_blank(),
            legend.text = element_text(size = 12)) + 
      ggtitle(paste0("PM2.5 CV folds on day ", day)) +
      scale_colour_brewer(palette = "Set3", name = "Fold") + 
      scale_shape_manual(values = c(16, 17), name = "AOD") +
      guides(shape = guide_legend(order = 2),col = guide_legend(order = 1)) +
      coord_fixed()
  )
  dev.off()
  
  # Create data.frame to save
  df_fold <- df2
  df_fold$fold_cv <- mod_10
  day_name <- paste0("day", day)
  df_list[[day_name]] <- df_fold
  rm(df_fold, df2, df, x, mod_10)
}

# 10-fold spatial cross-validation
# Use spatial blocking (random) using blockCVs
# The blocks can be pre-defined (b1) or defined daily (b2)
# The other main parameter is the spatial block width
# I propose setting this to 150km

## Tasks
# - Common spatial blocks across all days + plot
# - Plot for each day using the common spatial blocks
# - Distinct blocks for each day + plot
# - Save data.frame with fold IDs as well as list of folds by index
rm(list=setdiff(ls(), c("df_list")))

# Get x/y range for plot range
tmp_df <- readRDS(paste0("./rawrds/n4_day182.rds"))
x_range <- range(tmp_df$cmaq_x)
y_range <- range(tmp_df$cmaq_y)
rm(tmp_df)

days = 182:212
df_combine <- NULL
for (day in days) {
  x <- df_list[[paste0("day", day)]]
  testthat::expect_identical(day, x$day[1])
  df <- x %>% dplyr::select(cmaq_id, cmaq_x, cmaq_y)
  df_combine <- rbind(df_combine, df)
}
rm(df, x, day)
df_sum <- df_combine %>% 
  group_by(cmaq_id) %>%
  summarize(
    min_x = min(cmaq_x),
    max_x = max(cmaq_x),
    min_y = min(cmaq_y),
    max_y = max(cmaq_y)
  ) %>%
  ungroup(.)
testthat::expect_equivalent(df_sum$min_x, df_sum$max_x)
testthat::expect_equivalent(df_sum$min_y, df_sum$max_y)

df_sum2 <- df_sum %>%
  select(-max_x, -max_y) %>%
  rename(cmaq_x = min_x, cmaq_y = min_y) %>%
  arrange(cmaq_id) %>%
  as.data.frame(.)
nrow(df_sum2) # 1034 total locations in the U.S.

df_sp <- df_sum2
coordinates(df_sp) <- ~cmaq_x + cmaq_y
set.seed(2)
sbcv <- 
  spatialBlock(df_sp, theRange = 150000, k = 10L, selection = "random",
               showBlocks = T, progress = T, verbose = T)
df_sum2$fold_block <- sbcv$foldID
saveRDS(list(df=df_sum2, sbcv = sbcv), "./pm25/blockcv_overall.rds")

png("./pm25/cvplots/blockcv_overall.png",
    width = 800, height = 600)
print(
ggplot(data=df_sum2) +
  geom_point(aes(x = cmaq_x, y = cmaq_y, 
                 color = fct_inseq(as_factor(fold_block))),
             size = 3, alpha = 0.75) + 
  xlim(x_range) +
  ylim(y_range) + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        #legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #panel.background = element_blank(),
        legend.text = element_text(size = 12)) + 
  ggtitle(paste0("PM2.5 spatial block folds -- overall")) +
  scale_colour_brewer(palette = "Set3", name = "Fold") +
  coord_fixed()
)
dev.off()
rm(df_sp, df_combine, df_sum, sbcv)

b1_folds <- vector("list", length(days))
names(b1_folds) <- paste0("day", days)
cv_folds = b2_folds = b1_folds

# Loop through -- create two additional kinds of block CV folds
# -- fold_b1 using the overall spatial block assignment above
# -- fold_b2 using a daily random creation of spatial blocks

for (day in days) {
  print(day)
  day_name <- paste0("day", day)
  
  # Joint with the overall block assignment
  df <- df_list[[day_name]] %>%
    left_join(df_sum2 %>% dplyr::select(-cmaq_x, -cmaq_y),
              by = c("cmaq_id" = "cmaq_id")) %>%
    rename(fold_b1 = fold_block) %>%
    select(-fold)
 
  # Plot and save
  png(paste0("./pm25/cvplots/b1_", day, ".png"),
      width = 800, height = 600)
  print(
    ggplot(data=df) +
      geom_point(aes(x = cmaq_x, y = cmaq_y, 
                     color = fct_inseq(as_factor(fold_b1))),
                 size = 3, alpha = 0.75) + 
      xlim(x_range) +
      ylim(y_range) + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks = element_blank(),
            #legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            #panel.background = element_blank(),
            legend.text = element_text(size = 12)) + 
      ggtitle(paste0("PM2.5 overall spatial blocking for day ", day)) +
      scale_colour_brewer(palette = "Set3", name = "Fold") +
      coord_fixed()
  )
  dev.off()
  
  png(paste0("./pm25/cvplots/b1m_", day, ".png"),
      width = 800, height = 600)
  print(
    ggplot(data=df) +
      geom_point(aes(x = cmaq_x, y = cmaq_y, 
                     color = fct_inseq(as_factor(fold_b1)),
                     shape = aod_missing),
                 size = 3, alpha = 0.75) + 
      xlim(x_range) +
      ylim(y_range) + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks = element_blank(),
            #legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            #panel.background = element_blank(),
            legend.text = element_text(size = 12)) + 
      ggtitle(paste0("PM2.5 overall spatial blocking for day ", day)) +
      scale_colour_brewer(palette = "Set3", name = "Fold") +
      scale_shape_manual(values = c(16, 17), name = "AOD") +
      guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
      coord_fixed()
  )
  dev.off()
  
  # Use day as seed and create random spatial blocking for each day (b2)
  # -- range of blocks is 150km in width
  df_sp <- df
  coordinates(df_sp) <- ~cmaq_x + cmaq_y
  set.seed(day)
  suppressWarnings(
    sb2 <- 
      spatialBlock(df_sp, theRange = 150000, k = 10L, selection = "random",
                   showBlocks = F, progress = F, verbose = F)
  )
  # Warnings are just due to no CRS defined -- correctly assumes project
  df$fold_b2 <- sb2$foldID
  df_list[[day_name]] <- df # Re-assign the new df to the list
  
  png(paste0("./pm25/cvplots/b2_", day, ".png"),
      width = 800, height = 600)
  print(
    ggplot(data=df) +
      geom_point(aes(x = cmaq_x, y = cmaq_y, 
                     color = fct_inseq(as_factor(fold_b2))),
                 size = 3, alpha = 0.75) + 
      xlim(x_range) +
      ylim(y_range) + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks = element_blank(),
            #legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            #panel.background = element_blank(),
            legend.text = element_text(size = 12)) + 
      ggtitle(paste0("PM2.5 daily spatial blocking for day ", day)) +
      scale_colour_brewer(palette = "Set3", name = "Fold") +
      coord_fixed()
  )
  dev.off()
  
  png(paste0("./pm25/cvplots/b2m_", day, ".png"),
      width = 800, height = 600)
  print(
    ggplot(data=df) +
      geom_point(aes(x = cmaq_x, y = cmaq_y, 
                     color = fct_inseq(as_factor(fold_b2)),
                     shape = aod_missing),
                 size = 3, alpha = 0.75) + 
      xlim(x_range) +
      ylim(y_range) + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks = element_blank(),
            #legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            #panel.background = element_blank(),
            legend.text = element_text(size = 12)) + 
      ggtitle(paste0("PM2.5 daily spatial blocking for day ", day)) +
      scale_colour_brewer(palette = "Set3", name = "Fold") +
      scale_shape_manual(values = c(16, 17), name = "AOD") +
      guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
      coord_fixed()
  )
  dev.off()
  
  # Create lists of folds for each of the 10 folds
  b1_list = b2_list = cv_list = vector("list", 10L)
  for (fi in seq(length(cv_list))) {
    cv_list[[fi]][[1]] <- which(df$fold_cv != fi)
    cv_list[[fi]][[2]] <- which(df$fold_cv == fi)
    
    b1_list[[fi]][[1]] <- which(df$fold_b1 != fi)
    b1_list[[fi]][[2]] <- which(df$fold_b1 == fi)
    
    b2_list[[fi]][[1]] <- which(df$fold_b2 != fi)
    b2_list[[fi]][[2]] <- which(df$fold_b2 == fi)
  }
  cv_folds[[day_name]] <- cv_list
  b1_folds[[day_name]] <- b1_list
  b2_folds[[day_name]] <- b2_list
  
  rm(df, sb2, df_sp, cv_list, b1_list, b2_list, day_name)
}


saveRDS(list(df_list = df_list, 
             cv_folds = cv_folds, b1_folds = b1_folds, b2_folds = b2_folds), 
        "./pm25/daily_cv.rds")

