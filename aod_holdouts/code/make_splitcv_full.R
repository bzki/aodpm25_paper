# Break the full observed AOD data into CV folds

# 2 methods I consider for creating spatially clustered CVs
# 1 -- blockCV through spatialBlock method
# 2 -- kmeans with, say, 30 to 40 blocks, and then reducing to 10 blocks.
# Both work fine -- blockCV appears to generally split things up more
# Both (as I have constructed) are random and require seed setting.
# Use blockCV package

## Function does as follows
# Read full AOD data
# Set seed equal to day
# Use blockCV::spatialBlock to construct daily CV
# Create plot
# Save training data with fold information -- include nearest-neighbor distance
# Save data split into training/test as list

library(dplyr) # Easier data.frame manipulation
library(blockCV) # For generating folds
library(sp)
library(RColorBrewer) # For coloring plots
library(fields) # For calculating distance to nearest training point
library(sessioninfo)

sessionInfo()
session_info()

make_splitcv_full <- function(days, save.data = TRUE,
                              theRange = 400000, cv.k = 10,
                              prefix_name = "split_fcv_",
                              png_width = 800, png_height = 600) {
  
  # full range of xy in plots
  full_xy <- readRDS("./rawrds/n4_day182.rds") %>%
    dplyr::select(cmaq_x, cmaq_y)
  xlim_range <- range(full_xy$cmaq_x)
  ylim_range <- range(full_xy$cmaq_y)
  rm(full_xy)
  
  # Generate colors for CV plots
  my_brewer <- brewer.pal(n = cv.k, name = "Set3")
  
  # Ensure that the plot directory exists
  if (!dir.exists("./aod_holdouts/cvplots_full/")) {
    print("Making directory for cvplots")
    dir.create("./aod_holdouts/cvplots_full/")
  }
  
  for (day in days) {
    
    print(day)  
    
    day_name = paste0("day", day)
    # Read full observed data
    x <- readRDS(paste0("./rawrds/aod_", day, ".rds"))
    testthat::expect_identical(x$day[1], day)
    
    # Cross-validation folds based only on non-missing data
    x_nontest <- x
    xy_nontest <- x_nontest[, c("cmaq_x", "cmaq_y")]
    rm(x)
    
    # Use the day as the seed to construct validation holdout
    xy_sp <- xy_nontest
    coordinates(xy_sp) <- ~cmaq_x + cmaq_y
    set.seed(day) # Use the day as seed for reproducibility
    sb <- spatialBlock(xy_sp,
                       theRange = theRange,
                       species = NULL, 
                       k = cv.k, selection = "random",
                       verbose = FALSE, progress = FALSE, showBlocks = FALSE)
    warnings()
    print(warnings())
    rm(xy_sp)
    
    # Create plot and save -- use same ylim/xlim for all plots
    png(paste0("./aod_holdouts/cvplots_full/", prefix_name, day, ".png"),
        width = png_width, height = png_height)
    plot(xy_nontest[, c("cmaq_x", "cmaq_y")], cex = 0.2, type = "n",
         xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         xlim = xlim_range, ylim = ylim_range,
         main = paste0("CV Folds for day ", day))
    points(xy_nontest[, c("cmaq_x", "cmaq_y")],
           col = my_brewer[sb$foldID], cex = 0.4)
    dev.off()
    
    x_nontest$foldID <- sb$foldID
    
    # Construct distance between validation and non-validation for each fold
    x_nontest$cv_nndist <- NA
    for (fi in seq(cv.k)) {
      nn_dist <- 
        apply(rdist(x_nontest[sb$folds[[fi]][[2]], c("cmaq_x", "cmaq_y")],
                    x_nontest[sb$folds[[fi]][[1]], c("cmaq_x", "cmaq_y")]),
              1, min)
      x_nontest$cv_nndist[sb$folds[[fi]][[2]]] <- nn_dist
    }
    
    # Output
    res_list <- list(
      df = x_nontest,
      folds = sb$folds,
      records = sb$records,
      range = sb$range
    )
    if (save.data) {
      saveRDS(res_list, paste0("./rawrds/", prefix_name,  day, ".rds"))
    }
    rm(sb, x_nontest, xy_nontest, fi, nn_dist, res_list)
    gc()
  }
  warnings()
  print(warnings())
}

# 
# debugonce(make_splitcv_full)
# make_splitcv_full(days = 182L, save.data = F,
#                   theRange = 400000, cv.k = 10,
#                   prefix_name = "split_fcv_")

make_splitcv_full(days = 182:212, save.data = TRUE,
                  theRange = 400000, cv.k = 10,
                  prefix_name = "split_fcv_")
