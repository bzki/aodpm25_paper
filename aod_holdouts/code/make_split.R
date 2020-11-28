# Create easy to read and mess with data files based on train/validation/test split
# Single data.frame with indicator variable these 3 categories

source("./functions/hole.R", echo = T, max.deparse.length = Inf)
library(dplyr)
library(ranger)
library(sessioninfo)
sessionInfo()
session_info()

make_split <- function(days, holdout_list, val_radius = 175, num.holes = 10,
                       prefix_name = "./rawrds/split_si_") {

  for (day in days) {
    
    print(day)  
    
    day_name = paste0("day", day)
    x <- readRDS(paste0("./rawrds/aod_", day, ".rds"))
    testthat::expect_identical(x$day[1], day)
    
    # Exclude holdouts for parameter tuning
    test_index <- holdout_list[[day_name]]$holdout_index
    x_test <- x[test_index, ]
    x_nontest <- x[-test_index, ]
    xy_nontest <- x_nontest[, c("cmaq_x", "cmaq_y")]
    
    # Use the day as the seed to construct validation holdout
    set.seed(day)
    circle_info <- getHoles.fast(xy_nontest, num.holes = num.holes, 
                                 dist.hole = val_radius * 1000)
    val_holdout <- circle_info$omitted.index
    x_train <- x_nontest[-val_holdout, ]
    x_val <- x_nontest[val_holdout, ]
    
    # Save a data.frame with indicator for train/val/test split
    x_test$split <- "Test"
    x_train$split <- "Train"
    x_val$split <- "Val"
    
    x_comb <- rbind(x_test, x_train, x_val)
    testthat::expect_true(base::setequal(x_comb$cmaq_id, x$cmaq_id))
    saveRDS(x_comb, paste0(prefix_name,  day, ".rds"))
    rm(x_comb, x_nontest, x_test, x_train, x_val, xy_nontest, x, day_name,
       test_index, val_holdout)
  }
}

holdout_1 <- readRDS("./aod_holdouts/holdout_1_random10.rds")
holdout_2 <- readRDS("./aod_holdouts/holdout_2_random60.rds")
holdout_3 <- readRDS("./aod_holdouts/holdout_3_circleA.rds")
holdout_4 <- readRDS("./aod_holdouts/holdout_4_circleB.rds")

days = 182:212

make_split(days = days, holdout_list = holdout_1, 
           val_radius = 175, num.holes = 10,
           prefix_name = "./rawrds/split_1_")

make_split(days = days, holdout_list = holdout_2, 
           val_radius = 175, num.holes = 10,
           prefix_name = "./rawrds/split_2_")

make_split(days = days, holdout_list = holdout_3, 
           val_radius = 175, num.holes = 10,
           prefix_name = "./rawrds/split_3_")

make_split(days = days, holdout_list = holdout_4, 
           val_radius = 175, num.holes = 10,
           prefix_name = "./rawrds/split_4_")

# Plot and ensure it looks reasonable
# x <- readRDS("./rawrds/split_4_212.rds")
# x_train <- x %>% filter(split == "Train") %>% select(cmaq_x, cmaq_y)
# x_test <- x %>% filter(split == "Test") %>% select(cmaq_x, cmaq_y)
# x_val <- x %>% filter(split == "Val") %>% select(cmaq_x, cmaq_y)
# plot(x_train)
# points(x_val, col = "red", cex = 0.2)
# points(x_test, col = "blue", cex = 0.2)
