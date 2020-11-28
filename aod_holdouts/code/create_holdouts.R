# Create holdouts for various experiments testing AOD prediction performance
# Holdouts datasets are for each day separately

# Script was run using R 3.6
# NOTE: sample() behaves differently in R 3.6 vs. previous versions
# See: https://blog.revolutionanalytics.com/2019/05/whats-new-in-r-360.html

# For each experiment, save a list of data.frames for each day with two columns
# - index of testing observations
# - cmaq_id for testing observations

# Testing data created in the following ways:
# 1: Hold out a random 10% of AOD observed values.
# 2: Hold out a random 60% of AOD observed values.
# (NOTE: The mean/median missing # of AOD out of the 53K+ cells is about 60%)
# 3: Hold out 10 random holes excluding points within 125km -- holes A
# 4: Hold out 10 random holes excluding points within 250km -- holes B
# Only the last (#4) is used in subsequent analyses.

source("./functions/hole.R", echo = T, max.deparse.length = Inf)

# Create 2 directories if they do not already exist.
if (!dir.exists("aod_holdouts")) {
  dir.create("aod_holdouts")
}

if (!dir.exists("aod_holdouts/plots")) {
  dir.create("aod_holdouts/plots")
}

# May wish to detect the R version info to ensure it is 3.6.0 or above
# otherwise the sample() function may behave unexpectedly
sessioninfo::session_info()
sessioninfo::platform_info()
R.Version()$major
R.Version()$minor
# This information is relevant to the sample() function
print(RNGkind())
# R 3.6 and above should be: 
# [1] "Mersenne-Twister" "Inversion"        "Rejection" 
# R 3.5: 
# [1] "Mersenne-Twister" "Inversion"

days = 182:212

## 1: 10% Random Holdout
# Create holdout list: each entry contains the holdout for the specified day
random10_holdout = vector(mode = "list", length = length(days))
names(random10_holdout) <- paste0("day", days)
# Also save additional information -- seed, R version, platform, description
random10_holdout$seed <- 12221
random10_holdout$description <- "Random 10% holdout based on aod_XXX.rds"
random10_holdout$R.version <- R.version.string
random10_holdout$platform <- R.version$platform

set.seed(random10_holdout$seed)
for (day in days) {
  read.name <- paste0("./rawrds/aod_", day, ".rds")
  print(read.name)
  x <- readRDS(read.name)
  holdout <- sample.int(n = nrow(x), size = floor(0.10 * nrow(x)))
  day_df <- data.frame(day = day, holdout_index = holdout, 
                       holdout_id = x$cmaq_id[holdout])
  list_entry <- paste0("day", day)
  random10_holdout[[list_entry]] <- day_df
  rm(read.name, x, holdout, day_df, list_entry)
}
saveRDS(random10_holdout, "./aod_holdouts/holdout_1_random10.rds", version = 2)
# rm(random10_holdout, day)
rm(day)

# 2: 60% Random Holdout
# Create holdout list: each entry contains the holdout for the specified day
random60_holdout = vector(mode = "list", length = length(days))
names(random60_holdout) <- paste0("day", days)
# Also save additional information -- seed, R version, platform, description
random60_holdout$seed <- 12222
random60_holdout$description <- "Random 60% holdout based on aod_XXX.rds"
random60_holdout$R.version <- R.version.string
random60_holdout$platform <- R.version$platform

set.seed(random60_holdout$seed)
for (day in days) {
  read.name <- paste0("./rawrds/aod_", day, ".rds")
  print(read.name)
  x <- readRDS(read.name)
  holdout <- sample.int(n = nrow(x), size = floor(0.60 * nrow(x)))
  day_df <- data.frame(day = day, holdout_index = holdout, 
                       holdout_id = x$cmaq_id[holdout])
  list_entry <- paste0("day", day)
  random60_holdout[[list_entry]] <- day_df
  rm(read.name, x, holdout, day_df, list_entry)
}
saveRDS(random60_holdout, "./aod_holdouts/holdout_2_random60.rds", version = 2)
# rm(random60_holdout, day)
rm(day)


# 3: Holdout 10 random large circular areas with radius 125km
# Create holdout list: each entry contains the holdout for the specified day
circleA_holdout = vector(mode = "list", length = length(days))
names(circleA_holdout) <- paste0("day", days)
# Also save additional information -- seed, R version, platform, description
circleA_holdout$seed <- 12223
circleA_holdout$description <- paste0("10 random circles with radius 125km",
                                      " based on aod_XXX.rds")
circleA_holdout$R.version <- R.version.string
circleA_holdout$platform <- R.version$platform

set.seed(circleA_holdout$seed)
for (day in days) {
  read.name <- paste0("./rawrds/aod_", day, ".rds")
  print(read.name)
  x <- readRDS(read.name)
  xy <- x[, c("cmaq_x", "cmaq_y")]
  circle_info <- getHoles.fast(xy, num.holes = 10, dist.hole = 125 * 1000)
  holdout <- circle_info$omitted.index
  holdout_portion = length(holdout) / nrow(xy) * 100
  print(paste0(length(holdout), ", or ", round(holdout_portion, 2), "% of ",
               nrow(xy), " total observations are held out"))
  day_df <- data.frame(day = day, holdout_index = holdout, 
                       holdout_id = x$cmaq_id[holdout])
  list_entry <- paste0("day", day)
  circleA_holdout[[list_entry]] <- day_df
  rm(read.name, x, xy, holdout, day_df, list_entry, circle_info,
     holdout_portion)
}
saveRDS(circleA_holdout, "./aod_holdouts/holdout_3_circleA.rds", version = 2)
# rm(circleA_holdout, day)
rm(day)


# 4: Holdout 10 random large circular areas with radius 250km
# Create holdout list: each entry contains the holdout for the specified day
circleB_holdout = vector(mode = "list", length = length(days))
names(circleB_holdout) <- paste0("day", days)
# Also save additional information -- seed, R version, platform, description
circleB_holdout$seed <- 12224
circleB_holdout$description <- paste0("10 random circles with radius 250km",
                                      " based on aod_XXX.rds")
circleB_holdout$R.version <- R.version.string
circleB_holdout$platform <- R.version$platform

set.seed(circleB_holdout$seed)
for (day in days) {
  read.name <- paste0("./rawrds/aod_", day, ".rds")
  print(read.name)
  x <- readRDS(read.name)
  xy <- x[, c("cmaq_x", "cmaq_y")]
  circle_info <- getHoles.fast(xy, num.holes = 10, dist.hole = 250 * 1000)
  holdout <- circle_info$omitted.index
  holdout_portion = length(holdout) / nrow(xy) * 100
  print(paste0(length(holdout), ", or ", round(holdout_portion, 2), "% of ",
               nrow(xy), " total observations are held out"))
  day_df <- data.frame(day = day, holdout_index = holdout, 
                       holdout_id = x$cmaq_id[holdout])
  list_entry <- paste0("day", day)
  circleB_holdout[[list_entry]] <- day_df
  rm(read.name, x, xy, holdout, day_df, list_entry, circle_info,
     holdout_portion)
}
saveRDS(circleB_holdout, "./aod_holdouts/holdout_4_circleB.rds", version = 2)
# rm(circleB_holdout, day)
rm(day)

##########################################################
# - Create plots for each day under each holdout scenario
# - Use package magick to construct gifs
# - Save all to aod_holdouts/plots/ sub-folders

library(magick) # Note -- may require additional installation depending on OS
library(dplyr)

plot1_names <- paste0("./aod_holdouts/plots/holdout_1_", days,".png")
plot2_names <- paste0("./aod_holdouts/plots/holdout_2_", days,".png")
plot3_names <- paste0("./aod_holdouts/plots/holdout_3_", days,".png")
plot4_names <- paste0("./aod_holdouts/plots/holdout_4_", days,".png")

# Load in the a full day and save the full x/y range so that
#  the x/y lims are the same for all plots
df_names <- paste0("./rawrds/aod_", days, ".rds")
full.df <- readRDS("./rawrds/n4_day182.rds")
x.range = range(full.df$cmaq_x)
y.range = range(full.df$cmaq_y)
rm(full.df)

# Create lists for each holdout -- these will contain the images
maps_1 = maps_2 = maps_3 = maps_4 <- list()

# Loop through, create .png, and save in list
for (day in days) {
  counter = day - 181 # Create counter to index list starting from 1
  
  # Read in data for day
  readname = df_names[counter]
  print(paste0("Read in: ", readname))
  x <- readRDS(readname) %>% 
    dplyr::select(cmaq_x, cmaq_y)
  
  # Get holdout indices for the day
  day_name = paste0("day", day)
  holdout_current <- list(random10_holdout[[day_name]]$holdout_index,
                          random60_holdout[[day_name]]$holdout_index,
                          circleA_holdout[[day_name]]$holdout_index,
                          circleB_holdout[[day_name]]$holdout_index)

  savenames_current = list(plot1_names[counter],
                           plot2_names[counter],
                           plot3_names[counter],
                           plot4_names[counter])

  print(paste0("Day ", day))
  for (pi in seq(length(holdout_current))) {

    # Create map and save
    print(savenames_current[[pi]])
    png(filename = savenames_current[[pi]])
    # png(filename=savenames_current[[pi]], width=577, height=433)
    
    holdout_pi = holdout_current[[pi]]
    x_in = x[-holdout_pi, ]
    x_out = x[holdout_pi, ]
    
    plot(x_in, 
         xlab="", ylab="", xaxt = "n", yaxt = "n", bty="n",
         main = paste("Day", day),
         pch = 16, col = "red", cex=0.15, 
         xlim=x.range,
         ylim=y.range)
    points(x_out, 
           pch = 1, col = "blue", cex = 0.3)
    legend("bottomleft", c("Observed", "Holdout"), 
           pch=c(16, 1), col=c("red", "blue"))
    dev.off()
    
    rm(x_in, x_out, holdout_pi)
  }
  
  maps_1[[counter]] <- magick::image_read(path=plot1_names[counter])    
  maps_2[[counter]] <- magick::image_read(path=plot2_names[counter])    
  maps_3[[counter]] <- magick::image_read(path=plot3_names[counter])    
  maps_4[[counter]] <- magick::image_read(path=plot4_names[counter])
  
  rm(x, holdout_current, savenames_current, day_name, readname)
}


# Create gifs using image_animate in magick package
temp_1 <- image_animate(image_join( maps_1 ), fps=2, loop=0)
image_write(temp_1, path="./aod_holdouts/plots/holdout_1_all.gif", format="gif")

temp_2 <- image_animate(image_join( maps_2 ), fps=2, loop=0)
image_write(temp_2, path="./aod_holdouts/plots/holdout_2_all.gif", format="gif")

temp_3 <- image_animate(image_join( maps_3 ), fps=2, loop=0)
image_write(temp_3, path="./aod_holdouts/plots/holdout_3_all.gif", format="gif")

temp_4 <- image_animate(image_join( maps_4 ), fps=2, loop=0)
image_write(temp_4, path="./aod_holdouts/plots/holdout_4_all.gif", format="gif")

