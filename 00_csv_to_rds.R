# Reads in the raw data in .csv format and outputs into .rds format
# .rds will be much smaller in size
# For the time being, will add .rds to the .gitignore
# NOTE: Currently executing this script via R CMD BATCH

rm(list=ls())
sessioninfo::session_info()
print(getwd())
csv.directory <- "../" # Raw data assumed to be in parent directory
write.directory <- "./rawrds/"
core.name <- "n_aod_pred_2011_new8_"

# Set file name and read in each day, then save as .rds
for (day in 182:212) {
  file.name <- paste0(csv.directory, core.name, day, ".csv")
  print(file.name)
  df <- read.csv(file.name)
  write.name <- paste0(write.directory, core.name, day, ".rds")
  # Saving to version 2 in order for R prior to 3.5 to be able to read
  # See: https://blog.revolutionanalytics.com/2019/05/whats-new-in-r-360.html
  saveRDS(df, write.name, version = 2)
  rm(file.name, df, write.name)
}
