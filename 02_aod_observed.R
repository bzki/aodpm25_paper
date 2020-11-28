# Read in the .rds files that have NLDAS missing values filled in
# (n4 versions -- see 01_fill_nldas.R)
# Output datasets only with the rows where AOD is observed

rm(list=ls())
library(dplyr)
sessioninfo::session_info()
print(getwd())
file.directory <- "./rawrds/"
read.name <- "n4_day"

# Placeholder matrix for recording percent missing
aodobs_matrix <- matrix(NA, nrow = length(182:212), ncol = 4)
class(aodobs_matrix) <- "numeric"
colnames(aodobs_matrix) <- c("day", "total", "observed", "pct_observed")
aodobs_matrix[, 1] <- 182:212

# Set file name and read in each day, then save as .rds
for (day in 182:212) {
  file.name <- paste0(file.directory, read.name, day, ".rds")
  print(file.name)
  df <- readRDS(file.name)
  df2 <- df %>%
    filter(!is.na(aod_value))
  pct_observed = nrow(df2) / nrow(df) * 100
  counter = day - 181
  aodobs_matrix[counter, "total"] <- nrow(df) 
  aodobs_matrix[counter, "observed"] <- nrow(df2) 
  aodobs_matrix[counter, "pct_observed"] <- pct_observed
  print(paste0("Day ", day, ": ", 
               nrow(df2), "(", round(pct_observed, 2), "%)",
               " observations with AOD observed out of ",
               nrow(df), " total observations"))
  write.name <- paste0(file.directory, "aod_", day, ".rds")
  # Saving to version 2 in order for R prior to 3.5 to be able to read
  # See: https://blog.revolutionanalytics.com/2019/05/whats-new-in-r-360.html
  saveRDS(df2, write.name, version = 2)
  rm(file.name, df, df2, write.name)
}

# Save missing_matrix
saveRDS(aodobs_matrix, "./rawrds/observed_summary.rds", version = 2)
