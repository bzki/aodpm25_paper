# Read in the rawrds/ data files
# Use the nearest neighbor filling function and fill in NLDAS variables
# Although I do not plan on using them, I also fill in NARR variables
# (This is done to see what we might potentially lose by excluding them)

# For now, output with two methods
# - nearest 1 neighbor
# - nearest 4 neighbors, using weights proportional to 1/(dist^2) 
# Output to rawrds/ with names: n1_dayXXX.rds and n4_dayXXX.rds

source("./functions/nearest.R", echo=TRUE, max.deparse.length=Inf)
library(sessioninfo)
library(testthat)
library(dplyr)
library(stringr)
sessioninfo::session_info()
read_prefix <- "rawrds/n_aod_pred_2011_new8_"

# Helper functions to convert Inf/MAX values to NAs
replace_inf <- function(x) {
  y <- replace(x, is.infinite(x), NA)
  return(y)
}
replace_max <- function(x, max.val) {
  y <- replace(x, x == max.val, NA)
  return(y)
}

# Create the variable list automatically based on the data.frame
get_varnames <- function(df) {
  na.vec <- sapply(df %>% select(starts_with("narr")),
                   function(x) sum(is.na(x)))
  unique.vals <- sort(unique(na.vec))
  # After replacing Inf/Extremes to NAs, should be either 609 or 746 missing
  testthat::expect_identical(unique.vals, c(609L, 746L))
  set.1 <- names(na.vec)[which(na.vec == unique.vals[1])]
  set.2 <- names(na.vec)[which(na.vec == unique.vals[2])]
  return(list(set.1, set.2))
}

tmp.df <- readRDS("./rawrds/n_aod_pred_2011_new8_182.rds")
MAX.VALUE = max(tmp.df$narr_ugrd10m, na.rm=T)
print(MAX.VALUE)
rm(tmp.df)


for (day in 182:212) {
  
  read_name <- paste0(read_prefix, day, ".rds")
  # Read raw RDS file
  print(paste("Read name:", read_name))
  x <- readRDS(read_name)
  
  # Get variable names for nldas
  nldas.names <- x %>% dplyr::select(starts_with("nldas")) %>% names(.)
  
  # Get index of missing values
  na_index <- which(is.na(x[, nldas.names[1]]))
  expect_identical(length(na_index), 3L)
  expect_identical(na_index, c(2305L, 9217L, 29010L))

  # Separate data into two parts
  x_na <- x[na_index, ]
  x_obs <- x[-na_index, ]
  
  # Run function
  nn.nldas1 <- add_nearestN(x_na, x_obs, N.nn = 1, nth = 1, 
                            var.names = nldas.names)
  nn.nldas4 <- add_nearestN(x_na, x_obs, N.nn = 4, nth = 1, 
                            var.names = nldas.names)
  
  print(nn.nldas4 %>% dplyr::select(starts_with("nn_distance")))
  
  # Fill in the missing values with the output of the function
  x_n1 <- x
  x_n1[na_index, nldas.names] <- nn.nldas1[, paste0("nn_", nldas.names)]
  x_n4 <- x
  x_n4[na_index, nldas.names] <- nn.nldas4[, paste0("nn_", nldas.names)]
  
  # Make sure nothing else in the original data has changed
  expect_identical(x_n1[, !(colnames(x_n1) %in% nldas.names)], 
                   x[, !(colnames(x) %in% nldas.names)])
  expect_identical(x_n4[, !(colnames(x_n4) %in% nldas.names)], 
                   x[, !(colnames(x) %in% nldas.names)])
  expect_identical(x_n1[-na_index, nldas.names],
                   x[-na_index, nldas.names])
  expect_identical(x_n4[-na_index, nldas.names],
                   x[-na_index, nldas.names])
  
  rm(x_na, x_obs, nn.nldas1, nn.nldas4, nldas.names, na_index)
  
  # Fill in the NARR variables.
  expect_identical(max(x$narr_ugrd10m, na.rm=T), MAX.VALUE)
  # Replace Inf/Extremes to NAs
  x2 <- x %>% 
    mutate_at(
      vars(starts_with("narr")), 
      replace_max, max.val = MAX.VALUE
    ) %>%
    mutate_at(
      vars(starts_with("narr")),
      replace_inf
    ) 
  
  narr_list <- get_varnames(x2)
  narr1_names <- narr_list[[1]]
  narr2_names <- narr_list[[2]]
  
  narr1_nas <- which(is.na(x2[, narr1_names[1]]))
  
  x2.na1 <- x2[narr1_nas, ]
  x2.obs1 <- x2[-narr1_nas, ]
  
  nn1.narr1 <- add_nearestN(x2.na1, x2.obs1, N.nn = 1, nth = 1, 
                            var.names = narr1_names)
  nn4.narr1 <- add_nearestN(x2.na1, x2.obs1, N.nn = 4, nth = 1, 
                            var.names = narr1_names)
  
  # Test: name order matches
  expect_identical(
    str_sub(names(nn1.narr1), 4)[1:length(narr1_names)],
    narr1_names
  )
  
  # Replace missing values
  x_n1[narr1_nas, narr1_names] <- nn1.narr1[, paste0("nn_", narr1_names)]
  x_n4[narr1_nas, narr1_names] <- nn4.narr1[, paste0("nn_", narr1_names)]
  
  narr2_nas <- which(is.na(x2[, narr2_names[1]]))
  x2.na2 <- x2[narr2_nas, ]
  x2.obs2 <- x2[-narr2_nas, ]
  
  nn1.narr2 <- add_nearestN(x2.na2, x2.obs2, N.nn = 1, nth = 1, 
                            var.names = narr2_names)
  nn4.narr2 <- add_nearestN(x2.na2, x2.obs2, N.nn = 4, nth = 1, 
                            var.names = narr2_names)
  
  # Test: name order matches
  expect_identical(
    str_sub(names(nn1.narr2), 4)[1:length(narr2_names)],
    narr2_names
  )
  
  # Replace missing values
  x_n1[narr2_nas, narr2_names] <- nn1.narr2[, paste0("nn_", narr2_names)]
  x_n4[narr2_nas, narr2_names] <- nn4.narr2[, paste0("nn_", narr2_names)]
  
  # Test that there are now no NAs
  expect_identical(0L, sum(is.na(x_n1[, narr1_names])))
  expect_identical(0L, sum(is.na(x_n4[, narr2_names])))
  
  # Test that where there were no missing values, things are the same.
  narr_all_nas = union(narr1_nas, narr2_nas)
  expect_equal(
    x_n1[-narr_all_nas, ] %>% 
      dplyr::select(starts_with("narr")),
    x[-narr_all_nas, ] %>% 
      dplyr::select(starts_with("narr"))
  )
  expect_equal(
    x_n4[-narr_all_nas, ] %>% 
      dplyr::select(starts_with("narr")),
    x[-narr_all_nas, ] %>% 
      dplyr::select(starts_with("narr"))
  )
  expect_true(all_equal(
      x_n1[-narr_all_nas, ] %>% 
        dplyr::select(starts_with("narr")),
      x[-narr_all_nas, ] %>% 
        dplyr::select(starts_with("narr")),
      ignore_col_order = F, ignore_row_order = F
    )
  )
  expect_true(all_equal(
      x_n4[-narr_all_nas, ] %>% 
        dplyr::select(starts_with("narr")),
      x[-narr_all_nas, ] %>% 
        dplyr::select(starts_with("narr")),
      ignore_col_order = F, ignore_row_order = F
    )
  )
  
  # Save using saveRDS -- for now I specify version = 2 to ensure backwards
  #  compatibility with versions of R prior to 3.5.0
  write_n1 <- paste0("rawrds/n1_day", day, ".rds")
  print(paste("Write name:", write_n1))
  saveRDS(x_n1, write_n1, version = 2)
  
  write_n4 <- paste0("rawrds/n4_day", day, ".rds")
  print(paste("Write name:", write_n4))
  saveRDS(x_n4, write_n4, version = 2)

  # Remove objects at the end of each iteration
  rm(x, x_n1, x_n4, write_n1, write_n4, read_name,
     x2, narr1_names, narr1_nas, narr2_names, narr2_nas,
     x2.na1, x2.obs1, x2.na2, x2.obs2,
     nn1.narr1, nn4.narr1, nn1.narr2, nn4.narr2, narr_all_nas, narr_list)  
}
