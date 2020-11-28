context("Testing make_convolution function")
source("../functions/make_convolution.R")

library(testthat)
library(dplyr)
library(fields)


test_that("Hand example", {
  
  data_mat <- as.data.frame(
    cbind(matrix(c(1, 2, 0, 2, 3, 2, 4, 2), byrow = T, nrow = 4),
          c(1, 10, 20, 25)))
  names(data_mat) <- c("x", "y", "z")
  
  # First two observations as test observations, last 2 as training
  list_train <- list(
    list(c(3, 4))
  )
  convFromFunction <- 
    make_convolution(data_mat, x_name = "x", y_name = "y",
                     var_name = "z", fold_list = list_train)
  
  # Distance for obs 1 from obs 3 is: 2
  # Distance for obs 1 from obs 4 is: 3
  conv_1 <- (1/(2^2)*data_mat$z[3] + 1/(3^2)*data_mat$z[4]) /
    (1/(2^2) + 1/(3^2))
  # Distance for obs 2 from obs 3 is: 3
  # Distance for obs 2 from obs 4 is: 4
  conv_2 <- (1/(3^2)*data_mat$z[3] + 1/(4^2)*data_mat$z[4]) /
    (1/(3^2) + 1/(4^2))
  # Distance for obs 3 from obs 4 is: 1
  # Distance for obs 4 fro obs 3 is: 1
  conv_3 <- (data_mat$z[4])
  conv_4 <- (data_mat$z[3])
  conv_manual <- c(conv_1, conv_2, conv_3, conv_4)  
  testthat::expect_equivalent(
    conv_manual,
    convFromFunction[, 1]
  )  
  
})

test_that("Data example", {
  
  # Read in a day's worth of data and create x/y in kilometers
  x <- readRDS("../rawrds/n_aod_pred_2011_new8_194.rds") %>%
    filter(!is.na(pm25_value)) %>%
    select(cmaq_id, cmaq_x, cmaq_y, pm25_value) %>%
    mutate(x_km = cmaq_x/1000, y_km = cmaq_y/1000)
  
  # Create fold list
  list_1 <- list(
    list(c(11:nrow(x))), # Test observations: 1 through 10
    list(c(1:10, 21:nrow(x))) # Test observations: 11 through 20
  )
  
  # Output from function
  convFromFunction <- 
    make_convolution(x, x_name = "x_km", y_name = "y_km",
                     var_name = "pm25_value", fold_list = list_1)
  
  # Focus on first 10 observations in first column first
  dist_1 <- rdist(x[1:10, c("x_km", "y_km")],
                  x[11:nrow(x), c("x_km", "y_km")])
  invdist_1 <- 1 / (dist_1^2)
  
  expect_equivalent(
    sum(invdist_1[1, ] * x$pm25_value[11:nrow(x)]) / sum(invdist_1[1, ]),
    convFromFunction[1, 1]
  )
  
  expect_equivalent(
    as.vector(
      invdist_1 %*% matrix(x$pm25_value[11:nrow(x)], ncol = 1) / 
        rowSums(invdist_1)
    ),
    convFromFunction[1:10, 1]
  )
  
  dist_2 <- rdist(x[11:20, c("x_km", "y_km")],
                  x[c(1:10, 21:nrow(x)), c("x_km", "y_km")])
  invdist_2 <- 1 / (dist_2^2)
  expect_equivalent(
    as.vector(
      invdist_2 %*% matrix(x$pm25_value[c(1:10, 21:nrow(x))], ncol = 1) / 
        rowSums(invdist_2)
    ),
    convFromFunction[11:20, 2]
  )
  
})

