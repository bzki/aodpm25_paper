context("Testing getHole function")
source("../functions/hole.R")

library(testthat)
library(dplyr)
library(Rfast)

test_that("Make sure excluding points correctly", {
  
  xy <- readRDS("../rawrds/n_aod_pred_2011_new8_194.rds")  %>%
    dplyr::select(cmaq_x, cmaq_y)
  
  # dist_xy <- fields::rdist(as.matrix(xy[1:100, ]), as.matrix(xy))
  # Note: cmaq_x and cmaq_y are in meters -- I look within 100km distance
  set.seed(1212)
  function_output <- getHoles.fast(xy, num.holes = 10, dist.hole = 100 * 1000)
  hole.index <- function_output$omitted.index
  set.seed(1212)
  initial.index <- sample.int(n = nrow(xy), size = 10)
  expect_identical(initial.index, function_output$initial.index)
  # Check distance between the excluded and the initial index.
  # -- should be within 100km
  dist_1 <- fields::rdist(as.matrix(xy[hole.index, ]),
                          as.matrix(xy[initial.index, ]))
  # Get row minimums
  val_1 <- Rfast::rowMins(dist_1, value = T)
  # These should be less than 100*1000
  expect_lt(max(val_1), 100*1000)
  
  # The values left-in should be greater than 100*1000 in distance
  dist_2 <- fields::rdist(as.matrix(xy[-hole.index, ]),
                          as.matrix(xy[initial.index, ]))
  val_2 <- Rfast::rowMins(dist_2, value = T)
  expect_gt(max(val_2), 100*1000)
})
