context("Testing nearest-neighbor funcions")
source("../functions/nearest.R")

library(testthat)

test_that("Checking minN", {
  vec_1 <- c(3L, 1L, 2L)
  vec_2 <- c(1.000003, 1.000002, 1.000001) 

  expect_identical(minN(vec_1, N=1), 1L)
  expect_identical(minN(vec_1, N=2), 2L)
  expect_identical(minN(vec_1, N=3), 3L)
  
  expect_identical(match(minN(vec_1, N=1), vec_1), 2L)
  expect_identical(match(minN(vec_1, N=2), vec_1), 3L)
  expect_identical(match(minN(vec_1, N=3), vec_1), 1L)
  
  expect_identical(minN(vec_2, N=1), 1.000001)
  expect_identical(minN(vec_2, N=2), 1.000002)
  expect_identical(minN(vec_2, N=3), 1.000003)

  expect_identical(match(minN(vec_2, N=1), vec_2), 3L)
  expect_identical(match(minN(vec_2, N=2), vec_2), 2L)
  expect_identical(match(minN(vec_2, N=3), vec_2), 1L)
  
})


test_that("Check add_nearestN, where NA/Obs parts separated (nth = 1)", {

  test.df <- data.frame(
    id = seq(5), 
    x = rep(1, 5),
    y = seq(5),
    var1 = c(NA, 2, 4, 6, NA),
    var2 = c(1.2, 2.2, NA, 6.2, 7.2),
    var3 = c(NA, NA, 5.0, 6.0, 7.0),
    var4 = c(NA, NA, 2.3, 4.3, 2.0)
  )
  
  # Fill in var1 and var2 separately (different missingness patterns)
  var1.na = which(is.na(test.df$var1))
  df.na.1 <- test.df[var1.na, ]
  df.obs.1 <- test.df[-var1.na, ]
  
  var2.na = which(is.na(test.df$var2))
  df.na.2 <- test.df[var2.na, ]
  df.obs.2 <- test.df[-var2.na, ]
  
  rev.df = test.df %>% arrange(desc(id))
  var2r.na = which(is.na(rev.df$var2))
  df.na.2r <- rev.df[var2r.na, ]
  df.obs.2r <- rev.df[-var2r.na, ]
  
  # Test that N.nn = 1 correctly picks the values
  nn.var1 <- add_nearestN(df.na.1, df.obs.1, N.nn = 1, nth = 1,
                          xy.name = c("x", "y"), var.names = "var1")
  
  expect_identical(as.numeric(nn.var1$nn_var1), c(2, 6))
  expect_identical(as.numeric(nn.var1$nn_distance1), c(1, 1))
  expect_identical(as.numeric(nn.var1$weighted_distance), c(1,1))
  
  nn.var2 <- add_nearestN(df.na.2, df.obs.2, N.nn = 1, nth = 1,
                          xy.name = c("x", "y"), var.names = "var2")

  # Picks the earlier coordinate because it appears first.
  expect_identical(as.numeric(nn.var2$nn_var2), 2.2)
  expect_identical(as.numeric(nn.var2$nn_distance1), c(1))
  expect_identical(as.numeric(nn.var2$weighted_distance), c(1))
  
  nn.var2r <- add_nearestN(df.na.2r, df.obs.2r, N.nn = 1, nth = 1,
                           xy.name = c("x", "y"), var.names = "var2")
  
  # Picks 6.2 because it now appears first.
  expect_identical(as.numeric(nn.var2r$nn_var2), 6.2)
  
  # Try multiple variables.
  var34.na = which(is.na(test.df$var3))
  df.na = test.df[var34.na, ]
  df.obs = test.df[-var34.na, ]
  nn.var34 <- add_nearestN(df.na, df.obs, N.nn = 1, nth = 1,
                           xy.name = c("x", "y"), 
                           var.names = c("var3", "var4"))
  
  expect_identical(as.numeric(nn.var34$nn_var3), c(5.0, 5.0))
  expect_identical(as.numeric(nn.var34$nn_var4), c(2.3, 2.3))
  expect_identical(as.numeric(nn.var34$nn_distance1), c(2, 1))
  expect_identical(as.numeric(nn.var34$weighted_distance), c(2, 1))
  
  # Tests with N.nn = 2
  nn.var34.n2 <- add_nearestN(df.na, df.obs, N.nn = 2, nth = 1,
                              xy.name = c("x", "y"), 
                              var.names = c("var3", "var4"))
  # Weights should add to 1 and be proportional to (1/(d^2))
  tmp.sum1 = (1/(2^2)) + (1/(3^2))
  weights.1a = (1/(2^2))/tmp.sum1
  weights.1b = (1/(3^2))/tmp.sum1
  tmp.sum2 = (1/1) + (1/(2^2))
  weights.2a = (1/1)/tmp.sum2
  weights.2b = (1/(2^2))/tmp.sum2
  expect_identical(as.numeric(nn.var34.n2$nn_var3), 
                   c(5.0 * weights.1a + 6.0 * weights.1b,
                     5.0 * weights.2a + 6.0 * weights.2b))
  expect_identical(as.numeric(nn.var34.n2$nn_var4), 
                   c(2.3 * weights.1a + 4.3 * weights.1b,
                     2.3 * weights.2a + 4.3 * weights.2b))
  expect_identical(as.numeric(nn.var34.n2$nn_distance1), c(2, 1))
  expect_identical(as.numeric(nn.var34.n2$nn_distance2), c(3, 2))
  expect_identical(as.numeric(nn.var34.n2$weighted_distance),
                   c(weights.1a * 2 + weights.1b * 3,
                     weights.2a * 1 + weights.2b * 2))
  
  # Check for single variable when equidistant to 2 points
  nn.var2.n2 <- add_nearestN(df.na.2, df.obs.2, N.nn = 2, nth = 1,
                             xy.name = c("x", "y"), var.names = c("var2"))
  expect_identical(as.numeric(nn.var2.n2$nn_var2), 
                   c(1/2 * 2.2 + 1/2 * 6.2))
  expect_identical(as.numeric(nn.var2.n2$weighted_distance),
                   c(1/2 * 1 + 1/2 * 1))
  
})

test_that(
  "Check that add_nearestN and add_nearest_value produce same results", {
  x <- readRDS("../rawrds/n_aod_pred_2011_new8_182.rds")  %>%
    select(cmaq_id, cmaq_x, cmaq_y, aod_value)
  x_na <- x %>%
    dplyr::filter(is.na(aod_value)) %>%
    dplyr::slice(1:500)
  x_obs <- x %>%
    dplyr::filter(!is.na(aod_value))
  
  ret_1a <- add_nearest_value(x_na, x_obs, nth = 1, var.names = "aod_value")
  ret_2a <- add_nearestN(x_na, x_obs, N.nn = 1, nth = 1, var.names = "aod_value")
  
  expect_identical(ret_1a$nn_aod_value, ret_2a$nn_aod_value)
  expect_identical(ret_1a$nn_distance, ret_2a$nn_distance1)
  
  ret_1b <- add_nearest_value(x_obs[1:100, ], x_obs, nth = 2, 
                              var.names = "aod_value")
  ret_2b <- add_nearestN(x_obs[1:100, ], x_obs, N.nn = 1, nth = 2, 
                         var.names = "aod_value")
  expect_identical(ret_1b$nn_aod_value, ret_2b$nn_aod_value)
  expect_identical(ret_1b$nn_distance, ret_2b$nn_distance1)
  
  # Check for row 15 that you are getting the correct result.
  row_15 <- ret_1b[15, ]
  tmp_dist <- fields::rdist(as.matrix(x_obs[1:100, c("cmaq_x", "cmaq_y")]),
                            as.matrix(x_obs[, c("cmaq_x", "cmaq_y")]))[15, ]
  sorted_dist <- sort.list(tmp_dist)[2]
  expect_identical(row_15$nn_distance, tmp_dist[sorted_dist])
  expect_identical(row_15$nn_aod_value, x_obs$aod_value[sorted_dist])
})


test_that("Check add_nearestN for nth = 2", {

  # Distances will be unique for this data.frame -- avoid issue with ties
  df <- data.frame(
    id = seq(5), 
    x = c(1.0, -1.2, 3.7, 5.4, 18.2),
    y = seq(5),
    var3 = c(15.0, 7.2, 5.0, 6.0, 7.0),
    var4 = c(4.9, 0.9, 2.3, 4.3, 2.0)
  )
  
  dist.mat <- as.matrix(dist(as.matrix(df[, c("x", "y")])))
  sort.mat = dist.sort = matrix(NA, nrow = nrow(dist.mat), ncol = nrow(df) - 1)
  for (ri in seq(nrow(dist.mat))) {
    sort.mat[ri, ] <- sort.list(dist.mat[ri, ])[-1]
    dist.sort[ri, ] <- dist.mat[ri, sort.mat[ri, ]]
  }
  near.var3 <- df$var3[sort.mat[, 1]]
  near.var4 <- df$var4[sort.mat[, 1]]
  
  # Fit function with nth = 2
  nn_1 <- add_nearestN(df, df, N.nn = 1, nth = 2, xy.name = c("x", "y"),
                       var.names = c("var3", "var4"))
  expect_identical(as.numeric(nn_1$nn_distance1), dist.sort[, 1])
  expect_identical(as.numeric(nn_1$nn_var3), near.var3)
  expect_identical(as.numeric(nn_1$nn_var4), near.var4)
})

test_that("Test that multiple variables works as expected", {
  
  # Distances will be unique for this data.frame -- avoid issue with ties
  df <- data.frame(
    id = seq(5), 
    x = c(1.0, -1.2, 3.7, 5.4, 18.2),
    y = seq(5),
    var3 = c(15.0, 7.2, 5.0, 6.0, 7.0),
    var4 = c(4.9, 0.9, 2.3, 4.3, 2.0)
  )
  
  nn_4 <- add_nearestN(df, df, N.nn = 4, nth = 2, xy.name = c("x", "y"),
                       var.names = c("var3", "var4"))
  nn_4_alt <- add_nearestN(df, df, N.nn = 4, nth = 2, xy.name = c("x", "y"),
                           var.names = c("var4", "var3"))
  nn_1_var3 <- add_nearestN(df, df, N.nn = 4, nth = 2, xy.name = c("x", "y"),
                            var.names = c("var3"))
  nn_1_var4 <- add_nearestN(df, df, N.nn = 4, nth = 2, xy.name = c("x", "y"),
                            var.names = c("var4"))
  
  expect_identical(nn_4$nn_var3, nn_1_var3$nn_var3)  
  expect_identical(nn_4$nn_var4, nn_1_var4$nn_var4)  
  expect_identical(nn_4_alt$nn_var3, nn_1_var3$nn_var3)  
  expect_identical(nn_4_alt$nn_var4, nn_1_var4$nn_var4)  
  expect_identical(names(nn_4_alt)[1:2], c("nn_var4", "nn_var3"))
  
})

