library(testthat)
library(fields)
library(dplyr)

minN <- function(x, N=2){
  # Finds the Nth smallest value of vector x
  # I modify this from the following:
  # https://stackoverflow.com/a/53144760
  len <- length(x)
  if (N>len) {
    stop('N greater than length(x).')
  }
  sort(x, partial=N)[N]
}

add_nearest_value <- function(data.1, data.2, nth = 2,
                              var.names = "aod_value") {
  
  # Finds the nearest point's value for obs in data.1 in obs in data.2
  # By default, this function ASSUMES that data.1 is contained in data.2
  # Change nth = 1 if they are disjoint. 
  
  # NOTE: This function requires a lot of RAM
  # Brief estimations of memory needed:
  # -- 6.4 GB for 29,397 rows in matrix
  # -- 21 GB for 53,807 rows in matrix
  # Will crash computer if you are *near* the limit based on experience.
  # fields will throw an error beforehand if it tries to allocate a matrix 
  #  that is too large (prevents crash).
  
  # data.1 contains the observaions we want the information for
  # data.2 contains the full dataset from which we are judging distance
  
  # Get distance matrix between data.1 and data.2  
  dist.mat <- 
    fields::rdist(as.matrix(data.1 %>% dplyr::select(cmaq_x, cmaq_y)),
                  as.matrix(data.2 %>% dplyr::select(cmaq_x, cmaq_y)))
  
  # Find the index in data.2 for the 2nd nearest point
  #    (first nearest point is the point itself)
  dist.min2 <- apply(dist.mat, 1, 
                     function(x) {
                       return(minN(x, N=nth))
                     })
  index.min2 <- apply(dist.mat, 1, 
                      function(x) {
                        return(match(minN(x, N=nth), x))
                      })
  df.min2 <- data.2[index.min2, var.names, drop = F]
  
  # Output: data.frame with the stuff.
  # - Nearest neighbor AOD
  # - Distance to nearest neighbor AOD
  df.min2 <- cbind(df.min2, dist.min2)
  colnames(df.min2) <- c(paste0("nn_", var.names), "nn_distance")
  
  return(df.min2)      
}



add_nearestN <- function(data.1, data.2, N.nn = 4, nth = 2,
                         xy.name = c("cmaq_x", "cmaq_y"),
                         var.names = "aod_value") {
  
  # FIXME
  # Change nth to Nth or something else -- nth is a function in dplyr
  # Works but will avoid confusion if you can change
  
  # FIXME
  # Note -- if seeking nearest value, and multiple values match, only the 
  #   first value (in order of appearance) is used.
  # Similarly, if we ask for N.nn = 4, and the 4th and 5th nearest values 
  #   are the same, we use only the 4th (in order of appearance). 
  # May change this appearance later but will not likely impact our intended
  #  application with the CMAQ grid projected data. 
  
  # Function that averages nearest N points: weights are prop. to 1/(d^2)
  # If N.nn = 1, this simplifies to just using the nearest point. 
  
  # NOTE: This function requires a lot of RAM
  # Brief estimations of memory needed:
  # -- 6.4 GB for 29,397 rows in matrix
  # -- 21 GB for 53,807 rows in matrix
  # Will crash computer if you are *near* the limit based on experience.
  # fields will throw an error beforehand if it tries to allocate a matrix 
  #  that is too large (prevents crash).
  
  # data.1 contains the observaions we want the information for
  # data.2 contains the full dataset from which we are judging distance
  
  # Get distance matrix between data.1 and data.2  
  dist.mat <- fields::rdist(as.matrix(data.1[, xy.name]),
                            as.matrix(data.2[, xy.name]))
  
  # Loop and find the nearest N points
  # Record distances and indices in two matrices
  
  dist.nn = weight.nn = index.nn = matrix(NA, nrow = nrow(data.1), ncol = N.nn)
  class(dist.nn) = class(weight.nn) = "numeric"
  class(index.nn) <- "integer"

  # Use sorted.list to get the indices for N.nn values (nth to N.nn + nth)
  for (ri in seq(nrow(dist.mat))) {
    index.nn[ri, ] <- sort.list(dist.mat[ri, ])[nth:(nth+N.nn-1)]
    dist.nn[ri, ] <- dist.mat[ri, index.nn[ri, ]]
  }
  weight.nn <- 1 / (dist.nn ^ 2)

  # Weights are proportional to 1/dist^2
  weight.norm <- base::rowSums(weight.nn)
  weight.adjust <- weight.nn / weight.norm
  expect_equal(rowSums(weight.adjust), rep(1, nrow(data.1)))
  
  # Loop through variable names and construct the weighted average
  var.matrix <- matrix(0, nrow = nrow(data.1), ncol = length(var.names))
  class(var.matrix) <- "numeric"
  for (vi in seq(length(var.names))) {
    v.name = var.names[vi]
    for (ci in seq(ncol(weight.adjust))) {
      var.matrix[, vi] = var.matrix[, vi] + 
        weight.adjust[, ci]*data.2[index.nn[, ci], v.name, drop = T]
    }
  }
  
  # weighted.distance
  weighted.distance <- base::rowSums(weight.adjust * dist.nn)
  # Output: data.frame
  # - Nearest neighbor value of variables
  # - Distance to nearest neighbor
  df.nn <- data.frame(var.matrix, dist.nn, weighted.distance)
  colnames(df.nn) <- c(paste0("nn_", var.names), 
                       paste0("nn_distance", seq(N.nn)),
                       "weighted_distance")
  
  return(df.nn)      
}