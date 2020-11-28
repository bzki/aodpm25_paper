library(fields)

getHoles.fast <- function(xy, num.holes, dist.hole) {
  # This function creates holes within set of coordinates given
  # RETURN:
  #  - omitted.index: vector indexing which of the xy that we will omit.
  #  - initial.index: vector indexing the center points of the circles omitted
  # ARGUMENTS
  # - xy: matrix of xy coordinates (first column x, second column y)
  # - num.holes: Number of holes we want to exclude
  # - dist.hole: After picking the holes, we exclude any coordinate within
  #     dist.hole distance of the holes. 
  # Note -- sample() functions are impacted by the change in R 3.6
  initial.index <- sample.int(n = nrow(xy), size = num.holes)
  omit.list <- vector(mode = "list", length = num.holes)
  for (h in seq(num.holes)) {
    distances <- as.vector(fields::rdist( 
      matrix(xy[initial.index[h], ], nrow = 1), xy) )
    omit.index <- which(distances < dist.hole)
    omit.list[[h]] <- omit.index
  }
  # Obtain the union of all index values we want to omit using Reduce()
  omitted.index <- Reduce(union, omit.list)
  return(list(omitted.index = omitted.index,
              initial.index = initial.index))
}