# Chunked nearest neighbor distance calculation
# This is done to save space as the distance matrix can be extremely large
# TODO: Write tests

nn_chunked <- function(x1, x2, chunk_size = 5000, verbose = FALSE) {
  # Returns nearest neighbor distance in x2 for each row in x1
  # -- assumes coordinates are already projected (i.e., not great circle dist)
  start_row <- 1
  end_row <- min(chunk_size, nrow(x1))
  nn_dist <- NULL
  while (end_row <= nrow(x1)) {
    if (verbose) {
      print(paste0("Distance for ", start_row, " to ", end_row))
    }
    nn_dist <- 
      c(nn_dist, 
        apply(
          fields::rdist(x1[start_row:end_row, , drop = F], x2), 
          1, min)
      )
    if (end_row == nrow(x1)) {
      break
    } else {
      start_row <- start_row + chunk_size
      end_row <- min(end_row + chunk_size, nrow(x1))
    }
  }
  return(nn_dist)
}

