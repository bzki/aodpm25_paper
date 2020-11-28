# Chunked ranger prediction to deal with memory issues when predicting many obs
library(ranger)

rpred_chunked <- function(fit_object, new_data, chunk_size = 50000, 
                          verbose = TRUE) {
  
  start_row <- 1
  end_row <- min(chunk_size, nrow(new_data))
  pred_value <- as.numeric(rep(NA, nrow(new_data)))
  
  while (end_row <= nrow(new_data)) {
    if (verbose) {
      print(paste0("Prediction for ", start_row, " to ", end_row))
    }
    
    pred_value[start_row:end_row] <- 
      predict(fit_object, new_data[start_row:end_row, ])$predictions
    
    if (end_row == nrow(new_data)) {
      break
    } else {
      start_row <- start_row + chunk_size
      end_row <- min(end_row + chunk_size, nrow(new_data))
    }
  }
  
  return(pred_value)
}


