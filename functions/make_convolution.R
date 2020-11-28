library(fields)

make_convolution <- function(df_input, x_name, y_name, var_name, 
                             fold_list, scale_num = 100000) {
  
  # Return vector with convolution of PM2.5 (or any other variable)
  # fold_list: list of k lists for k-fold CV
  # -- each sublist has [[1]] training 
  
  # Prepare a numeric matrix that will contain the k columns
  ret_mat <- matrix(NA, nrow = nrow(df_input), ncol = length(fold_list))
  colnames(ret_mat) <- paste0("f", seq(length(fold_list)))
  class(ret_mat) <- "numeric"
  
  for (fi in seq(length(fold_list))) {
  
    train_index <- fold_list[[fi]][[1]]
    train_value <- df_input[train_index, var_name]
    
    # Distance between all data and the training set only
    # (train, test) x (train)
    dist_mat <- 
      rdist(as.matrix(df_input[, c(x_name, y_name)]),
            as.matrix(df_input[train_index, c(x_name, y_name)]))
    # Replace 0s with NA -- should only impact training rows
    dist_mat[dist_mat == 0] <- NA
    testthat::expect_equal(sum(is.na(dist_mat)), length(train_index))
    # Use scaling number to preserve numerical stability
    invdist_mat <- scale_num * 1 / (dist_mat^2) 
    sum_weights <- rowSums(invdist_mat, na.rm = T)
    # Create the weights which sum to 1 for each row
    weight_mat <- invdist_mat/sum_weights
    testthat::expect_equivalent(rowSums(weight_mat, na.rm = T),
                                rep(1, nrow(df_input)))

    # Create convolution layer
    for (ri in seq(nrow(df_input))) {
      ret_mat[ri, fi] <- 
        sum(weight_mat[ri, ] * train_value, na.rm = T) /
        sum(weight_mat[ri, ], na.rm = T)
    } 
  }
  return(ret_mat)
}
