
#########################  Generating Null Statistics   ##################

estimate_propensity <- function (grna_assignment_matrix, covariate) {

  estimated_prop <- matrix(NA, nrow = nrow(grna_assignment_matrix), ncol = ncol(grna_assignment_matrix))

  # Loop through each row of the grna_assignment_matrix
  for (i in 1:nrow(grna_assignment_matrix)) {

    # Extract the current row from the grna_assignment_matrix
    grna_row <- as.numeric(grna_assignment_matrix[i, ])

    # Fit logistic regression model using covariates and the current row
    glm_fit <- glm(grna_row ~ ., family = binomial(link = "logit"), data = covariate)

    # Predict the fitted means (estimated probabilities) for this row
    estimated_prop[i, ] <- predict(glm_fit, type = "response")
  }
  return(estimated_prop)
}

resample_and_combine<- function(grna_assignment_matrix, n_permute) {
  propensity_matrix <- estimate_propensity(grna_assignment_matrix , covariate)

  propensity <- do.call(rbind, replicate(n_permute, propensity_matrix, simplify = FALSE))


  resampled_mat <- matrix(rbinom(n = length(propensity),
                                 size = 1,
                                 prob = as.vector(propensity)),
                          nrow = nrow(propensity),
                          ncol = ncol(propensity))

  rownames(resampled_mat) <- rep(rownames(grna_assignment_matrix),n_permute)
  combined_mat <- rbind(grna_assignment_matrix, resampled_mat)

  return(combined_mat)
}

permute_and_combine <- function(grna_assignment_matrix, n_permute) {

  permute_list <- list()

  for (i in 1:n_permute) {
    # Generate a permutation of the column indices
    permuted_indices <- sample(ncol(grna_assignment_matrix))

    # Apply the permutation to the columns
    permute_list[[i]] <- grna_assignment_matrix[, permuted_indices]

  }
  combined_mat <-do.call(rbind, permute_list)
  combined_matrix_with_original <- rbind(grna_assignment_matrix, combined_mat)
  return(combined_matrix_with_original)
}

permutation_or_resampling <- function(grna_assignment_matrix, n_permute, use_resample=F){
  if(use_resample){
    return(resample_and_combine(grna_assignment_matrix, n_permute))
  } else{
    return(permute_and_combine(grna_assignment_matrix, n_permute))
  }
}


################### Get test statistics #############################

t_test <- function(response_matrix, grna_assignment_matrix) {

  response_matrix <- as.matrix(response_matrix)
  grna_assignment_matrix <- as.matrix(grna_assignment_matrix)

  n_cells <- ncol(response_matrix)
  n_genes <- nrow(response_matrix)
  all_cells <- seq_len(n_cells)
  n_rows <- nrow(grna_assignment_matrix)

  response_matrix_t <- t(response_matrix)
  response_matrix_t_squared <-response_matrix_t^2

  grna_infected_cell_count <- rowSums(grna_assignment_matrix)

  mean_expression_trt <- grna_assignment_matrix %*% response_matrix_t / grna_infected_cell_count
  mean_expression_squared_trt <- grna_assignment_matrix %*% response_matrix_t_squared / grna_infected_cell_count
  var_expression_trt <- mean_expression_squared_trt - mean_expression_trt^2

  grna_control_cell_count <- n_cells - grna_infected_cell_count
  mean_expression_control <- (1 - grna_assignment_matrix) %*% response_matrix_t / grna_control_cell_count
  mean_expression_squared_control <- (1 - grna_assignment_matrix) %*% response_matrix_t_squared / grna_control_cell_count
  var_expression_control <- mean_expression_squared_control - mean_expression_control^2

  t_statistics <- (mean_expression_trt - mean_expression_control) / sqrt(var_expression_trt / grna_infected_cell_count + var_expression_control / grna_control_cell_count)

  return(abs(t_statistics))
}

################### Use gene modules ###########################

convert_gene_list_to_matrix <-function(response_matrix, gene_module){

  genes_in_response <- rownames(response_matrix)

  # Initialize a matrix to store the result
  result_matrix <- matrix(0, nrow = nrow(response_matrix), ncol = length(gene_module))

  # Loop over each gene module (i.e., each list of genes in gene_module)
  for (i in seq_along(gene_module)) {
    # Get the current gene module (a vector of genes)
    current_module <- gene_module[[i]]

    # For each gene in the response_matrix, check if it's in the current gene module
    result_matrix[, i] <- as.numeric(genes_in_response %in% current_module)
  }

  colnames(result_matrix) <- names(gene_module)

  # Return the resulting matrix
  return(result_matrix)
}

combine_columns_from_matrices <- function(matrices,gene_module_matrix) {
  # Combine columns from the input matrices column by column
  combined <- do.call(cbind, lapply(1:ncol(matrices[[1]]), function(i) {
    # Extract the i-th column from each matrix and combine them
    cbind(sapply(matrices, function(mat) mat[, i]))
  }))
  colnames(combined) <- rep(colnames(gene_module_matrix),each=length(matrices))
  return(combined)
}


combine_test_statistics <- function(test_result, gene_module_matrix){

  test_result <- abs(test_result)

  power_sum_1 <- test_result%*%gene_module_matrix
  power_sum_2 <- test_result^2%*%gene_module_matrix
  power_sum_3 <- test_result^3%*%gene_module_matrix
  power_sum_4 <- test_result^4%*%gene_module_matrix
  power_sum_5 <- test_result^5%*%gene_module_matrix
  power_sum_6 <- test_result^6%*%gene_module_matrix
  max_stat <- apply(gene_module_matrix, 2, function(module) {
    apply(test_result * module, 1, max)
  })


  combined_matrix <- combine_columns_from_matrices(list(power_sum_1,power_sum_2,power_sum_3,power_sum_4,power_sum_5,power_sum_6, max_stat), gene_module_matrix)
  # combined_matrix <- combine_columns_from_matrices(list(max_stat), gene_module_matrix)

  return(combined_matrix)
}

############################### Convert test matrix to tensor ###############################

convert_matrix_to_tensor <- function(matrix, n_permute, grna_assignment_matrix) {

  # Get the number of rows in grna_assignment_matrix, this will be the number of slices (depth)
  depth <- nrow(grna_assignment_matrix)

  # Initialize an empty tensor with dimensions (n_permute x ncol(matrix) x depth)
  tensor <- array(0, dim = c((n_permute+1), ncol(matrix), depth))

  # Loop over the number of slices (depth)
  for (i in 1:depth) {
    # For each slice, extract the corresponding rows for each group (NTC-1, NTC-2, ..., NTC-n_permute)
    row_indices <- seq(i, nrow(matrix), by = depth)

    # Assign the extracted rows to the i-th slice of the tensor
    tensor[,,i] <- matrix[row_indices, ]
  }

  dimnames(tensor) <- list(

    c("Original", paste0("Permutation", 1:n_permute)),
    # Second dimension: Column names from the input matrix
    colnames(matrix),

    # Third dimension: Slice names based on grna_assignment_matrix (e.g., the row names or indices)
    rownames(grna_assignment_matrix)
  )
  return(tensor)
}
#################################### Compute p value ##############################


rank_tensor <- function(tensor) {
  # Get the dimensions of the input tensor
  dim_tensor <- dim(tensor)

  # Initialize an empty tensor to store the ranked values with the same dimensions
  ranked_tensor <- array(NA, dim = dim_tensor, dimnames = dimnames(tensor))

  # Use mclapply to parallelize the operation across the 3rd dimension (slices)
  ranked_slices <- mclapply(1:dim_tensor[3], function(i) {
    # Rank each slice (2D matrix) of the tensor
    ranked_slice <- rank_matrix(t(tensor[,,i]))$X_ranked
    return(ranked_slice)
  }, mc.cores = parallel::detectCores() - 1)  # Using all available cores except one

  # Assign the ranked slices back to the ranked_tensor
  for (i in 1:dim_tensor[3]) {
    ranked_tensor[,,i] <- ranked_slices[[i]]
  }

  return(ranked_tensor)
}


compute_p_value_from_tensor <- function(tensor,n_permute,side_code){

  rank<- rank_tensor(tensor)
  if (side_code==1){ #right side
    rank <- (n_permute+2-rank)
  }
  rank[rank == (n_permute + 1)] <- n_permute
  return(rank/(n_permute+1))
}


min_p_val <- function(p_val_tensor, gene_module_index_matrix) {
  # Get the dimensions of the input tensor
  dim_tensor <- dim(p_val_tensor)

  n_groups <- ncol(gene_module_index_matrix)

  # Calculate the number of columns per group
  cols_per_group <- dim_tensor[2] / n_groups

  # Initialize an empty tensor to store the results with the new dimensions
  min_p_tensor <- array(NA, dim = c(dim_tensor[1], n_groups, dim_tensor[3]),
                        dimnames = list(dimnames(p_val_tensor)[[1]],
                                        colnames(gene_module_index_matrix),
                                        dimnames(p_val_tensor)[[3]]))

  # Use lapply to process the 3rd dimension (slices)
  results <- lapply(1:dim_tensor[3], function(i) {
    current_slice <- p_val_tensor[,,i]
    result_slice <- matrix(NA, nrow = dim_tensor[1], ncol = n_groups)

    # Loop over each group
    for (j in 1:n_groups) {
      # Define start and end columns for the current group
      start_col <- (j - 1) * cols_per_group + 1
      end_col <- j * cols_per_group

      # Ensure start_col and end_col are integers
      start_col <- as.integer(start_col)
      end_col <- as.integer(end_col)

      # Slice the appropriate columns for the current group
      current_group <- current_slice[, start_col:end_col, drop = FALSE]

      # Calculate the row-wise minimum for this group
      result_slice[, j] <- apply(current_group, 1, min, na.rm = TRUE)
    }

    return(result_slice)
  })

  # Combine the results from lapply into the min_p_tensor
  for (i in 1:dim_tensor[3]) {
    min_p_tensor[,,i] <- results[[i]]
  }

  return(min_p_tensor)
}

compute_final_p_val <- function(min_p_tensor, n_permute) {
  # Get the dimensions of the tensor
  dim_tensor <- dim(min_p_tensor)

  # Initialize an empty matrix to store the p-value counts
  final_p_vals <- matrix(0,
                         nrow = dim_tensor[3],
                         ncol = dim_tensor[2],
                         dimnames = list(dimnames(min_p_tensor)[[3]], dimnames(min_p_tensor)[[2]]))

  # Use a loop to iterate over each slice of the tensor (rows)
  for (i in 1:dim_tensor[3]) {
    # Extract the current slice (2D matrix)
    current_slice <- min_p_tensor[,,i]

    # Apply rank with ties.method set to "min" and calculate p-values
    final_p_vals[i, ] <- apply(current_slice, 2, function(col) rank(col, ties.method = "min")[1] / (n_permute + 1))
  }

  return(final_p_vals)
}
