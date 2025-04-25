
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

resample_and_combine<- function(grna_assignment_matrix, n_permute, covariate) {
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

permute_and_combine <- function(matrix_input, n_permute = 1, seed = 123) {
  set.seed(seed)  # Set seed for reproducibility

  # Initialize a list to store the original and shuffled matrices
  permuted_matrices <- list(matrix_input)

  # Generate n_permute shuffled versions
  for (i in 1:n_permute) {
    shuffled_matrix <- matrix_input[, sample(ncol(matrix_input))]  # Shuffle columns
    permuted_matrices[[i + 1]] <- shuffled_matrix  # Store the shuffled matrix
  }

  # Combine the original and shuffled matrices using rbind
  final_matrix <- do.call(rbind, permuted_matrices)

  return(final_matrix)
}


permutation_or_resampling <- function(grna_assignment_matrix, n_permute,covariate, use_resample=F){
  if(use_resample){
    return(resample_and_combine(grna_assignment_matrix, n_permute,covariate))
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

  return(t_statistics)
}





t_test_pooled_variance <- function(response_matrix, grna_assignment_matrix) {

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

  grna_control_cell_count <- n_cells - grna_infected_cell_count
  mean_expression_control <- (1 - grna_assignment_matrix) %*% response_matrix_t / grna_control_cell_count

  mean_expression <- colSums(response_matrix_t)/nrow(response_matrix_t)
  mean_expression_squared <- colSums(response_matrix_t^2)/nrow(response_matrix_t)

  var_expression <- mean_expression_squared-mean_expression^2


  t_statistics <- t(t(mean_expression_trt - mean_expression_control) / sqrt(var_expression))

  return(t_statistics)
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


test_normalizer <- function(gene_wise_test_tensor) {
  # Get the dimensions of the tensor
  dims <- dim(gene_wise_test_tensor)

  # Loop over all slices along the third dimension (j)
  for (j in seq_len(dims[3])) {
    for (i in seq_len(dims[2])) {
      # Extract the vector along the first dimension for gene_wise_test_tensor[, i, j]
      vec <- gene_wise_test_tensor[, i, j]

      # Calculate the median and MAD (robust measures)
      median_val <- median(vec, na.rm = TRUE)
      mad_val <- mad(vec, na.rm = TRUE)

      # Normalize the vector: (value - median) / MAD
      if (mad_val != 0) {
        gene_wise_test_tensor[, i, j] <- (vec - median_val) / mad_val
      } else {
        # If MAD is zero, set to NA to avoid division by zero
        gene_wise_test_tensor[, i, j] <- NA
      }
    }
  }

  # Return the normalized tensor
  return(gene_wise_test_tensor)
}

test_normalizer_vectorized <- function(gene_wise_test_tensor) {
  # 使用 apply 按基因维度计算
  normalized_tensor <- apply(gene_wise_test_tensor, c(2,3), function(x) {
    (x - median(x, na.rm=TRUE)) / mad(x, na.rm=TRUE)
  })
  # 恢复维度顺序
  aperm(normalized_tensor, c(2,1,3))
}




test_normalizer_sd <- function(gene_wise_test_tensor) {
  # Get the dimensions of the tensor
  dims <- dim(gene_wise_test_tensor)

  # Loop over all slices along the third dimension (j)
  for (j in seq_len(dims[3])) {
    for (i in seq_len(dims[2])) {
      # Extract the vector along the first dimension for gene_wise_test_tensor[, i, j]
      vec <- gene_wise_test_tensor[, i, j]

      # Calculate the mean and SD (standard deviation)
      mean_val <- mean(vec, na.rm = TRUE)
      sd_val <- sd(vec, na.rm = TRUE)

      # Standardize the vector: (value - mean) / SD
      if (sd_val != 0) {
        gene_wise_test_tensor[, i, j] <- (vec - mean_val) / sd_val
      } else {
        # If SD is zero, set to NA to avoid division by zero
        gene_wise_test_tensor[, i, j] <- NA
      }
    }
  }

  # Return the standardized tensor
  return(gene_wise_test_tensor)
}

##### Cpp:test_normalizer_sd_rcpp #######


cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// Function to compute mean
double compute_mean(NumericVector x) {
    int n = x.size();
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        sum += x[i];
    }

    return sum / n;
}

// Function to compute standard deviation
double compute_sd(NumericVector x, double mean_val) {
    int n = x.size();
    double sum_sq = 0.0;

    for (int i = 0; i < n; i++) {
        sum_sq += pow(x[i] - mean_val, 2);
    }

    return sqrt(sum_sq / (n - 1));  // Sample standard deviation
}

// Rcpp function for tensor standardization
// [[Rcpp::export]]
NumericVector test_normalizer_sd_rcpp(NumericVector tensor) {
    IntegerVector dims = tensor.attr("dim");  // Automatically get dimensions
    int d1 = dims[0], d2 = dims[1], d3 = dims[2];

    NumericVector result(clone(tensor));  // Make a copy to avoid modifying the original tensor

    for (int j = 0; j < d3; j++) {
        for (int i = 0; i < d2; i++) {
            // Extract slice manually
            NumericVector vec(d1);
            for (int k = 0; k < d1; k++) {
                vec[k] = tensor[k + i * d1 + j * d1 * d2];
            }

            // Compute mean and standard deviation
            double mean_val = compute_mean(vec);
            double sd_val = compute_sd(vec, mean_val);

            // Standardize and update result tensor
            for (int k = 0; k < d1; k++) {
                if (sd_val != 0) {
                    result[k + i * d1 + j * d1 * d2] = (vec[k] - mean_val) / sd_val;
                } else {
                    result[k + i * d1 + j * d1 * d2] = NA_REAL;  // Avoid division by zero
                }
            }
        }
    }

    result.attr("dim") = dims;  // Restore dimension attributes
    return result;
}')


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

  test_result[is.na(test_result)] <- 0

  test_result <- abs(test_result)


  power_sum_1 <- test_result%*%gene_module_matrix
  power_sum_2 <- (test_result^2%*%gene_module_matrix)^(1/2)
  power_sum_3 <- (test_result^3%*%gene_module_matrix)^(1/3)
  power_sum_4 <- (test_result^4%*%gene_module_matrix)^(1/4)
  power_sum_5 <- (test_result^5%*%gene_module_matrix)^(1/5)
  power_sum_6 <- (test_result^6%*%gene_module_matrix)^(1/6)
  power_sum_7 <- (test_result^7%*%gene_module_matrix)^(1/7)
  power_sum_8 <- (test_result^8%*%gene_module_matrix)^(1/8)
  max_stat <- apply(gene_module_matrix, 2, function(module) {
    apply(t(test_result) * module, 2, max)
  })


  combined_matrix <- combine_columns_from_matrices(list(power_sum_1,power_sum_2,power_sum_3,power_sum_4,power_sum_5,power_sum_6,power_sum_7, power_sum_8,max_stat), gene_module_matrix)
  # combined_matrix <- combine_columns_from_matrices(list(max_stat), gene_module_matrix)

  return(combined_matrix)
}


extract_and_combine_test_result_columns <- function(test_result, gene_module_matrix) {

  extracted_matrices <- list()

  # Loop through each column of gene_module_matrix
  for (j in seq_len(ncol(gene_module_matrix))) {
    # Find the indices of columns in test_result corresponding to "1" in gene_module_matrix[, j]
    selected_indices <- which(gene_module_matrix[, j] == 1)

    # Extract the corresponding columns from test_result
    extracted_matrix <- test_result[, selected_indices, drop = FALSE]

    # Set the column names to the current column name of gene_module_matrix
    colnames(extracted_matrix) <- rep(colnames(gene_module_matrix)[j], length(selected_indices))

    # Store the matrix in the list
    extracted_matrices[[j]] <- extracted_matrix
  }

  # Combine all extracted matrices into a single matrix using cbind
  combined_matrix <- do.call(cbind, extracted_matrices)

  return(combined_matrix)
}



combine_test_statistics_tensor <- function(test_tensor, gene_module_matrix) {
  # Combine all slices into a single matrix using rbind
  combined_matrix <- do.call(rbind, lapply(seq_len(dim(test_tensor)[3]), function(j) test_tensor[, , j]))

  # Apply combine_test_statistics() to the combined matrix
  combined_result <- combine_test_statistics(combined_matrix, gene_module_matrix)

  # Determine the new dimensions for the output tensor
  new_dims <- c(dim(test_tensor)[1], ncol(combined_result), dim(test_tensor)[3])

  # Initialize the output tensor
  combined_test_tensor <- array(NA, dim = new_dims)

  # Set dimnames for the output tensor
  dimnames(combined_test_tensor) <- list(
    dimnames(test_tensor)[[1]],  # Row names unchanged
    colnames(combined_result),  # Updated column names
    dimnames(test_tensor)[[3]]  # Slice names unchanged
  )

  # Transfer the combined result back to the tensor
  for (j in seq_len(dim(test_tensor)[3])) {
    combined_test_tensor[, , j] <- combined_result[
      (1 + (j - 1) * nrow(test_tensor)):(j * nrow(test_tensor)), ]
  }

  return(combined_test_tensor)
}


#################### Optional: Keep original gene_wise test statistics within gene module #################


combine_test_statistics_advanced <- function(test_result, gene_module_matrix) {
  # Replace NA values with 0 and take the absolute values
  test_result[is.na(test_result)] <- 0
  test_result <- abs(test_result)

  # Calculate power sums
  power_sum_1 <- test_result %*% gene_module_matrix
  power_sum_2 <- (test_result^2%*%gene_module_matrix)^(1/2)
  power_sum_3 <- (test_result^3%*%gene_module_matrix)^(1/3)
  power_sum_4 <- (test_result^4%*%gene_module_matrix)^(1/4)
  power_sum_5 <- (test_result^5%*%gene_module_matrix)^(1/5)
  power_sum_6 <- (test_result^6%*%gene_module_matrix)^(1/6)
  power_sum_7 <- (test_result^7%*%gene_module_matrix)^(1/7)
  power_sum_8 <- (test_result^8%*%gene_module_matrix)^(1/8)
  # Calculate max_stat
  max_stat <- apply(gene_module_matrix, 2, function(module) {
    apply(t(test_result) * module, 2, max)
  })

  # Extract and combine gene-wise results within modules
  gene_wise_within_module <- extract_and_combine_test_result_columns(test_result, gene_module_matrix)

  # Row-bind all results
  combined_matrix <- cbind(
    power_sum_1,
    power_sum_2,
    power_sum_3,
    power_sum_4,
    power_sum_5,
    power_sum_6,
    power_sum_7,
    power_sum_8,
    max_stat,
    gene_wise_within_module
  )

  # Order columns based on colnames in gene_module_matrix
  col_order <- colnames(gene_module_matrix)
  combined_matrix <- combined_matrix[, order(match(colnames(combined_matrix), col_order)), drop = FALSE]

  return(combined_matrix)
}

combine_test_statistics_tensor_advanced <- function(test_tensor, gene_module_matrix) {
  # Combine all slices into a single matrix using rbind
  combined_matrix <- do.call(rbind, lapply(seq_len(dim(test_tensor)[3]), function(j) test_tensor[, , j]))

  # Apply combine_test_statistics() to the combined matrix
  combined_result <- combine_test_statistics_advanced(combined_matrix, gene_module_matrix)

  # Determine the new dimensions for the output tensor
  new_dims <- c(dim(test_tensor)[1], ncol(combined_result), dim(test_tensor)[3])

  # Initialize the output tensor
  combined_test_tensor <- array(NA, dim = new_dims)

  # Set dimnames for the output tensor
  dimnames(combined_test_tensor) <- list(
    dimnames(test_tensor)[[1]],  # Row names unchanged
    colnames(combined_result),  # Updated column names
    dimnames(test_tensor)[[3]]  # Slice names unchanged
  )

  # Transfer the combined result back to the tensor
  for (j in seq_len(dim(test_tensor)[3])) {
    combined_test_tensor[, , j] <- combined_result[
      (1 + (j - 1) * nrow(test_tensor)):(j * nrow(test_tensor)), ]
  }

  return(combined_test_tensor)
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
  } else if (side_code ==3) # two side
  {
    rank <- 2*pmin(rank, (n_permute+2-rank))
  }

  rank[rank == (n_permute + 1)] <- n_permute

  dimnames(rank)<-dimnames(tensor)

  return(rank/(n_permute+1))
}

compute_p_value_from_tensor_precise <- function(tensor, n_permute, side_code) {
  dims<-dim(tensor)

  p_values <- array(NA, dim = dims, dimnames = dimnames(tensor))

  # Loop through each slice along the second and third dimensions
  for (i in seq_len(dims[2])) {
    for (j in seq_len(dims[3])) {
      # Extract the vector along the first dimension for tensor[, i, j]
      vec <- tensor[, i, j]

      # Count the number of non-NA values
      non_na_count <- sum(!is.na(vec))

      if (side_code==1){ #right tail
        ranks <- rank(-vec, ties.method = "average", na.last = "keep")

      } else if (side_code==2){
        ranks <- rank(vec, ties.method = "average", na.last = "keep")

      } else if (side_code==3){
        ranks <- 2 * pmin(rank(vec, ties.method = "average", na.last = "keep"),
                         non_na_count+1-rank(vec, ties.method = "average", na.last = "keep"))
      }

      else {
        stop("Unknown side code")
      }

      # Divide ranks by the number of non-NA values
      p_values[, i, j] <- ranks / non_na_count
    }
  }

  # Return the array of p-values
  return(p_values)
}


min_p_val <- function(p_val_tensor, gene_module_index_matrix) {
  # Get the dimensions of the input tensor
  dim_tensor <- dim(p_val_tensor)

  # Row-bind each slice of p_val_tensor
  combined_matrix <- do.call(rbind, lapply(seq_len(dim_tensor[3]), function(slice_idx) {
    p_val_tensor[, , slice_idx]
  }))

  # Find unique column names
  unique_colnames <- unique(colnames(combined_matrix))

  # Calculate element-wise minimum for columns with the same name
  reduced_matrix <- do.call(cbind, lapply(unique_colnames, function(col_name) {
    columns_with_same_name <- combined_matrix[, colnames(combined_matrix) == col_name, drop = FALSE]
    apply(columns_with_same_name, 1, function(row) {
      if (all(is.na(row))) {
        return(NA)  # If all values in the row are NA, return NA
      } else {
        return(min(row, na.rm = TRUE))  # Calculate the minimum ignoring NA
      }
    })  # Element-wise minimum
  }))

  # Assign the new column names
  colnames(reduced_matrix) <- unique_colnames

  # Recover the tensor
  min_p_tensor <- array(NA, dim = c(dim_tensor[1], ncol(gene_module_index_matrix)
, dim_tensor[3]),
                        dimnames = list(
                          dimnames(p_val_tensor)[[1]],
                          colnames(gene_module_index_matrix),
                          dimnames(p_val_tensor)[[3]]
                        ))

  # Split the reduced matrix back into slices
  for (slice_idx in seq_len(dim_tensor[3])) {
    start_row <- (slice_idx - 1) * dim_tensor[1] + 1
    end_row <- slice_idx * dim_tensor[1]
    min_p_tensor[, , slice_idx] <- reduced_matrix[start_row:end_row, ]
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
    final_p_vals[i, ] <- apply(current_slice, 2, function(col) rank(col, ties.method = "min", na.last="keep")[1] / sum(!is.na(col)))
  }

  return(final_p_vals)
}



resample_columns <- function(mat, n) {
  resampled_list <- lapply(1:n, function(i) {
    apply(mat, 2, sample)  # Permute each column
  })
  # Convert list of matrices into a single matrix with renamed columns
  resampled_matrix <- do.call(cbind, resampled_list)
  colnames(resampled_matrix) <- as.vector(sapply(colnames(mat), function(name) {
    paste0(name, "_", 1:n)
  }))
  return(resampled_matrix)
}


fit_and_evaluate_skew_t <- function(x,TStat1,side){

  x.mean = mean(x)
  x.sd = sd(x)
  x.norm = (x - x.mean)/x.sd
  mod = sgt.mle(X.f = ~x.norm,
                start = list(mu = 0, sigma = 1, lambda = 0, p = 2, q = 10))
  para = mod$estimate

  if(side==1){
  p = 1-psgt((TStat1-x.mean)/x.sd, mu = para['mu'], sigma = para['sigma'],
           lambda = para['lambda'], p = para['p'], q = para['q'])
  }
  return(p)
}
