setwd("/project/xuanyao/jiaming/Getting_started")

library(presto)
library(sceptre)
library(sceptredata)
library(readr)
library(dplyr)
library(Matrix)
library(parallel)
library(BH)
library(Rcpp)
library(peakRAM)
library(sgt)


sourceCpp("code/fit_skew_normal.cpp")
Morris_data <- T

setwd("/project/xuanyao/jiaming/Getting_started")

if (Morris_data) {
  #directories <- "data/STINGseq-v1_GDO"
  directories <- "data/GDO_threshold"
  grna_target_data_frame <- data.frame(read_csv(paste0(directories,"/grna_target_data_frame.csv")))
  grna_target_data_frame_omit  <- grna_target_data_frame %>%
    filter(!(grna_id %in% c("NTC-3", "NTC-12")))

    sceptre_object <- import_data_from_cellranger(
    directories = directories,
    moi = "high",
    grna_target_data_frame = grna_target_data_frame_omit
  )
} else {
  directories <- paste0(system.file("extdata", package = "sceptredata"), "/highmoi_example/gem_group_", 1:2)
  data(grna_target_data_frame_highmoi)
  sceptre_object <- import_data_from_cellranger(
    directories = directories,
    moi = "high",
    grna_target_data_frame = grna_target_data_frame_highmoi
  )
}


positive_control_pairs <- construct_positive_control_pairs(sceptre_object)

discovery_pairs_cis <- construct_cis_pairs(
  sceptre_object = sceptre_object,
  positive_control_pairs = positive_control_pairs,
  response_position_data_frame=gene_position_data_frame_grch37,
  distance_threshold = 5e5
)

discovery_pairs_trans <- construct_trans_pairs(
  sceptre_object = sceptre_object,
  positive_control_pairs = positive_control_pairs,
  pairs_to_exclude = "pairs_containing_pc_targets"
)
# sceptre_object <- sceptre_object |> # |> is R's base pipe, similar to %>%
#   set_analysis_parameters(
#     discovery_pairs=discovery_pairs_trans,
#     positive_control_pairs=positive_control_pairs,
#     side="both") |>
#   assign_grnas()

sceptre_object <- sceptre_object |>  # |> is R's base pipe, similar to %>%
  set_analysis_parameters(
    discovery_pairs = discovery_pairs_trans,
    positive_control_pairs = positive_control_pairs,
    side = "both",
    grna_integration_strategy="singleton"
  ) |>
  assign_grnas(
    method = "thresholding",
    threshold = 1
  )


compute_ustat <- function(Xr, cols_list, n1n2_list, group_size_list) {
  total_loops <- length(cols_list) # Total number of loops
  loop_counter <- 1 # Initialize loop counter

  mapply(function(cols, n1n2, group_size) {
    cat("Compute wilcox test statisitcs for gRNA #", loop_counter, "of", total_loops,".    ", round(loop_counter / total_loops * 100, 2)
        ,"% of total progress" , "\n")
    loop_counter <<- loop_counter + 1 # Increment the counter

    grs <- sumGroups(Xr, cols)

    if (is(Xr, "dgCMatrix")) {
      gnz <- (group_size - nnzeroGroups(Xr, cols))
      zero_ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
      ustat <- t((t(gnz) * zero_ranks)) + grs - group_size * (group_size + 1) / 2
    } else {
      ustat <- grs - group_size * (group_size + 1) / 2
    }

    return(ustat)
  }, cols_list, n1n2_list, group_size_list, SIMPLIFY = FALSE)
}



compute_pval <- function(ustat_list, ties, N, n1n2_list) {
  total_loops <- length(ustat_list) # Total number of loops
  loop_counter <- 1 # Initialize loop counter
  mapply(function(ustat, n1n2) {
    cat("Compute wilcox p value for gRNA #", loop_counter, "of", total_loops,".    ", round(loop_counter / total_loops * 100, 2)
        ,"% of total progress" , "\n")
    loop_counter <<- loop_counter + 1 # Increment the counter
    z <- ustat - 0.5 * n1n2
    z <- z - sign(z) * 0.5
    .x1 <- N ^ 3 - N
    .x2 <- 1 / (12 * (N^2 - N))

    rhs <- lapply(ties, function(tvals) {
      (.x1 - sum(tvals ^ 3 - tvals)) * .x2
    }) %>% unlist

    usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
    z <- t(z / usigma)

    pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    return(pvals)
  }, ustat_list, n1n2_list, SIMPLIFY = FALSE)
}



tidy_results_reduced <- function(wide_res, features, groups) {
  res <- Reduce(cbind, lapply(wide_res, as.numeric)) %>% data.frame()
  colnames(res) <- names(wide_res)

  # Create the 'feature' and 'group' columns
  res$feature <- rep(features, times = length(groups))
  res$group <- rep(groups, each = length(features))

  # Filter only the first group
  first_group <- groups[1]
  res <- res %>% dplyr::filter(group == first_group)

  # Select the relevant columns
  res %>% dplyr::select(
    feature,
    statistic,
    pval
  )
}


wilcoxauc <- function(X, ...) {
  UseMethod("wilcoxauc")
}

wilcoxauc.default <- function(X, y, groups_use = NULL, verbose = TRUE, ...) {
  ## Check and possibly correct input values
  if (is(X, "dgeMatrix")) X <- as.matrix(X)
  if (is(X, "data.frame")) X <- as.matrix(X)
  if (is(X, "dgTMatrix")) X <- as(X, "dgCMatrix")
  if (is(X, "TsparseMatrix")) X <- as(X, "dgCMatrix")
  if (ncol(X) != length(y[[1]])) stop("number of columns of X does not
                                match length of y")
  if (!is.null(groups_use)) {
    idx_use <- which(y %in% intersect(groups_use, y))
    y <- y[idx_use]
    X <- X[, idx_use]
  }

  y  <- lapply(y, factor)


  group.size <- lapply(y, function(vec) as.numeric(table(vec)))


  if (is.null(row.names(X))) {
    row.names(X) <- paste0("Feature", seq_len(nrow(X)))
  }

  ## Compute primary statistics

  n1n2 <- lapply(group.size, function(gs) gs * (ncol(X) - gs))

  if (is(X, "dgCMatrix")) {
    rank_res <- rank_matrix(Matrix::t(X))
  } else {
    rank_res <- rank_matrix(X)
  }

  ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
  pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2)



  results_list <- lapply(seq_along(ustat), function(i) {
    # Create a res_list for each i-th vector in ustat and pvals
    res_list <- list(
      pval = pvals[[i]],
      statistic = t(ustat[[i]])
    )

    tidy_results_reduced(res_list, row.names(X), levels(y[[i]]))
  })

  # results_list is now a list of the tidy_results_reduced outputs
  return(results_list)
}
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
slice_to_list <- function(large_matrix, n_permute) {
  # Initialize an empty list to store the smaller matrices
  sliced_matrices <- list()

  # Calculate the number of rows for each smaller matrix
  slice_size <- n_rows

  for (i in 1:(n_permute)) {
    start_row <- (i - 1) * slice_size + 1
    end_row <- i * slice_size

    sliced_matrices[[i]] <- large_matrix[start_row:end_row, ]
  }

  return(sliced_matrices)
}

# t_test <- function(response_matrix, grna_assignment_matrix){
#
#   n_cells <- ncol(response_matrix)
#   n_genes <- nrow(response_matrix)
#   all_cells <- seq_len(n_cells)
#   n_rows <- nrow(grna_assignment_matrix)
#
#   grna_infected_cell_count <- grna_assignment_matrix%*%matrix(rep(1,n_cells),nrow=n_cells)
#
#   mean_expression_trt <- grna_assignment_matrix%*%t(response_matrix)/as.vector(grna_infected_cell_count)
#
#   mean_expression_squared_trt <- grna_assignment_matrix%*%t(response_matrix^2)/as.vector(grna_infected_cell_count)
#
#   var_expression_trt <- mean_expression_squared_trt - mean_expression_trt^2
#
#   grna_control_cell_count <- n_cells-grna_infected_cell_count
#
#   mean_expression_control <- (1-grna_assignment_matrix)%*%t(response_matrix)/as.vector(grna_control_cell_count)
#
#   mean_expression_squared_control <-(1-grna_assignment_matrix)%*%t(response_matrix^2)/as.vector(grna_control_cell_count)
#
#   var_expression_control <- mean_expression_squared_control - mean_expression_control^2
#
#   t_statistics <- (mean_expression_trt-mean_expression_control)/sqrt(var_expression_trt/as.vector(grna_infected_cell_count)+var_expression_control/as.vector(grna_control_cell_count))
#
#   return(t_statistics)
# }

t_test <- function(response_matrix, grna_assignment_matrix) {

  response_matrix <- as.matrix(response_matrix)
  grna_assignment_matrix <- as.matrix(grna_assignment_matrix)

  n_cells <- ncol(response_matrix)
  n_genes <- nrow(response_matrix)
  all_cells <- seq_len(n_cells)
  n_rows <- nrow(grna_assignment_matrix)

  response_matrix_t <- t(response_matrix)
  response_matrix_t_squared <-response_matrix_t^2

  grna_infected_cell_count <- as.vector(grna_assignment_matrix %*% matrix(rep(1, n_cells), nrow = n_cells))

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

wilcox_test <- function(response_matrix, grna_assignment_matrix){

  n_rows <- nrow(grna_assignment_matrix)

  response_dgC <- as(as.matrix(response_matrix),"dgCMatrix")

  grna_assignment_list <- split(grna_assignment_matrix, row(grna_assignment_matrix))

  wilcox_result <- wilcoxauc(response_dgC,grna_assignment_list)

  names(wilcox_result) <- rownames(grna_assignment_matrix)

  statistic_list <- lapply(wilcox_result, function(df) df$pval)


  statistic_matrix <- do.call(rbind, statistic_list)
  colnames(statistic_matrix) <- rownames(response_matrix)

  return(statistic_matrix)
}
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

  power_sum_1 <- abs(test_result)%*%gene_module_matrix
  power_sum_2 <- abs(test_result)^2%*%gene_module_matrix
  power_sum_3 <- abs(test_result)^3%*%gene_module_matrix
  power_sum_4 <- abs(test_result)^4%*%gene_module_matrix
  power_sum_5 <- abs(test_result)^5%*%gene_module_matrix
  power_sum_6 <- abs(test_result)^6%*%gene_module_matrix

  combined_matrix <- combine_columns_from_matrices(list(power_sum_1,power_sum_2,power_sum_3,power_sum_4,power_sum_5,power_sum_6), gene_module_matrix)

  return(combined_matrix)

}

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
  return(rank/(n_permute+1))
}

compute_p_value_from_skew_normal_2 <- function(tensor, side_code) {
  # Get dimensions and dimnames of the tensor
  dim_tensor <- dim(tensor)
  dimnames_tensor <- dimnames(tensor)

  # Initialize a tensor to store p-values with the same dimensions and dimnames
  p_value_tensor <- array(NA, dim = c(dim_tensor[1], dim_tensor[2], dim_tensor[3]), dimnames = dimnames_tensor)

  # Apply over each slice (i.e., the 3rd dimension)
  result <- apply(tensor, 3, function(slice) {
    # Parallelize over the second dimension (columns) for each slice
    mclapply(1:dim_tensor[2], function(i) {
      # Apply the function to the column vector slice[,i]
      fit_and_evaluate_skew_normal(slice[,i], side_code, F)
    }, mc.cores = parallel::detectCores() - 1)  # Adjust the number of cores used
  })

  # Rearrange the results back into the p_value_tensor
  for (j in 1:dim_tensor[3]) {
    for (i in 1:dim_tensor[2]) {
      p_value_tensor[, i, j] <- result[[j]][[i]]
    }
  }

  return(p_value_tensor)
}

compute_p_value_from_skew_normal <- function(tensor, side_code) {
  # Get dimensions and dimnames of the tensor
  dim_tensor <- dim(tensor)
  dimnames_tensor <- dimnames(tensor)

  # Initialize an empty tensor with the same dimensions and dimnames
  p_value_tensor <- array(NA, dim = dim_tensor, dimnames = dimnames_tensor)

  # Use mclapply to apply the function in parallel across all slices
  result <- mclapply(1:dim_tensor[3], function(j) {
    # Apply the function across the first dimension (i.e., columns of tensor[,i,j])
    apply(tensor[,,j], 2, function(column_vector) {
      tryCatch({
        # Apply the fit_and_evaluate_skew_normal function to each column vector
        fit_and_evaluate_skew_normal(column_vector, side_code, F)
      }, error = function(e) {
        # Catch and handle any errors, print message, return NA
        message("Error in slice ", j, " for column vector: ", e)
        return(rep(NA, length(column_vector)))  # Return NA for the column if an error occurs
      })
    })
  }, mc.cores = parallel::detectCores() - 1)  # Adjust the number of cores used

  # Combine the results back into the tensor
  for (j in 1:dim_tensor[3]) {
    p_value_tensor[,,j] <- result[[j]]
  }

  return(p_value_tensor)
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

  # Use mclapply to parallelize the processing over the 3rd dimension (slices)
  results <- mclapply(1:dim_tensor[3], function(i) {
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
  }, mc.cores = parallel::detectCores() - 1)  # Adjust the number of cores used

  # Combine the results from mclapply into the min_p_tensor
  for (i in 1:dim_tensor[3]) {
    min_p_tensor[,,i] <- results[[i]]
  }

  return(min_p_tensor)
}
compute_final_p_val <- function(min_p_tensor, n_permute) {
  # Get the dimensions of the tensor
  dim_tensor <- dim(min_p_tensor)

  # Initialize an empty matrix to store the p-value counts
  final_p_vals <- matrix(0, nrow = dim_tensor[3], ncol = dim_tensor[2],
                         dimnames = list(dimnames(min_p_tensor)[[3]], dimnames(min_p_tensor)[[2]]))

  # Use mclapply to parallelize over each slice of the tensor (rows)
  results <- mclapply(1:dim_tensor[3], function(i) {
    # Extract the current slice (2D matrix)
    current_slice <- min_p_tensor[,,i]

    # Extract the original values (the first row)
    original_values <- current_slice[1,]

    # Use matrix broadcasting and colSums to count the number of permutation values smaller than original
    colSums(current_slice <= matrix(original_values, nrow = nrow(current_slice),
                                    ncol = ncol(current_slice), byrow = TRUE))
  }, mc.cores = parallel::detectCores() - 1)  # Adjust based on available cores

  # Assign the results to the final_p_vals matrix
  for (i in 1:dim_tensor[3]) {
    final_p_vals[i, ] <- results[[i]]
  }

  return(final_p_vals / (n_permute + 1))
}

compute_final_p_val_skew_normal <- function(tensor) {
  # Get the dimensions of the tensor
  dim_tensor <- dim(tensor)

  p_value_matrix <- matrix(NA, nrow = dim_tensor[3], ncol = dim_tensor[2])

  # Apply mclapply in parallel over the 3rd dimension (slices)
  result <- mclapply(1:dim_tensor[3], function(j) {
    # Apply fit_and_evaluate_skew_normal to each column vector in the 2nd dimension for the current slice
    sapply(1:dim_tensor[2], function(i) {
      fit_and_evaluate_skew_normal(tensor[,i,j], 2, T)
    })
  }, mc.cores = parallel::detectCores() - 1)  # Use all but 1 core for parallel processing

  # Combine results into p_value_matrix
  for (j in 1:dim_tensor[3]) {
    p_value_matrix[j,] <- result[[j]]
  }
  rownames(p_value_matrix)<-dimnames(tensor)[[3]]
  colnames(p_value_matrix)<-dimnames(tensor)[[2]]

  return(p_value_matrix)
}
fit_skew_t <- function(x,TStat1){
  x.mean = mean(x)
  x.sd = sd(x)
  x.norm = (x - x.mean)/x.sd
  mod = sgt.mle(X.f = ~x.norm,
                start = list(mu = 0, sigma = 1, lambda = 0, p = 2, q = 10))
  para = mod$estimate
  p = psgt((TStat1-x.mean)/x.sd, mu = para['mu'], sigma = para['sigma'],
           lambda = para['lambda'], p = para['p'], q = para['q'])
  return(p)
}

compute_final_p_val_skew_normal_2 <- function(tensor) {

  dim_tensor <- dim(tensor)
  p_value_matrix <- matrix(NA, nrow = dim_tensor[3], ncol = dim_tensor[2])

  # Apply over the third dimension
  result <- lapply(1:dim_tensor[3], function(j) {
    sapply(1:dim_tensor[2], function(i) {
      fit_skew_t(tensor[-1, i, j],tensor[1,i,j])  # Apply function across the first dimension for each (i,j)
    })
  })

  # Combine results into p_value_matrix
  for (j in 1:dim_tensor[3]) {
    p_value_matrix[j, ] <- result[[j]]
  }

  rownames(p_value_matrix) <- dimnames(tensor)[[3]]
  colnames(p_value_matrix) <- dimnames(tensor)[[2]]

  return(p_value_matrix)
}

analysis_gene_module <- function (grna_assignment_matrix, response_matrix, gene_module_index_matrix, test_type, n_permute, use_resample, use_approximation) {

  func_start_time <- Sys.time()

  # permuted grna assignment matrix
  grna_assignment_matrix_combined <- permutation_or_resampling(grna_assignment_matrix, n_permute, use_resample)

  # Test statistics calculation
  step_start_time <- Sys.time()
  if (test_type == "t_test") {
    test_result <- t_test(response_matrix, grna_assignment_matrix_combined)
    step_end_time <- Sys.time()
    print(paste("Time for test statistics calculation:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))

    # Combined statistics
    step_start_time <- Sys.time()
    combined_result <- combine_test_statistics(test_result, gene_module_index_matrix)
    step_end_time <- Sys.time()
    print(paste("Time for combine test statistics:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))

    # Change matrix to tensor
    step_start_time <- Sys.time()
    test_tensor <- convert_matrix_to_tensor(combined_result, n_permute, grna_assignment_matrix)
    step_end_time <- Sys.time()
    print(paste("Time for matrix to tensor conversion:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))

    # Calculate p-value
    step_start_time <- Sys.time()
    p_tensor <- compute_p_value_from_tensor(test_tensor, n_permute, 1)
    step_end_time <- Sys.time()
    print(paste("Time for p-value calculation:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))


  } else if (test_type == "wilcox_test") {
    test_result <- wilcox_test(response_matrix, grna_assignment_matrix_combined)
    step_end_time <- Sys.time()
    print(paste("Time for test statistics calculation:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))

    # Combined statistics
    step_start_time <- Sys.time()
    combined_result <- combine_test_statistics(test_result, gene_module_index_matrix)
    step_end_time <- Sys.time()
    print(paste("Time for combine test statistics:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))

    # Change matrix to tensor
    step_start_time <- Sys.time()
    test_tensor <- convert_matrix_to_tensor(combined_result, n_permute, grna_assignment_matrix)
    step_end_time <- Sys.time()
    print(paste("Time for matrix to tensor conversion:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))

    # Calculate p-value
    step_start_time <- Sys.time()
    p_tensor <- compute_p_value_from_tensor(test_tensor, n_permute, 2)
    step_end_time <- Sys.time()
    print(paste("Time for p-value calculation:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))
  } else {
    stop("Unknown test type")
  }

  min_pval <- min_p_val(p_tensor, gene_module_index_matrix)

  precompute_p <- compute_final_p_val(min_pval, n_permute)

  if (use_approximation) {

    n_test <- dim(precompute_p)[1]
    threshold <- 0.1/n_test + sqrt(0.1/n_test/n_permute)
    indices <- which(precompute_p <= threshold, arr.ind = TRUE)


    # Use apply to loop through indices
    apply(indices, 1, function(idx) {
      i <- idx[2]
      j <- idx[1]

      # Perform the fit_skew_t calculation
      result <- fit_skew_t(min_pval[-1, i, j], min_pval[1, i, j])

      # Update final_p inside the function
      precompute_p[j, i] <<- result
    })
  }

  func_end_time <- Sys.time()
  print(paste("Total epoch execution time:", as.numeric(difftime(func_end_time, func_start_time, units = "secs"))))

  return(precompute_p)
}

analyze_chunk <- function(chunk,response_matrix) {
  analysis_gene_module(
    chunk,
    response_matrix,
    gene_module_index_matrix,
    test_type,
    n_permute,
    use_resample,
    use_approximation
  )
}


setwd("/project/xuanyao/jiaming/Getting_started")
gene_names<-rownames(sceptre_object@response_matrix[[1]])

gene_modules <- readRDS("data/gene_modules_id.rds")
gene_in_modules <- vector("logical", length(gene_names))
covariate <- sceptre_object@covariate_data_frame


for (i in seq_along(gene_names)) {
  gene <- gene_names[i]
  gene_in_modules[i] <- any(sapply(gene_modules, function(module) gene %in% module))
}

response_matrix_raw<-sceptre_object@response_matrix[[1]][gene_in_modules,]


# response_matrix_normalized<-t(t(response_matrix_raw)/covariate$response_n_umis)
response_matrix_normalized <- response_matrix_raw

row_sums <- rowSums(response_matrix_normalized!=0)

gene_to_exclude <- rownames(response_matrix_normalized)[row_sums <= 0.05*ncol(response_matrix_normalized)]

response_matrix <- response_matrix_normalized[!rownames(response_matrix_normalized) %in% gene_to_exclude, ]

gene_module_index_matrix <- convert_gene_list_to_matrix(response_matrix, gene_modules)

grna_assignment_matrix <-get_grna_assignments(
  sceptre_object = sceptre_object
)
n_rows <- nrow(grna_assignment_matrix)


use_resample <- T
use_approximation <- T
rank <- T
test_type="t_test"
n_permute <- 2000
chunk_size <- 50

if (rank) {
  response_matrix <-  t(rank_matrix(as.matrix(response_matrix))$X_ranked)
}

# grna_assignment_matrix <-get_grna_assignments(
#   sceptre_object = sceptre_object
# )[177:200,]



# peak_mem_usage <- peakRAM(
#   result <- analysis_gene_module(grna_assignment_matrix, response_matrix,
#                                  gene_module_index_matrix,test_type, n_permute, use_resample, use_approximation)
#   )

# peak_mem_usage <- peakRAM(
#   result <- analysis_gene_module(grna_assignment_matrix, t(rank_matrix(as.matrix(response_matrix))$X_ranked),
#                                  gene_module_index_matrix,test_type, n_permute, use_resample, use_approximation)
# )






chunks <- split(1:n_rows, rep(1:ceiling(n_rows / chunk_size), each = chunk_size, length.out = n_rows))

peak_mem_usage <- peakRAM(

  result <- do.call(rbind, lapply(chunks, function(rows) {
    chunk_matrix <- grna_assignment_matrix[rows, , drop = FALSE]
    analyze_chunk(chunk_matrix,response_matrix)
  }))

)

print(peak_mem_usage)

write.csv(result, "output/t_test_on_whole_dataset_author_threshold.csv", row.names = TRUE)
