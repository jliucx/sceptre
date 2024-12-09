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
source("code/functions.R")



author_threshold <- T

if (author_threshold) {
  #directories <- "data/STINGseq-v1_GDO"
  directories <- "data/gene_module_author_threshold"
  grna_target_data_frame <- data.frame(read_csv(paste0(directories,"/grna_target_data_frame.csv")))
  grna_target_data_frame_omit  <- grna_target_data_frame %>%
    filter(!(grna_id %in% c("NTC-3", "NTC-12")))

  sceptre_object <- import_data_from_cellranger(
    directories = directories,
    moi = "high",
    grna_target_data_frame = grna_target_data_frame_omit
  )

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

} else{

  directories <- "data/gene_module"
  grna_target_data_frame <- data.frame(read_csv(paste0(directories,"/grna_target_data_frame.csv")))
  sceptre_object <- import_data_from_cellranger(
    directories = directories,
    moi = "high",
    grna_target_data_frame = grna_target_data_frame
  )
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
  sceptre_object <- sceptre_object |> # |> is R's base pipe, similar to %>%
    set_analysis_parameters(
      discovery_pairs=discovery_pairs_trans,
      positive_control_pairs=positive_control_pairs,
      side="both") |>
    assign_grnas()

}


analysis_gene_module <- function (grna_assignment_matrix, response_matrix, gene_module_index_matrix, test_type, n_permute, use_resample, use_approximation) {
  set.seed(0)
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
    p_tensor <- compute_p_value_from_tensor(test_tensor, n_permute, 1)
    step_end_time <- Sys.time()
    print(paste("Time for p-value calculation:", as.numeric(difftime(step_end_time, step_start_time, units = "secs"))))
  } else {
    stop("Unknown test type")
  }

  step_start_time <- Sys.time()

  min_pval <- min_p_val(p_tensor, gene_module_index_matrix)

  precompute_p <- compute_final_p_val(min_pval, n_permute)

  if (use_approximation) {
    threshold <- 0.1

    # Find indices of elements in precompute_p below the threshold
    indices <- which(precompute_p < threshold, arr.ind = TRUE)

    if (length(indices) > 0) {
      for (k in seq_len(nrow(indices))) {
        # Extract row (j) and column (i) indices
        i <- indices[k, 2]
        j <- indices[k, 1]

        # Calculate p-values
        p_value_1 <- qnorm(min_pval[-1, i, j])
        p_value_2 <- qnorm(min_pval[1, i, j])



        if (all(is.finite(p_value_1)) && is.finite(p_value_2)) {
          # Fit and evaluate the skew-normal model
          result <- fit_and_evaluate_skew_normal(
            p_value_1,
            p_value_2,
            2,
            TRUE
          )

          # Update precompute_p with the result
          precompute_p[j, i] <- result

        }
      }
    } else {
      print("No indices meet the condition precompute_p < threshold")
    }
  }
  output_dir <- paste0("output/startover/", key)
  output_file <- paste0(output_dir, "/our_method.csv")
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

  if (!file.exists(output_file)) {
    write.table(precompute_p, file = output_file, sep = ",", row.names = TRUE, col.names = TRUE)
  }
  else {
    write.table(precompute_p, file = output_file, sep = ",", row.names = TRUE, col.names = FALSE, append = TRUE)
  }



  ############################ print gene wise p value ###################################
  if (T) {



    gene_wise_test_tensor <- convert_matrix_to_tensor(test_result, n_permute, grna_assignment_matrix)
    if (!use_approximation){
      gene_wise_p_tensor <- compute_p_value_from_tensor(gene_wise_test_tensor, n_permute, 1)
    } else {
      gene_wise_p_tensor <- array(NA, dim = dim(gene_wise_test_tensor),dimnames = dimnames(gene_wise_test_tensor))

      # Loop through all i and j
      dim1 <- dim(gene_wise_test_tensor)[1]
      dim2 <- dim(gene_wise_test_tensor)[2]
      dim3 <- dim(gene_wise_test_tensor)[3]

      for (i in seq_len(dim2)) {
        for (j in seq_len(dim3)) {
          # Exclude the first row for the first argument
          # Use the first row for the second argument
          gene_wise_p_tensor[, i, j] <- fit_and_evaluate_skew_normal(
            gene_wise_test_tensor[-1, i, j],
            gene_wise_test_tensor[1, i, j],
            1,
            FALSE
          )
        }
      }
    }

    gene_wise_p <- t(gene_wise_p_tensor[1,,])
    colnames(gene_wise_p) <- rownames(response_matrix)



    output_file <- paste0(output_dir, "/gene_wise_p.csv")

    # Create the directory if it doesn't exist

    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

    if (!file.exists(output_file)) {
      write.table(gene_wise_p, file = output_file, sep = ",", row.names = TRUE, col.names = TRUE)
      }
    else {
      write.table(gene_wise_p, file = output_file, sep = ",", row.names = TRUE, col.names = FALSE, append = TRUE)
      }


    ######################## print gene module minp value #####################################
    minp_list <- lapply(1:dim(gene_wise_p_tensor)[3], function(slice_idx) {
      apply(gene_module_index_matrix, 2, function(module) {

        module_genes <- which(module > 0)

        apply(gene_wise_p_tensor[,,slice_idx], 1, function(row) {
          if (length(module_genes) == 0) {
            return(NA)
          }
          return(min(row[module_genes], na.rm = TRUE))
        })
      })
    })

    minp <- array(unlist(minp_list), dim = c(dim(gene_wise_p_tensor)[1], ncol(gene_module_index_matrix), dim(gene_wise_p_tensor)[3]))
    final_minp <- compute_final_p_val(minp, n_permute)

    if (use_approximation) {
      threshold <- 0.1

      # Find indices of elements in precompute_p below the threshold
      indices <- which(final_minp < threshold, arr.ind = TRUE)

      if (length(indices) > 0) {
        for (k in seq_len(nrow(indices))) {
          # Extract row (j) and column (i) indices
          i <- indices[k, 2]
          j <- indices[k, 1]

          # Calculate p-values
          p_value_1 <- qnorm(minp[-1, i, j])
          p_value_2 <- qnorm(minp[1, i, j])

          if (all(is.finite(p_value_1)) && is.finite(p_value_2)) {
            # Fit and evaluate the skew-normal model
            result <- fit_and_evaluate_skew_normal(
              p_value_1,
              p_value_2,
              2,
              TRUE
            )

            # Update precompute_p with the result
            final_minp[j, i] <- result

          }
        }
      } else {
        print("No indices meet the condition precompute_p < threshold")
      }
    }

    rownames(final_minp) <- rownames(gene_wise_p)
    colnames(final_minp) <- colnames(gene_module_index_matrix)

    output_file <- paste0(output_dir, "/minp.csv")

    if (!file.exists(output_file)) {
      # Write header if the file does not exist
      write.table(final_minp, file = output_file, sep = ",", row.names = TRUE, col.names = TRUE)
    } else {
      # Append to the file if it already exists
      write.table(final_minp, file = output_file, sep = ",", row.names = TRUE, col.names = FALSE, append = TRUE)
    }

  }

  func_end_time <- Sys.time()
  print(paste("Total epoch execution time:", as.numeric(difftime(func_end_time, func_start_time, units = "secs"))))

  return(0)
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


gene_modules <- readRDS("data/gene_modules_id.rds")
covariate <- sceptre_object@covariate_data_frame


response_matrix_raw<-sceptre_object@response_matrix[[1]]
umis_quantiles <- quantile(covariate$response_n_umis, probs = c(0.01, 0.99))
nonzero_quantiles <- quantile(covariate$response_n_nonzero, probs = c(0.01, 0.99))

out_of_range_indices <- which(
  covariate$response_n_umis < umis_quantiles[1] |
    covariate$response_n_umis > umis_quantiles[2] |
    covariate$response_n_nonzero < nonzero_quantiles[1] |
    covariate$response_n_nonzero > nonzero_quantiles[2]
)

# response_matrix_raw <- response_matrix_raw[,-out_of_range_indices]
# covariate <- covariate[-out_of_range_indices,]

response_matrix_normalized<-t(t(response_matrix_raw)/covariate$response_n_umis)
# response_matrix_normalized <- response_matrix_raw


response_matrix <- response_matrix_normalized

gene_module_index_matrix <- convert_gene_list_to_matrix(response_matrix, gene_modules)

grna_assignment_matrix <-get_grna_assignments(
  sceptre_object = sceptre_object
)[151:208,]


n_rows <- nrow(grna_assignment_matrix)

key <- "no_approximation_20000"
use_resample <- T
use_approximation <- F
rank <- T
test_type="t_test"
n_permute <- 20000
chunk_size <- 50


if (rank) {
  response_matrix <-  t(rank_matrix(as.matrix(response_matrix))$X_ranked)
}
rownames(response_matrix)<- readLines("docs/row_names.txt")

# grna_assignment_matrix <-get_grna_assignments(
#   sceptre_object = sceptre_object
# )[177:200,]



# peak_mem_usage <- peakRAM(
#   result <- analysis_gene_module(grna_assignment_matrix, response_matrix,
#                                  gene_module_index_matrix,test_type, n_permute, use_resample, use_approximation)
#   )


chunks <- split(1:n_rows, rep(1:ceiling(n_rows / chunk_size), each = chunk_size, length.out = n_rows))

peak_mem_usage <- peakRAM({
  for (rows in chunks) {
    chunk_matrix <- grna_assignment_matrix[rows, , drop = FALSE]
    analyze_chunk(chunk_matrix, response_matrix)
  }
})
print(peak_mem_usage)
