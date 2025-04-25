#' Filter miRNAs based on expression in a minimum percentage of samples in any group
#'
#' This function filters miRNAs that have reads in at least a specified percentage of samples in any group.
#' The comparison operator for both the threshold and the read count can be customized (>, >=, <, <=).
#'
#' @param count_matrix A count matrix (miRNAs in rows, samples in columns).
#' @param metadata A data frame containing sample information, including a 'Condition' column with the group for each sample.
#' @param threshold Minimum percentage of samples with reads (value between 0 and 1, default is 0.5).
#' @param min_reads Minimum number of reads to consider a miRNA expressed in a sample (default is 1000).
#' @param threshold_comparison Comparison operator for the threshold. Options are ">", ">=", "<", or "<=" (default is ">=").
#' @param read_comparison Comparison operator for the read count. Options are ">", ">=", "<", or "<=" (default is ">").
#' @return A filtered count matrix containing miRNAs that meet the criteria in at least one group.
#' @export
#' @examples
#' # Example usage:
#' # count_matrix <- matrix(c(1000, 0, 2000, 1500, 0, 3000), nrow = 2, dimnames = list(c("miRNA1", "miRNA2"), c("Sample1", "Sample2", "Sample3")))
#' # metadata <- data.frame(Sample = c("Sample1", "Sample2", "Sample3"), Condition = c("A", "A", "B"))
#' # filtered_matrix <- filter_mirnas(count_matrix, metadata, threshold = 0.5, min_reads = 1000, threshold_comparison = ">=", read_comparison = ">")
filter_mirnas <- function(count_matrix, metadata, threshold = 0.5, min_reads = 1000, threshold_comparison = ">=", read_comparison = ">") {
  # Check if the number of samples in count_matrix matches the number of rows in metadata
  if (ncol(count_matrix) != nrow(metadata)) {
    stop("The number of columns in count_matrix must match the number of rows in metadata.")
  }

  # Check if the comparison operators are valid
  valid_comparisons <- c(">", ">=", "<", "<=")
  if (!(threshold_comparison %in% valid_comparisons)) {
    stop("Invalid threshold comparison operator. Use '>', '>=', '<', or '<='.")
  }
  if (!(read_comparison %in% valid_comparisons)) {
    stop("Invalid read comparison operator. Use '>', '>=', '<', or '<='.")
  }

  # Create a vector with the group names corresponding to each sample
  group_vector <- metadata$Condition

  # Get the unique group names
  groups <- unique(group_vector)

  # Initialize a vector to store the miRNAs that meet the condition
  filtered_miRNAs <- c()

  # Iterate over each miRNA
  for (i in 1:nrow(count_matrix)) {
    miRNA <- rownames(count_matrix)[i]
    passes_filter <- FALSE  # Variable to check if the miRNA passes the filter in at least one group

    for (group in groups) {
      group_samples <- which(group_vector == group)

      # Calculate the minimum number of samples with reads required
      min_samples <- ceiling(length(group_samples) * threshold)

      # Count how many samples have reads for this miRNA based on the read comparison operator
      if (read_comparison == ">") {
        samples_with_reads <- sum(count_matrix[i, group_samples] > min_reads)
      } else if (read_comparison == ">=") {
        samples_with_reads <- sum(count_matrix[i, group_samples] >= min_reads)
      } else if (read_comparison == "<") {
        samples_with_reads <- sum(count_matrix[i, group_samples] < min_reads)
      } else if (read_comparison == "<=") {
        samples_with_reads <- sum(count_matrix[i, group_samples] <= min_reads)
      }

      # Apply the threshold comparison operator
      if (threshold_comparison == ">") {
        passes_group <- samples_with_reads > min_samples
      } else if (threshold_comparison == ">=") {
        passes_group <- samples_with_reads >= min_samples
      } else if (threshold_comparison == "<") {
        passes_group <- samples_with_reads < min_samples
      } else if (threshold_comparison == "<=") {
        passes_group <- samples_with_reads <= min_samples
      }

      # If the miRNA passes the filter in at least one group, accept it and break the loop
      if (passes_group) {
        passes_filter <- TRUE
        break  # No need to check other groups
      }
    }

    # If the miRNA passes the filter in at least one group, add it to the list
    if (passes_filter) {
      filtered_miRNAs <- c(filtered_miRNAs, miRNA)
    }
  }

  # Filter the count matrix to include only the filtered miRNAs
  filtered_matrix <- count_matrix[rownames(count_matrix) %in% filtered_miRNAs, , drop = FALSE]

  return(filtered_matrix)
}
