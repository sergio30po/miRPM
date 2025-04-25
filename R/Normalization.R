#' Normalize a count matrix to RPM without scientific notation
#'
#' This function normalizes a count matrix to Reads Per Million (RPM) using the total reads per sample.
#' The resulting matrix is rounded to 2 decimal places and avoids scientific notation.
#'
#' @param count_matrix A count matrix (genes in rows, samples in columns).
#' @param metrics A dataframe containing sample information, including total reads.
#' @param reads_column The name of the column in metrics that contains the total reads.
#' @return A normalized matrix in RPM without scientific notation.
#' @export
#' @examples
#' # Example usage:
#' # count_matrix <- matrix(c(100, 200, 300, 400), nrow = 2, dimnames = list(c("Gene1", "Gene2"), c("Sample1", "Sample2")))
#' # metrics <- data.frame(Sample = c("Sample1", "Sample2"), TotalReads = c(1000, 2000))
#' # rpm_matrix <- normalize_rpm(count_matrix, metrics, "TotalReads")
normalize_rpm <- function(count_matrix, metrics, reads_column) {

  # Check if the reads column exists in metrics
  if (!(reads_column %in% colnames(metrics))) {
    stop("Error: The specified column does not exist in the 'metrics' dataset.")
  }

  # Automatically detect the column containing sample names
  possible_columns <- colnames(metrics)
  sample_column <- possible_columns[sapply(possible_columns, function(col) {
    any(metrics[[col]] %in% colnames(count_matrix))
  })]

  # Ensure only one valid column was found
  if (length(sample_column) != 1) {
    stop("Error: No unique column in 'metrics' contains the sample names.")
  }

  # Extract the total reads per sample
  reads_per_sample <- setNames(metrics[[reads_column]], metrics[[sample_column]])

  # Verify that all samples in the count matrix are present in metrics
  if (!all(colnames(count_matrix) %in% names(reads_per_sample))) {
    stop("Error: Not all samples in the count matrix have total reads in 'metrics'.")
  }

  # Normalize to RPM
  rpm_matrix <- sweep(count_matrix, 2, reads_per_sample[colnames(count_matrix)], FUN = "/") * 1e6

  # Round values to 2 decimal places
  rpm_matrix <- round(rpm_matrix, 2)

  return(rpm_matrix)
}
