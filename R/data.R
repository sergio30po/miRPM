#' Synthetic miRNA-Seq Example Data
#'
#' A small synthetic dataset designed to demonstrate the normalization,
#' filtering, statistical-analysis, and plotting functions provided by
#' `miRPM`.
#'
#' The dataset is intended exclusively for package examples and testing.
#' It does not represent measurements from biological samples and must not
#' be used for scientific benchmarking or biological interpretation.
#'
#' @format A list with two components:
#' \describe{
#'   \item{counts}{
#'     An integer matrix containing six miRNAs in rows and eight samples
#'     in columns.
#'   }
#'   \item{metadata}{
#'     A data frame containing one row per sample and the variables
#'     `Sample` and `Condition`.
#'   }
#' }
#'
#' @source Synthetic data generated for the `miRPM` package.
"miRPM_example"
