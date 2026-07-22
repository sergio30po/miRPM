#' Normalize a count matrix to reads per million
#'
#' Converts a complete count matrix to reads per million (RPM) using the
#' library size of each sample.
#'
#' @param count_matrix A numeric matrix or data frame with features in rows
#'   and samples in columns.
#' @param metrics Deprecated. A data frame containing a `Sample` column and
#'   the library-size column indicated by `reads_column`.
#' @param reads_column Deprecated. Name of the library-size column in
#'   `metrics`.
#' @param library_sizes Optional named numeric vector containing one positive
#'   library size per sample. Names must match the column names of
#'   `count_matrix`.
#'
#' @details
#' When `library_sizes` is `NULL`, library sizes are calculated from the
#' complete input matrix using `colSums(count_matrix)`.
#'
#' Normalization should therefore be performed before filtering miRNAs.
#' Applying the function to an already filtered matrix changes the denominator
#' and the biological meaning of the resulting RPM values.
#'
#' The deprecated `metrics` and `reads_column` arguments are retained
#' temporarily for compatibility with miRPM 0.1.0 pipelines.
#'
#' The returned matrix stores normalization information in the
#' `"miRPM_normalization"` attribute. This includes the library sizes,
#' their source and the approximate number of reads represented by one RPM
#' in each sample.
#'
#' @return A numeric matrix with the same dimensions and dimnames as
#'   `count_matrix`. The `"miRPM_normalization"` attribute records the
#'   normalization method, library sizes, their source and reads per RPM.
#'
#' @examples
#' count_matrix <- matrix(
#'   c(100, 200, 300, 400),
#'   nrow = 2,
#'   dimnames = list(
#'     c("miR-1", "miR-2"),
#'     c("Sample1", "Sample2")
#'   )
#' )
#'
#' rpm_matrix <- normalize_rpm(count_matrix)
#'
#' colSums(rpm_matrix)
#'
#' attr(rpm_matrix, "miRPM_normalization")
#'
#' @export
normalize_rpm <- function(
    count_matrix,
    metrics = NULL,
    reads_column = NULL,
    library_sizes = NULL
) {
  if (is.data.frame(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }

  if (!is.matrix(count_matrix)) {
    stop(
      "`count_matrix` must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }

  if (!is.numeric(count_matrix)) {
    stop(
      "`count_matrix` must contain only numeric values.",
      call. = FALSE
    )
  }

  if (nrow(count_matrix) == 0L || ncol(count_matrix) == 0L) {
    stop(
      "`count_matrix` must contain at least one row and one column.",
      call. = FALSE
    )
  }

  feature_names <- rownames(count_matrix)
  sample_names <- colnames(count_matrix)

  if (
    is.null(feature_names) ||
    anyNA(feature_names) ||
    any(feature_names == "")
  ) {
    stop(
      "`count_matrix` must have non-empty row names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(feature_names)) {
    stop(
      "`count_matrix` row names must be unique.",
      call. = FALSE
    )
  }

  if (
    is.null(sample_names) ||
    anyNA(sample_names) ||
    any(sample_names == "")
  ) {
    stop(
      "`count_matrix` must have non-empty column names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(sample_names)) {
    stop(
      "`count_matrix` column names must be unique.",
      call. = FALSE
    )
  }

  if (anyNA(count_matrix) || any(!is.finite(count_matrix))) {
    stop(
      "`count_matrix` cannot contain missing or infinite values.",
      call. = FALSE
    )
  }

  if (any(count_matrix < 0)) {
    stop(
      "`count_matrix` cannot contain negative values.",
      call. = FALSE
    )
  }

  using_legacy_interface <-
    !is.null(metrics) || !is.null(reads_column)

  external_library_sizes_supplied <- !is.null(library_sizes)

  if (using_legacy_interface) {
    if (!is.null(library_sizes)) {
      stop(
        paste0(
          "Use either `library_sizes` or the deprecated metrics interface, ",
          "not both."
        ),
        call. = FALSE
      )
    }

    if (is.null(metrics) || is.null(reads_column)) {
      stop(
        "`metrics` and `reads_column` must be supplied together.",
        call. = FALSE
      )
    }

    warning(
      paste0(
        "`metrics` and `reads_column` are deprecated. ",
        "Use `library_sizes` or omit both arguments."
      ),
      call. = FALSE
    )

    if (!is.data.frame(metrics)) {
      stop(
        "`metrics` must be a data frame.",
        call. = FALSE
      )
    }

    if (!"Sample" %in% names(metrics)) {
      stop(
        "`metrics` must contain a `Sample` column.",
        call. = FALSE
      )
    }

    if (
      !is.character(reads_column) ||
      length(reads_column) != 1L ||
      is.na(reads_column) ||
      reads_column == ""
    ) {
      stop(
        "`reads_column` must be a single non-empty column name.",
        call. = FALSE
      )
    }

    if (!reads_column %in% names(metrics)) {
      stop(
        "The column specified by `reads_column` does not exist.",
        call. = FALSE
      )
    }

    if (!is.numeric(metrics[[reads_column]])) {
      stop(
        "The library-size column in `metrics` must be numeric.",
        call. = FALSE
      )
    }

    metric_samples <- as.character(metrics$Sample)

    if (
      anyNA(metric_samples) ||
      any(metric_samples == "") ||
      anyDuplicated(metric_samples)
    ) {
      stop(
        "`metrics$Sample` must contain unique, non-missing sample names.",
        call. = FALSE
      )
    }

    missing_samples <- setdiff(sample_names, metric_samples)

    if (length(missing_samples) > 0L) {
      stop(
        paste0(
          "Not all count-matrix samples are present in `metrics`: ",
          paste(missing_samples, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    sample_order <- match(sample_names, metric_samples)

    library_sizes <- metrics[[reads_column]][sample_order]
    names(library_sizes) <- sample_names
  } else if (is.null(library_sizes)) {
    library_sizes <- colSums(count_matrix)
    names(library_sizes) <- sample_names
  } else {
    if (!is.numeric(library_sizes) || !is.vector(library_sizes)) {
      stop(
        "`library_sizes` must be a named numeric vector.",
        call. = FALSE
      )
    }

    if (
      is.null(names(library_sizes)) ||
      anyNA(names(library_sizes)) ||
      any(names(library_sizes) == "")
    ) {
      stop(
        "`library_sizes` must have non-empty sample names.",
        call. = FALSE
      )
    }

    if (anyDuplicated(names(library_sizes))) {
      stop(
        "`library_sizes` names must be unique.",
        call. = FALSE
      )
    }

    missing_samples <- setdiff(sample_names, names(library_sizes))

    if (length(missing_samples) > 0L) {
      stop(
        paste0(
          "`library_sizes` is missing samples: ",
          paste(missing_samples, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    library_sizes <- library_sizes[sample_names]
  }

  if (anyNA(library_sizes) || any(!is.finite(library_sizes))) {
    stop(
      "`library_sizes` cannot contain missing or infinite values.",
      call. = FALSE
    )
  }

  invalid_samples <- names(library_sizes)[library_sizes <= 0]

  if (length(invalid_samples) > 0L) {
    stop(
      paste0(
        "Library sizes must be greater than zero. Invalid samples: ",
        paste(invalid_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  library_size_source <- if (using_legacy_interface) {
    "legacy_metrics"
  } else if (external_library_sizes_supplied) {
    "external_vector"
  } else {
    "matrix_column_sums"
  }

  rpm_matrix <- sweep(
    count_matrix,
    MARGIN = 2,
    STATS = library_sizes,
    FUN = "/"
  ) * 1e6

  dimnames(rpm_matrix) <- dimnames(count_matrix)

  attr(rpm_matrix, "miRPM_normalization") <- list(
    method = "RPM",
    scale = 1e6,
    library_sizes = library_sizes,
    reads_per_rpm = library_sizes / 1e6,
    library_size_source = library_size_source,
    input_features = nrow(count_matrix),
    input_samples = ncol(count_matrix)
  )

  rpm_matrix
}
