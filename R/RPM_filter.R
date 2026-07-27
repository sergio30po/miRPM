#' Filter miRNAs by abundance and prevalence
#'
#' Filters a normalized RPM matrix using configurable abundance, prevalence,
#' sample-count and group-specific criteria.
#'
#' @param rpm_matrix A numeric matrix or data frame containing RPM values,
#'   with miRNAs in rows and samples in columns.
#' @param metadata Optional data frame containing sample information. Grouped
#'   filtering requires a `Sample` column and the column named by
#'   `group_column`.
#' @param threshold Deprecated. Minimum prevalence used by the miRPM 0.1.0
#'   interface.
#' @param min_reads Deprecated. Minimum expression value used by the miRPM
#'   0.1.0 interface. In the historical workflow this was normally applied to
#'   an RPM-normalized matrix despite the argument name.
#' @param threshold_comparison Deprecated comparison operator for
#'   `threshold`.
#' @param read_comparison Deprecated comparison operator for `min_reads`.
#' @param count_matrix Deprecated named alias for the normalized matrix used
#'   by miRPM 0.1.0.
#' @param raw_counts Optional raw count matrix corresponding to `rpm_matrix`.
#'   When supplied with `min_count`, a sample must satisfy both the RPM and
#'   raw-count thresholds.
#' @param min_rpm Minimum RPM required for a miRNA to be considered detected
#'   in a sample. The default is 5.
#' @param min_count Optional minimum raw count required in addition to
#'   `min_rpm`. When `raw_counts` is omitted, counts are inferred from the
#'   normalization metadata created by `normalize_rpm()`.
#' @param min_prevalence Minimum proportion of samples in which a miRNA must
#'   satisfy the abundance criterion. Must be between 0 and 1.
#' @param min_samples Optional minimum absolute number of samples that must
#'   satisfy the abundance criterion. When both `min_prevalence` and
#'   `min_samples` are supplied, both requirements are enforced.
#' @param group_column Optional name of the metadata column defining groups.
#' @param prevalence_scope Filtering scope. One of `"overall"`,
#'   `"any_group"` or `"all_groups"`.
#' @param min_groups Minimum number of groups in which the criterion must be
#'   satisfied when `prevalence_scope = "any_group"`.
#' @param return_diagnostics If `TRUE`, return a list containing the filtered
#'   matrix, diagnostics, summary and parameters. If `FALSE`, return only the
#'   filtered matrix and store the same information as attributes.
#'
#' @details
#' The recommended workflow is to normalize the complete count matrix with
#' `normalize_rpm()` before filtering. Filtering thresholds are configurable
#' because their interpretation depends on sequencing depth, sample type,
#' library complexity and study objectives.
#'
#' A provisional starting profile for differential-expression analysis is
#' `min_rpm = 5`, `min_count = 15`, `min_prevalence = 0.5`,
#' `min_samples = 3` and `prevalence_scope = "any_group"`. These values are
#' recommendations rather than universal requirements.
#'
#' With `prevalence_scope = "overall"`, the criterion is evaluated across all
#' samples. With `"any_group"`, it must be satisfied in at least `min_groups`
#' groups. With `"all_groups"`, it must be satisfied in every group.
#'
#' The historical miRPM 0.1.0 interface is retained temporarily for backward
#' compatibility.
#'
#' @return A filtered RPM matrix, or a list containing the filtered matrix,
#'   diagnostics, summary and parameters when `return_diagnostics = TRUE`.
#'
#' @examples
#' rpm_matrix <- matrix(
#'   c(
#'     5, 6, 0, 0,
#'     0, 0, 10, 0,
#'     0, 0, 0, 0
#'   ),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(
#'     c("miR-1", "miR-2", "miR-3"),
#'     c("S1", "S2", "S3", "S4")
#'   )
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("S1", "S2", "S3", "S4"),
#'   Group = c("Control", "Control", "Disease", "Disease")
#' )
#'
#' filtered <- filter_mirnas(
#'   rpm_matrix = rpm_matrix,
#'   metadata = metadata,
#'   min_rpm = 5,
#'   min_prevalence = 0.5,
#'   group_column = "Group",
#'   prevalence_scope = "any_group"
#' )
#'
#' @export
filter_mirnas <- function(
    rpm_matrix = NULL,
    metadata = NULL,
    threshold = NULL,
    min_reads = NULL,
    threshold_comparison = NULL,
    read_comparison = NULL,
    count_matrix = NULL,
    raw_counts = NULL,
    min_rpm = 5,
    min_count = NULL,
    min_prevalence = 0.5,
    min_samples = NULL,
    group_column = NULL,
    prevalence_scope = c("overall", "any_group", "all_groups"),
    min_groups = 1,
    return_diagnostics = FALSE
) {
  validate_expression_matrix <- function(
      x,
      argument_name,
      require_integer_like = FALSE
  ) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }

    if (!is.matrix(x)) {
      stop(
        paste0("`", argument_name, "` must be a numeric matrix or data frame."),
        call. = FALSE
      )
    }

    if (!is.numeric(x)) {
      stop(
        paste0("`", argument_name, "` must contain only numeric values."),
        call. = FALSE
      )
    }

    if (nrow(x) == 0L || ncol(x) == 0L) {
      stop(
        paste0(
          "`", argument_name,
          "` must contain at least one row and one column."
        ),
        call. = FALSE
      )
    }

    feature_names <- rownames(x)
    sample_names <- colnames(x)

    if (is.null(feature_names) || anyNA(feature_names) || any(feature_names == "")) {
      stop(
        paste0("`", argument_name, "` must have non-empty row names."),
        call. = FALSE
      )
    }

    if (anyDuplicated(feature_names)) {
      stop(
        paste0("`", argument_name, "` row names must be unique."),
        call. = FALSE
      )
    }

    if (is.null(sample_names) || anyNA(sample_names) || any(sample_names == "")) {
      stop(
        paste0("`", argument_name, "` must have non-empty column names."),
        call. = FALSE
      )
    }

    if (anyDuplicated(sample_names)) {
      stop(
        paste0("`", argument_name, "` column names must be unique."),
        call. = FALSE
      )
    }

    if (anyNA(x) || any(!is.finite(x))) {
      stop(
        paste0(
          "`", argument_name,
          "` cannot contain missing or infinite values."
        ),
        call. = FALSE
      )
    }

    if (any(x < 0)) {
      stop(
        paste0("`", argument_name, "` cannot contain negative values."),
        call. = FALSE
      )
    }

    if (
      require_integer_like &&
        any(abs(x - round(x)) > sqrt(.Machine$double.eps))
    ) {
      stop(
        paste0(
          "`", argument_name,
          "` must contain integer-like raw counts."
        ),
        call. = FALSE
      )
    }

    x
  }

  compare_values <- function(x, y, operator) {
    switch(
      operator,
      ">" = x > y,
      ">=" = x >= y,
      "<" = x < y,
      "<=" = x <= y,
      stop("Invalid comparison operator.", call. = FALSE)
    )
  }

  validate_single_number <- function(
      x,
      argument_name,
      minimum = -Inf,
      maximum = Inf,
      allow_null = FALSE
  ) {
    if (is.null(x) && allow_null) {
      return(invisible(NULL))
    }

    if (
      !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x < minimum ||
        x > maximum
    ) {
      stop(
        paste0("Invalid value supplied for `", argument_name, "`."),
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  validate_positive_integer <- function(
      x,
      argument_name,
      allow_null = FALSE
  ) {
    if (is.null(x) && allow_null) {
      return(invisible(NULL))
    }

    if (
      !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x < 1 ||
        x != floor(x)
    ) {
      stop(
        paste0("`", argument_name, "` must be a positive whole number."),
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  align_metadata <- function(metadata, sample_names, group_column) {
    if (!is.data.frame(metadata)) {
      stop("`metadata` must be a data frame.", call. = FALSE)
    }

    if (!"Sample" %in% names(metadata)) {
      stop("`metadata` must contain a `Sample` column.", call. = FALSE)
    }

    if (!group_column %in% names(metadata)) {
      stop(
        paste0(
          "The group column `", group_column,
          "` does not exist in `metadata`."
        ),
        call. = FALSE
      )
    }

    metadata_samples <- as.character(metadata$Sample)

    if (
      anyNA(metadata_samples) ||
        any(metadata_samples == "") ||
        anyDuplicated(metadata_samples)
    ) {
      stop(
        "`metadata$Sample` must contain unique, non-missing sample names.",
        call. = FALSE
      )
    }

    missing_samples <- setdiff(sample_names, metadata_samples)

    if (length(missing_samples) > 0L) {
      stop(
        paste0(
          "The following matrix samples are missing from `metadata`: ",
          paste(missing_samples, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    metadata <- metadata[
      match(sample_names, metadata_samples),
      ,
      drop = FALSE
    ]

    group_values <- as.character(metadata[[group_column]])

    if (anyNA(group_values) || any(group_values == "")) {
      stop(
        paste0(
          "`metadata$", group_column,
          "` cannot contain missing or empty values."
        ),
        call. = FALSE
      )
    }

    metadata
  }

  align_raw_counts <- function(raw_counts, rpm_matrix) {
    raw_counts <- validate_expression_matrix(
      raw_counts,
      "raw_counts",
      require_integer_like = TRUE
    )

    missing_features <- setdiff(rownames(rpm_matrix), rownames(raw_counts))
    missing_samples <- setdiff(colnames(rpm_matrix), colnames(raw_counts))

    if (length(missing_features) > 0L) {
      stop(
        paste0(
          "`raw_counts` is missing miRNAs: ",
          paste(missing_features, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (length(missing_samples) > 0L) {
      stop(
        paste0(
          "`raw_counts` is missing samples: ",
          paste(missing_samples, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    raw_counts[
      rownames(rpm_matrix),
      colnames(rpm_matrix),
      drop = FALSE
    ]
  }

  run_legacy_filter <- function(
      expression_matrix,
      metadata,
      threshold,
      min_reads,
      threshold_comparison,
      read_comparison,
      return_diagnostics
  ) {
    expression_matrix <- validate_expression_matrix(
      expression_matrix,
      "rpm_matrix"
    )

    if (!is.data.frame(metadata)) {
      stop("`metadata` must be a data frame.", call. = FALSE)
    }

    if (!"Condition" %in% names(metadata)) {
      stop(
        "The deprecated interface requires `metadata$Condition`.",
        call. = FALSE
      )
    }

    if (ncol(expression_matrix) != nrow(metadata)) {
      stop(
        paste0(
          "The deprecated interface requires the number of matrix columns ",
          "to equal the number of metadata rows."
        ),
        call. = FALSE
      )
    }

    if (is.null(threshold)) {
      threshold <- 0.5
    }
    if (is.null(min_reads)) {
      min_reads <- 1000
    }
    if (is.null(threshold_comparison)) {
      threshold_comparison <- ">="
    }
    if (is.null(read_comparison)) {
      read_comparison <- ">"
    }

    valid_comparisons <- c(">", ">=", "<", "<=")

    if (!threshold_comparison %in% valid_comparisons) {
      stop(
        "`threshold_comparison` must be one of '>', '>=', '<' or '<='.",
        call. = FALSE
      )
    }

    if (!read_comparison %in% valid_comparisons) {
      stop(
        "`read_comparison` must be one of '>', '>=', '<' or '<='.",
        call. = FALSE
      )
    }

    validate_single_number(
      threshold,
      "threshold",
      minimum = 0,
      maximum = 1
    )
    validate_single_number(
      min_reads,
      "min_reads",
      minimum = 0
    )

    group_values <- as.character(metadata$Condition)

    if (anyNA(group_values) || any(group_values == "")) {
      stop(
        "`metadata$Condition` cannot contain missing or empty values.",
        call. = FALSE
      )
    }

    groups <- unique(group_values)

    group_pass_list <- lapply(groups, function(current_group) {
      group_samples <- which(group_values == current_group)
      required_samples <- ceiling(length(group_samples) * threshold)

      expressed <- compare_values(
        expression_matrix[, group_samples, drop = FALSE],
        min_reads,
        read_comparison
      )

      detected_samples <- rowSums(expressed)

      compare_values(
        detected_samples,
        required_samples,
        threshold_comparison
      )
    })

    group_pass_matrix <- do.call(cbind, group_pass_list)
    rownames(group_pass_matrix) <- rownames(expression_matrix)
    colnames(group_pass_matrix) <- groups

    keep <- rowSums(group_pass_matrix) > 0L
    filtered_matrix <- expression_matrix[keep, , drop = FALSE]

    normalization_info <- attr(
      expression_matrix,
      "miRPM_normalization",
      exact = TRUE
    )

    if (!is.null(normalization_info)) {
      attr(filtered_matrix, "miRPM_normalization") <- normalization_info
    }

    diagnostics <- data.frame(
      miRNA = rownames(expression_matrix),
      groups_passing = rowSums(group_pass_matrix),
      passes_filter = keep,
      reason = ifelse(keep, "retained", "failed_legacy_filter"),
      row.names = NULL,
      check.names = FALSE
    )

    filter_summary <- list(
      input_mirnas = nrow(expression_matrix),
      retained_mirnas = sum(keep),
      removed_mirnas = sum(!keep),
      legacy_interface = TRUE
    )

    filter_parameters <- list(
      threshold = threshold,
      min_reads = min_reads,
      threshold_comparison = threshold_comparison,
      read_comparison = read_comparison
    )

    attr(filtered_matrix, "filter_diagnostics") <- diagnostics
    attr(filtered_matrix, "filter_summary") <- filter_summary
    attr(filtered_matrix, "filter_parameters") <- filter_parameters

    if (isTRUE(return_diagnostics)) {
      return(list(
        filtered_matrix = filtered_matrix,
        diagnostics = diagnostics,
        summary = filter_summary,
        parameters = filter_parameters
      ))
    }

    filtered_matrix
  }

  new_interface_argument_supplied <-
    !missing(raw_counts) ||
    !missing(min_rpm) ||
    !missing(min_count) ||
    !missing(min_prevalence) ||
    !missing(min_samples) ||
    !missing(group_column) ||
    !missing(prevalence_scope) ||
    !missing(min_groups)

  historical_default_call <-
    !is.null(metadata) &&
    !new_interface_argument_supplied

  legacy_arguments_used <-
    !is.null(count_matrix) ||
    !is.null(threshold) ||
    !is.null(min_reads) ||
    !is.null(threshold_comparison) ||
    !is.null(read_comparison) ||
    historical_default_call

  if (!is.null(rpm_matrix) && !is.null(count_matrix)) {
    stop(
      "Supply either `rpm_matrix` or deprecated `count_matrix`, not both.",
      call. = FALSE
    )
  }

  if (is.null(rpm_matrix) && is.null(count_matrix)) {
    stop("Supply `rpm_matrix`.", call. = FALSE)
  }

  if (legacy_arguments_used) {
    warning(
      paste0(
        "The miRPM 0.1.0 filtering interface is deprecated. ",
        "Use `rpm_matrix`, `min_rpm`, `min_prevalence`, `min_samples`, ",
        "`group_column` and `prevalence_scope`."
      ),
      call. = FALSE
    )

    legacy_matrix <- if (!is.null(count_matrix)) {
      count_matrix
    } else {
      rpm_matrix
    }

    return(run_legacy_filter(
      expression_matrix = legacy_matrix,
      metadata = metadata,
      threshold = threshold,
      min_reads = min_reads,
      threshold_comparison = threshold_comparison,
      read_comparison = read_comparison,
      return_diagnostics = return_diagnostics
    ))
  }

  rpm_matrix <- validate_expression_matrix(rpm_matrix, "rpm_matrix")
  prevalence_scope <- match.arg(prevalence_scope)

  validate_single_number(min_rpm, "min_rpm", minimum = 0)
  validate_single_number(
    min_count,
    "min_count",
    minimum = 0,
    allow_null = TRUE
  )
  validate_single_number(
    min_prevalence,
    "min_prevalence",
    minimum = 0,
    maximum = 1
  )
  validate_positive_integer(
    min_samples,
    "min_samples",
    allow_null = TRUE
  )
  validate_positive_integer(min_groups, "min_groups")

  if (
    !is.logical(return_diagnostics) ||
      length(return_diagnostics) != 1L ||
      is.na(return_diagnostics)
  ) {
    stop("`return_diagnostics` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!is.null(raw_counts) && is.null(min_count)) {
    stop(
      "`raw_counts` is only used when `min_count` is supplied.",
      call. = FALSE
    )
  }

  normalization_info <- attr(
    rpm_matrix,
    "miRPM_normalization",
    exact = TRUE
  )

  count_matrix_for_filter <- NULL
  count_requirement_source <- "not_used"

  if (!is.null(min_count)) {
    if (!is.null(raw_counts)) {
      count_matrix_for_filter <- align_raw_counts(raw_counts, rpm_matrix)
      count_requirement_source <- "raw_counts"
    } else {
      if (
        is.null(normalization_info) ||
          is.null(normalization_info$library_sizes)
      ) {
        stop(
          paste0(
            "`min_count` requires `raw_counts` or normalization metadata ",
            "created by `normalize_rpm()`."
          ),
          call. = FALSE
        )
      }

      library_sizes <- normalization_info$library_sizes

      if (
        is.null(names(library_sizes)) ||
          !all(colnames(rpm_matrix) %in% names(library_sizes))
      ) {
        stop(
          paste0(
            "The normalization metadata does not contain library sizes ",
            "for all samples."
          ),
          call. = FALSE
        )
      }

      library_sizes <- library_sizes[colnames(rpm_matrix)]

      count_matrix_for_filter <- sweep(
        rpm_matrix,
        MARGIN = 2,
        STATS = library_sizes / 1e6,
        FUN = "*"
      )

      count_requirement_source <- "inferred_from_normalization"
    }
  }

  detected <- rpm_matrix >= min_rpm

  if (!is.null(min_count)) {
    count_tolerance <- sqrt(.Machine$double.eps)
    detected <- detected & (
      count_matrix_for_filter >= min_count - count_tolerance
    )
  }

  if (prevalence_scope == "overall") {
    group_indices <- list(overall = seq_len(ncol(rpm_matrix)))
  } else {
    if (
      is.null(group_column) ||
        !is.character(group_column) ||
        length(group_column) != 1L ||
        is.na(group_column) ||
        group_column == ""
    ) {
      stop(
        paste0(
          "`group_column` must be supplied when `prevalence_scope` is ",
          "'any_group' or 'all_groups'."
        ),
        call. = FALSE
      )
    }

    metadata <- align_metadata(
      metadata = metadata,
      sample_names = colnames(rpm_matrix),
      group_column = group_column
    )

    group_values <- as.character(metadata[[group_column]])
    group_indices <- split(seq_along(group_values), group_values)
  }

  number_of_groups <- length(group_indices)

  if (
    prevalence_scope == "any_group" &&
      min_groups > number_of_groups
  ) {
    stop(
      paste0(
        "`min_groups` cannot exceed the number of available groups (",
        number_of_groups,
        ")."
      ),
      call. = FALSE
    )
  }

  group_results <- lapply(names(group_indices), function(current_group) {
    sample_indices <- group_indices[[current_group]]
    group_size <- length(sample_indices)

    required_by_prevalence <- ceiling(group_size * min_prevalence)
    required_by_samples <- if (is.null(min_samples)) 0 else min_samples
    required_samples <- max(required_by_prevalence, required_by_samples)

    detected_samples <- rowSums(
      detected[, sample_indices, drop = FALSE]
    )

    list(
      group = current_group,
      group_size = group_size,
      required_samples = required_samples,
      detected_samples = detected_samples,
      prevalence = detected_samples / group_size,
      passes = detected_samples >= required_samples
    )
  })

  group_pass_matrix <- do.call(
    cbind,
    lapply(group_results, `[[`, "passes")
  )

  if (is.null(dim(group_pass_matrix))) {
    group_pass_matrix <- matrix(
      group_pass_matrix,
      ncol = 1L
    )
  }

  rownames(group_pass_matrix) <- rownames(rpm_matrix)
  colnames(group_pass_matrix) <- names(group_indices)

  groups_passing <- rowSums(group_pass_matrix)

  keep <- switch(
    prevalence_scope,
    overall = group_pass_matrix[, 1L],
    any_group = groups_passing >= min_groups,
    all_groups = groups_passing == number_of_groups
  )

  diagnostics <- data.frame(
    miRNA = rownames(rpm_matrix),
    mean_rpm = rowMeans(rpm_matrix),
    median_rpm = apply(rpm_matrix, 1L, stats::median),
    maximum_rpm = apply(rpm_matrix, 1L, max),
    overall_detected_samples = rowSums(detected),
    overall_prevalence = rowMeans(detected),
    groups_passing = groups_passing,
    passes_filter = keep,
    reason = ifelse(
      keep,
      "retained",
      "below_abundance_or_prevalence_requirement"
    ),
    row.names = NULL,
    check.names = FALSE
  )

  group_names <- names(group_indices)
  safe_group_names <- make.names(group_names, unique = TRUE)
  names(safe_group_names) <- group_names

  for (group_result in group_results) {
    safe_group_name <- safe_group_names[[group_result$group]]

    diagnostics[[paste0("group_size_", safe_group_name)]] <-
      group_result$group_size

    diagnostics[[paste0("required_samples_", safe_group_name)]] <-
      group_result$required_samples

    diagnostics[[paste0("detected_samples_", safe_group_name)]] <-
      group_result$detected_samples

    diagnostics[[paste0("prevalence_", safe_group_name)]] <-
      group_result$prevalence

    diagnostics[[paste0("passes_", safe_group_name)]] <-
      group_result$passes
  }

  filtered_matrix <- rpm_matrix[keep, , drop = FALSE]

  if (!is.null(normalization_info)) {
    attr(filtered_matrix, "miRPM_normalization") <- normalization_info
  }

  rpm_threshold_read_equivalents <- NULL

  if (
    !is.null(normalization_info) &&
      !is.null(normalization_info$reads_per_rpm)
  ) {
    reads_per_rpm <- normalization_info$reads_per_rpm

    if (
      !is.null(names(reads_per_rpm)) &&
        all(colnames(rpm_matrix) %in% names(reads_per_rpm))
    ) {
      reads_per_rpm <- reads_per_rpm[colnames(rpm_matrix)]
      rpm_threshold_read_equivalents <- min_rpm * reads_per_rpm
    }
  }

  filter_summary <- list(
    input_mirnas = nrow(rpm_matrix),
    retained_mirnas = sum(keep),
    removed_mirnas = sum(!keep),
    retained_fraction = sum(keep) / nrow(rpm_matrix),
    min_rpm = min_rpm,
    min_count = min_count,
    min_prevalence = min_prevalence,
    min_samples = min_samples,
    prevalence_scope = prevalence_scope,
    min_groups = min_groups,
    number_of_groups = number_of_groups,
    count_requirement_source = count_requirement_source,
    rpm_threshold_read_equivalents = rpm_threshold_read_equivalents,
    legacy_interface = FALSE
  )

  filter_parameters <- list(
    min_rpm = min_rpm,
    min_count = min_count,
    min_prevalence = min_prevalence,
    min_samples = min_samples,
    group_column = group_column,
    prevalence_scope = prevalence_scope,
    min_groups = min_groups
  )

  attr(filtered_matrix, "filter_diagnostics") <- diagnostics
  attr(filtered_matrix, "filter_summary") <- filter_summary
  attr(filtered_matrix, "filter_parameters") <- filter_parameters

  if (isTRUE(return_diagnostics)) {
    return(list(
      filtered_matrix = filtered_matrix,
      diagnostics = diagnostics,
      summary = filter_summary,
      parameters = filter_parameters
    ))
  }

  filtered_matrix
}
