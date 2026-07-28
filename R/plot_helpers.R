.validate_plot_expression_matrix <- function(
    x,
    argument_name = "rpm_matrix"
) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x)) {
    stop(
      paste0(
        "`",
        argument_name,
        "` must be a numeric matrix or data frame."
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(x)) {
    stop(
      paste0(
        "`",
        argument_name,
        "` must contain only numeric values."
      ),
      call. = FALSE
    )
  }

  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop(
      paste0(
        "`",
        argument_name,
        "` must contain at least one row and one column."
      ),
      call. = FALSE
    )
  }

  feature_names <- rownames(x)
  sample_names <- colnames(x)

  if (
    is.null(feature_names) ||
      anyNA(feature_names) ||
      any(feature_names == "")
  ) {
    stop(
      paste0(
        "`",
        argument_name,
        "` must have non-empty miRNA row names."
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(feature_names)) {
    stop(
      paste0(
        "`",
        argument_name,
        "` row names must be unique."
      ),
      call. = FALSE
    )
  }

  if (
    is.null(sample_names) ||
      anyNA(sample_names) ||
      any(sample_names == "")
  ) {
    stop(
      paste0(
        "`",
        argument_name,
        "` must have non-empty sample column names."
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(sample_names)) {
    stop(
      paste0(
        "`",
        argument_name,
        "` column names must be unique."
      ),
      call. = FALSE
    )
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop(
      paste0(
        "`",
        argument_name,
        "` cannot contain missing or infinite values."
      ),
      call. = FALSE
    )
  }

  if (any(x < 0)) {
    stop(
      paste0(
        "`",
        argument_name,
        "` cannot contain negative values."
      ),
      call. = FALSE
    )
  }

  x
}

.validate_plot_column_name <- function(
    x,
    argument_name
) {
  if (
    !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      x == ""
  ) {
    stop(
      paste0(
        "`",
        argument_name,
        "` must be a single non-empty column name."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.validate_plot_metadata <- function(
    metadata,
    sample_names,
    sample_column,
    condition_column
) {
  if (!is.data.frame(metadata)) {
    stop(
      "`metadata` must be a data frame.",
      call. = FALSE
    )
  }

  .validate_plot_column_name(
    sample_column,
    "sample_column"
  )

  .validate_plot_column_name(
    condition_column,
    "condition_column"
  )

  reserved_columns <- c(
    "miRNA",
    "RPM",
    ".miRPM_tooltip"
  )

  if (identical(sample_column, condition_column)) {
    stop(
      "`sample_column` and `condition_column` must be different.",
      call. = FALSE
    )
  }

  if (sample_column %in% reserved_columns) {
    stop(
      paste0(
        "`sample_column` cannot be one of: ",
        paste(reserved_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (condition_column %in% reserved_columns) {
    stop(
      paste0(
        "`condition_column` cannot be one of: ",
        paste(reserved_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  required_columns <- c(
    sample_column,
    condition_column
  )

  missing_columns <- setdiff(
    required_columns,
    names(metadata)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "Missing metadata columns: ",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  metadata_samples <- as.character(
    metadata[[sample_column]]
  )

  if (
    anyNA(metadata_samples) ||
      any(metadata_samples == "")
  ) {
    stop(
      paste0(
        "`metadata$",
        sample_column,
        "` cannot contain missing or empty sample identifiers."
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(metadata_samples)) {
    stop(
      paste0(
        "`metadata$",
        sample_column,
        "` contains duplicated sample identifiers."
      ),
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    sample_names,
    metadata_samples
  )

  if (length(missing_samples) > 0L) {
    stop(
      paste0(
        "Samples missing from `metadata`: ",
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

  group_values <- as.character(
    metadata[[condition_column]]
  )

  if (
    anyNA(group_values) ||
      any(group_values == "")
  ) {
    stop(
      paste0(
        "`metadata$",
        condition_column,
        "` cannot contain missing or empty group values."
      ),
      call. = FALSE
    )
  }

  metadata[[sample_column]] <- as.character(
    metadata[[sample_column]]
  )

  metadata
}

.resolve_plot_mirnas <- function(
    mirnas,
    available_mirnas
) {
  if (is.data.frame(mirnas)) {
    if (ncol(mirnas) == 0L) {
      stop(
        "`mirnas` data frame must contain at least one column.",
        call. = FALSE
      )
    }

    mirnas <- mirnas[[1]]
  } else if (is.matrix(mirnas)) {
    if (ncol(mirnas) == 0L) {
      stop(
        "`mirnas` matrix must contain at least one column.",
        call. = FALSE
      )
    }

    mirnas <- mirnas[, 1]
  }

  if (is.factor(mirnas)) {
    mirnas <- as.character(mirnas)
  }

  if (!is.character(mirnas)) {
    stop(
      paste0(
        "`mirnas` must be a character vector or a table ",
        "whose first column contains miRNA names."
      ),
      call. = FALSE
    )
  }

  if (
    length(mirnas) == 0L ||
      anyNA(mirnas) ||
      any(mirnas == "")
  ) {
    stop(
      "`mirnas` must contain at least one non-missing miRNA name.",
      call. = FALSE
    )
  }

  mirnas <- unique(mirnas)

  missing_mirnas <- setdiff(
    mirnas,
    available_mirnas
  )

  selected_mirnas <- mirnas[
    mirnas %in% available_mirnas
  ]

  if (length(selected_mirnas) == 0L) {
    stop(
      "None of the requested miRNAs were found in `rpm_matrix`.",
      call. = FALSE
    )
  }

  list(
    selected = selected_mirnas,
    missing = missing_mirnas
  )
}

.resolve_plot_groups <- function(
    groups,
    observed_groups
) {
  observed_groups <- unique(
    as.character(observed_groups)
  )

  if (is.null(groups)) {
    return(observed_groups)
  }

  if (is.factor(groups)) {
    groups <- as.character(groups)
  }

  if (
    !is.character(groups) ||
      length(groups) == 0L ||
      anyNA(groups) ||
      any(groups == "")
  ) {
    stop(
      "`groups` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  groups <- unique(groups)

  missing_groups <- setdiff(
    groups,
    observed_groups
  )

  if (length(missing_groups) > 0L) {
    stop(
      paste0(
        "Groups missing from `metadata`: ",
        paste(missing_groups, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  groups
}

.validate_condition_colors <- function(
    colors,
    groups
) {
  if (is.null(colors)) {
    return(NULL)
  }

  if (
    !is.character(colors) ||
      length(colors) == 0L ||
      is.null(names(colors)) ||
      anyNA(names(colors)) ||
      any(names(colors) == "")
  ) {
    stop(
      "`colors` must be a named character vector.",
      call. = FALSE
    )
  }

  if (anyDuplicated(names(colors))) {
    stop(
      "`colors` names must be unique.",
      call. = FALSE
    )
  }

  missing_colors <- setdiff(
    groups,
    names(colors)
  )

  if (length(missing_colors) > 0L) {
    stop(
      paste0(
        "`colors` is missing values for: ",
        paste(missing_colors, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  colors[groups]
}

.prepare_mirna_plot_data <- function(
    rpm_matrix,
    metadata,
    mirnas,
    sample_column,
    condition_column,
    groups = NULL
) {
  rpm_matrix <- .validate_plot_expression_matrix(
    rpm_matrix,
    "rpm_matrix"
  )

  metadata <- .validate_plot_metadata(
    metadata = metadata,
    sample_names = colnames(rpm_matrix),
    sample_column = sample_column,
    condition_column = condition_column
  )

  selected_groups <- .resolve_plot_groups(
    groups = groups,
    observed_groups = metadata[[condition_column]]
  )

  sample_is_selected <- as.character(
    metadata[[condition_column]]
  ) %in% selected_groups

  selected_metadata <- metadata[
    sample_is_selected,
    ,
    drop = FALSE
  ]

  selected_samples <- selected_metadata[[sample_column]]

  selected_matrix <- rpm_matrix[
    ,
    selected_samples,
    drop = FALSE
  ]

  resolved_mirnas <- .resolve_plot_mirnas(
    mirnas = mirnas,
    available_mirnas = rownames(selected_matrix)
  )

  selected_matrix <- selected_matrix[
    resolved_mirnas$selected,
    ,
    drop = FALSE
  ]

  long_data <- data.frame(
    miRNA = rep(
      rownames(selected_matrix),
      times = ncol(selected_matrix)
    ),
    sample_identifier = rep(
      colnames(selected_matrix),
      each = nrow(selected_matrix)
    ),
    RPM = as.vector(selected_matrix),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  names(long_data)[2] <- sample_column

  metadata_index <- match(
    long_data[[sample_column]],
    selected_metadata[[sample_column]]
  )

  metadata_columns <- setdiff(
    names(selected_metadata),
    sample_column
  )

  long_data <- cbind(
    long_data,
    selected_metadata[
      metadata_index,
      metadata_columns,
      drop = FALSE
    ]
  )

  long_data$miRNA <- factor(
    long_data$miRNA,
    levels = resolved_mirnas$selected
  )

  long_data[[condition_column]] <- factor(
    as.character(long_data[[condition_column]]),
    levels = selected_groups
  )

  list(
    data = long_data,
    rpm_matrix = selected_matrix,
    metadata = selected_metadata,
    selected_mirnas = resolved_mirnas$selected,
    missing_mirnas = resolved_mirnas$missing,
    groups = selected_groups
  )
}

.escape_plot_tooltip_value <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"

  x <- gsub(
    "&",
    "&amp;",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "<",
    "&lt;",
    x,
    fixed = TRUE
  )

  x <- gsub(
    ">",
    "&gt;",
    x,
    fixed = TRUE
  )

  x
}

.format_plot_tooltip_value <- function(x) {
  if (is.numeric(x)) {
    x <- format(
      round(x, digits = 4),
      trim = TRUE,
      scientific = FALSE
    )
  }

  .escape_plot_tooltip_value(x)
}

.build_mirna_tooltip <- function(
    data,
    sample_column,
    condition_column,
    tooltip_columns = NULL
) {
  if (!is.null(tooltip_columns)) {
    if (
      !is.character(tooltip_columns) ||
        anyNA(tooltip_columns) ||
        any(tooltip_columns == "")
    ) {
      stop(
        "`tooltip_columns` must be `NULL` or a character vector.",
        call. = FALSE
      )
    }

    missing_columns <- setdiff(
      tooltip_columns,
      names(data)
    )

    if (length(missing_columns) > 0L) {
      stop(
        paste0(
          "Tooltip columns missing from `metadata`: ",
          paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  tooltip_columns <- unique(
    c(
      sample_column,
      "miRNA",
      "RPM",
      condition_column,
      tooltip_columns
    )
  )

  tooltip_parts <- lapply(
    tooltip_columns,
    function(column_name) {
      values <- .format_plot_tooltip_value(
        data[[column_name]]
      )

      label <- switch(
        column_name,
        "miRNA" = "miRNA",
        "RPM" = "RPM",
        gsub(
          "_",
          " ",
          column_name,
          fixed = TRUE
        )
      )

      paste0(
        "<b>",
        .escape_plot_tooltip_value(label),
        ":</b> ",
        values
      )
    }
  )

  Reduce(
    function(x, y) {
      paste(
        x,
        y,
        sep = "<br>"
      )
    },
    tooltip_parts
  )
}

.sanitize_plot_filename <- function(
    x,
    extension
) {
  if (
    !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      trimws(x) == ""
  ) {
    stop(
      "`output_name` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  x <- trimws(x)

  extension_pattern <- paste0(
    "\\.",
    extension,
    "$"
  )

  x <- sub(
    extension_pattern,
    "",
    x,
    ignore.case = TRUE
  )

  x <- gsub(
    "[<>:\"/\\\\|?*]",
    "_",
    x
  )

  x <- gsub(
    "[. ]+$",
    "",
    x
  )

  if (x == "") {
    stop(
      "`output_name` does not contain a valid file name.",
      call. = FALSE
    )
  }

  paste0(
    x,
    ".",
    extension
  )
}
