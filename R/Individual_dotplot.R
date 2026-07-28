#' Create individual expression plots for selected miRNAs
#'
#' Creates one dot plot per selected miRNA across one or more groups. Plots
#' can be returned without writing files or optionally saved as PNG images.
#'
#' @param filtered_results Optional data frame, matrix or character vector
#'   defining the miRNAs to plot. For a table, the miRNA names are taken from
#'   `mirna_column`; when that column is absent, the first column is used for
#'   backward compatibility with miRPM 0.1.0.
#' @param rpm_matrix Numeric matrix or data frame containing normalized RPM
#'   values, with miRNAs in rows and samples in columns.
#' @param metadata Data frame containing sample metadata.
#' @param condition_column Name of the metadata column defining the groups.
#' @param sample_column Name of the metadata column containing sample
#'   identifiers. Default is `"Sample"`.
#' @param groups_to_include Optional character vector defining the groups to
#'   display and their order. When `NULL`, all observed groups are included.
#' @param condition_colors Optional named character vector assigning a color
#'   to every selected group. When `NULL`, the default `ggplot2` palette is
#'   used.
#' @param adjusted_pvalue_column Optional name of the adjusted p-value column
#'   in `filtered_results`. When `NULL`, no adjusted p-value is shown.
#'   Default is `"FDR"`.
#' @param output_dir Main directory used when `save_png = TRUE`.
#' @param save_png Logical. If `TRUE`, save one PNG file per miRNA.
#' @param mirna_column Name of the miRNA column in `filtered_results`.
#'   Default is `"miRNA"`.
#' @param mirnas Preferred direct character vector of miRNA names. Supply
#'   either `mirnas` or `filtered_results`, not both.
#' @param point_size Numeric point size.
#' @param point_alpha Numeric point transparency between 0 and 1.
#' @param jitter_width Numeric horizontal jitter width. Use zero to disable
#'   jitter.
#' @param width Numeric PNG width in inches.
#' @param height Numeric PNG height in inches.
#' @param dpi Numeric PNG resolution.
#'
#' @details
#' Metadata are aligned to the RPM matrix using sample identifiers, not row
#' position. Missing requested miRNAs are reported and omitted.
#'
#' The historical miRPM 0.1.0 positional interface remains supported:
#'
#' `individual_miRNA_plot(filtered_results, rpm_matrix, metadata,
#' condition_column, sample_column, groups_to_include, condition_colors,
#' adjusted_pvalue_column, output_dir)`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `plots`: named list of `ggplot2` objects.
#'   \item `files`: named character vector of saved PNG paths, or an empty
#'   character vector when `save_png = FALSE`.
#'   \item `data`: aligned long-format data used for plotting.
#'   \item `selected_mirnas`: miRNAs included in the plots.
#'   \item `missing_mirnas`: requested miRNAs absent from the matrix.
#'   \item `groups`: groups included in the plots.
#'   \item `adjusted_pvalues`: named adjusted p-values used in captions.
#'   \item `output_folder`: output directory, or `NULL`.
#' }
#'
#' @export
#'
#' @examples
#' rpm_matrix <- matrix(
#'   c(
#'     5, 6, 10, 11,
#'     2, 3, 8, 9
#'   ),
#'   nrow = 2,
#'   byrow = TRUE,
#'   dimnames = list(
#'     c("miR-1", "miR-2"),
#'     c("S1", "S2", "S3", "S4")
#'   )
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("S3", "S1", "S4", "S2"),
#'   Condition = c("Disease", "Control", "Disease", "Control")
#' )
#'
#' results <- data.frame(
#'   miRNA = c("miR-1", "miR-2"),
#'   FDR = c(0.01, 0.04)
#' )
#'
#' output <- individual_miRNA_plot(
#'   filtered_results = results,
#'   rpm_matrix = rpm_matrix,
#'   metadata = metadata,
#'   condition_column = "Condition",
#'   groups_to_include = c("Control", "Disease"),
#'   save_png = FALSE
#' )
#'
#' output$plots[["miR-1"]]
individual_miRNA_plot <- function(
    filtered_results = NULL,
    rpm_matrix = NULL,
    metadata = NULL,
    condition_column = NULL,
    sample_column = "Sample",
    groups_to_include = NULL,
    condition_colors = NULL,
    adjusted_pvalue_column = "FDR",
    output_dir = "Individual_plots",
    save_png = TRUE,
    mirna_column = "miRNA",
    mirnas = NULL,
    point_size = 3,
    point_alpha = 0.7,
    jitter_width = 0.12,
    width = 6,
    height = 6,
    dpi = 300
) {
  validate_single_number <- function(
      x,
      argument_name,
      minimum = -Inf,
      maximum = Inf,
      strictly_positive = FALSE
  ) {
    invalid <- !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < minimum ||
      x > maximum

    if (strictly_positive) {
      invalid <- invalid || x <= 0
    }

    if (invalid) {
      stop(
        paste0(
          "`",
          argument_name,
          "` must be a single finite ",
          if (strictly_positive) "positive " else "",
          "number."
        ),
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  resolve_filtered_results <- function(
      filtered_results,
      mirna_column,
      adjusted_pvalue_column
  ) {
    if (is.factor(filtered_results)) {
      filtered_results <- as.character(
        filtered_results
      )
    }

    if (is.character(filtered_results)) {
      if (
        length(filtered_results) == 0L ||
          anyNA(filtered_results) ||
          any(filtered_results == "")
      ) {
        stop(
          paste0(
            "`filtered_results` must contain at least one ",
            "non-missing miRNA name."
          ),
          call. = FALSE
        )
      }

      return(
        list(
          mirnas = unique(filtered_results),
          adjusted_pvalues = numeric()
        )
      )
    }

    if (is.matrix(filtered_results)) {
      filtered_results <- as.data.frame(
        filtered_results,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    if (!is.data.frame(filtered_results)) {
      stop(
        paste0(
          "`filtered_results` must be a data frame, matrix ",
          "or character vector."
        ),
        call. = FALSE
      )
    }

    if (ncol(filtered_results) == 0L) {
      stop(
        "`filtered_results` must contain at least one column.",
        call. = FALSE
      )
    }

    .validate_plot_column_name(
      mirna_column,
      "mirna_column"
    )

    selected_mirna_column <- mirna_column

    if (!selected_mirna_column %in% names(filtered_results)) {
      selected_mirna_column <- names(filtered_results)[1]
    }

    selected_mirnas <- filtered_results[[
      selected_mirna_column
    ]]

    if (is.factor(selected_mirnas)) {
      selected_mirnas <- as.character(
        selected_mirnas
      )
    }

    if (
      !is.character(selected_mirnas) ||
        length(selected_mirnas) == 0L ||
        anyNA(selected_mirnas) ||
        any(selected_mirnas == "")
    ) {
      stop(
        paste0(
          "The selected miRNA column in `filtered_results` ",
          "must contain non-missing character values."
        ),
        call. = FALSE
      )
    }

    if (anyDuplicated(selected_mirnas)) {
      stop(
        paste0(
          "The selected miRNA column in `filtered_results` ",
          "contains duplicated miRNA names."
        ),
        call. = FALSE
      )
    }

    adjusted_pvalues <- numeric()

    if (!is.null(adjusted_pvalue_column)) {
      .validate_plot_column_name(
        adjusted_pvalue_column,
        "adjusted_pvalue_column"
      )

      if (
        adjusted_pvalue_column %in%
          names(filtered_results)
      ) {
        adjusted_values <- filtered_results[[
          adjusted_pvalue_column
        ]]

        if (
          !is.numeric(adjusted_values) ||
            anyNA(adjusted_values) ||
            any(!is.finite(adjusted_values)) ||
            any(adjusted_values < 0) ||
            any(adjusted_values > 1)
        ) {
          stop(
            paste0(
              "`filtered_results$",
              adjusted_pvalue_column,
              "` must contain finite values between 0 and 1."
            ),
            call. = FALSE
          )
        }

        adjusted_pvalues <- adjusted_values
        names(adjusted_pvalues) <- selected_mirnas
      } else if (
        !identical(adjusted_pvalue_column, "FDR")
      ) {
        stop(
          paste0(
            "`adjusted_pvalue_column` was not found in ",
            "`filtered_results`: ",
            adjusted_pvalue_column
          ),
          call. = FALSE
        )
      }
    }

    list(
      mirnas = selected_mirnas,
      adjusted_pvalues = adjusted_pvalues
    )
  }

  if (
    !is.null(filtered_results) &&
      !is.null(mirnas)
  ) {
    stop(
      "Supply either `mirnas` or `filtered_results`, not both.",
      call. = FALSE
    )
  }

  if (
    is.null(filtered_results) &&
      is.null(mirnas)
  ) {
    stop(
      "Supply `mirnas` or `filtered_results`.",
      call. = FALSE
    )
  }

  if (is.null(rpm_matrix)) {
    stop(
      "Supply `rpm_matrix`.",
      call. = FALSE
    )
  }

  if (is.null(metadata)) {
    stop(
      "Supply `metadata`.",
      call. = FALSE
    )
  }

  .validate_plot_column_name(
    condition_column,
    "condition_column"
  )

  .validate_plot_column_name(
    sample_column,
    "sample_column"
  )

  if (
    !is.logical(save_png) ||
      length(save_png) != 1L ||
      is.na(save_png)
  ) {
    stop(
      "`save_png` must be `TRUE` or `FALSE`.",
      call. = FALSE
    )
  }

  validate_single_number(
    point_size,
    "point_size",
    strictly_positive = TRUE
  )

  validate_single_number(
    point_alpha,
    "point_alpha",
    minimum = 0,
    maximum = 1
  )

  validate_single_number(
    jitter_width,
    "jitter_width",
    minimum = 0
  )

  validate_single_number(
    width,
    "width",
    strictly_positive = TRUE
  )

  validate_single_number(
    height,
    "height",
    strictly_positive = TRUE
  )

  validate_single_number(
    dpi,
    "dpi",
    strictly_positive = TRUE
  )

  resolved_results <- if (!is.null(mirnas)) {
    list(
      mirnas = mirnas,
      adjusted_pvalues = numeric()
    )
  } else {
    resolve_filtered_results(
      filtered_results = filtered_results,
      mirna_column = mirna_column,
      adjusted_pvalue_column = adjusted_pvalue_column
    )
  }

  prepared <- .prepare_mirna_plot_data(
    rpm_matrix = rpm_matrix,
    metadata = metadata,
    mirnas = resolved_results$mirnas,
    sample_column = sample_column,
    condition_column = condition_column,
    groups = groups_to_include
  )

  selected_colors <- .validate_condition_colors(
    colors = condition_colors,
    groups = prepared$groups
  )

  if (length(prepared$missing_mirnas) > 0L) {
    warning(
      paste0(
        "The following miRNAs were not found and were omitted: ",
        paste(prepared$missing_mirnas, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  adjusted_pvalues <- resolved_results$adjusted_pvalues

  if (length(adjusted_pvalues) > 0L) {
    adjusted_pvalues <- adjusted_pvalues[
      prepared$selected_mirnas
    ]
  }

  output_folder <- NULL

  if (isTRUE(save_png)) {
    if (
      !is.character(output_dir) ||
        length(output_dir) != 1L ||
        is.na(output_dir) ||
        trimws(output_dir) == ""
    ) {
      stop(
        "`output_dir` must be a single non-empty path.",
        call. = FALSE
      )
    }

    comparison_name <- paste(
      prepared$groups,
      collapse = "_vs_"
    )

    comparison_name <- sub(
      "\\.png$",
      "",
      .sanitize_plot_filename(
        comparison_name,
        extension = "png"
      ),
      ignore.case = TRUE
    )

    output_folder <- file.path(
      output_dir,
      comparison_name
    )

    dir.create(
      output_folder,
      recursive = TRUE,
      showWarnings = FALSE
    )

    output_folder <- normalizePath(
      output_folder,
      winslash = "/",
      mustWork = TRUE
    )
  }

  plots <- stats::setNames(
    vector(
      mode = "list",
      length = length(prepared$selected_mirnas)
    ),
    prepared$selected_mirnas
  )

  files <- stats::setNames(
    character(),
    character()
  )

  for (current_mirna in prepared$selected_mirnas) {
    plot_data <- prepared$data[
      as.character(prepared$data$miRNA) ==
        current_mirna,
      ,
      drop = FALSE
    ]

    plot_data$.miRPM_group <- factor(
      as.character(
        plot_data[[condition_column]]
      ),
      levels = prepared$groups
    )

    adjusted_value <- if (
      current_mirna %in% names(adjusted_pvalues)
    ) {
      unname(
        adjusted_pvalues[[current_mirna]]
      )
    } else {
      NA_real_
    }

    caption <- if (is.finite(adjusted_value)) {
      paste0(
        adjusted_pvalue_column,
        ": ",
        format(
          adjusted_value,
          digits = 4,
          scientific = adjusted_value < 0.0001
        )
      )
    } else {
      NULL
    }

    position_object <- if (jitter_width > 0) {
      ggplot2::position_jitter(
        width = jitter_width,
        height = 0,
        seed = 1
      )
    } else {
      "identity"
    }

    plot_object <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .miRPM_group,
        y = RPM,
        color = .miRPM_group
      )
    ) +
      ggplot2::geom_point(
        size = point_size,
        alpha = point_alpha,
        position = position_object
      ) +
      ggplot2::labs(
        title = current_mirna,
        x = "Group",
        y = "RPM",
        color = condition_column,
        caption = caption
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(
          fill = "white",
          color = NA
        ),
        panel.background = ggplot2::element_rect(
          fill = "white",
          color = NA
        ),
        panel.grid.major = ggplot2::element_line(
          color = "gray90"
        ),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(
          color = "black"
        ),
        axis.text = ggplot2::element_text(
          color = "black",
          size = 12
        ),
        axis.title = ggplot2::element_text(
          color = "black",
          size = 14,
          face = "bold"
        ),
        plot.title = ggplot2::element_text(
          color = "black",
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        plot.caption = ggplot2::element_text(
          color = "gray50",
          size = 10,
          hjust = 1
        ),
        legend.position = "none"
      )

    if (!is.null(selected_colors)) {
      plot_object <- plot_object +
        ggplot2::scale_color_manual(
          values = selected_colors,
          breaks = prepared$groups,
          drop = FALSE
        )
    }

    plots[[current_mirna]] <- plot_object

    if (isTRUE(save_png)) {
      file_name <- .sanitize_plot_filename(
        current_mirna,
        extension = "png"
      )

      output_file <- file.path(
        output_folder,
        file_name
      )

      ggplot2::ggsave(
        filename = output_file,
        plot = plot_object,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white"
      )

      files[[current_mirna]] <- normalizePath(
        output_file,
        winslash = "/",
        mustWork = TRUE
      )
    }
  }

  list(
    plots = plots,
    files = files,
    data = prepared$data,
    selected_mirnas = prepared$selected_mirnas,
    missing_mirnas = prepared$missing_mirnas,
    groups = prepared$groups,
    adjusted_pvalues = adjusted_pvalues,
    output_folder = output_folder
  )
}
