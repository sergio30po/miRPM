#' Plot expression for selected miRNAs
#'
#' Creates a joint dot plot for selected miRNAs across one or more groups.
#' The function returns both a static `ggplot2` object and an interactive
#' `plotly` object. Saving the interactive plot as HTML is optional.
#'
#' @param miRNAs_DE Deprecated-compatible input for the selected miRNAs.
#'   It may be a character vector or a data frame or matrix whose first
#'   column contains miRNA names. Use `mirnas` in new code.
#' @param miRNA_ftd Deprecated-compatible input for the normalized RPM
#'   matrix. Use `rpm_matrix` in new code.
#' @param metadata A data frame containing sample metadata.
#' @param condition_column Name of the metadata column defining the groups.
#' @param groups Optional character vector defining the groups to display
#'   and their order. When `NULL`, all observed groups are included.
#' @param colors Optional named character vector assigning a color to every
#'   selected group. When `NULL`, the default `ggplot2` palette is used.
#' @param output_name Optional name for the HTML file. The `.html` extension
#'   is added automatically.
#' @param plot_title Character string used as the plot title.
#' @param sample_column Name of the metadata column containing sample
#'   identifiers. Default is `"Sample"`.
#' @param tooltip_columns Optional character vector naming additional
#'   metadata columns to include in the interactive tooltip.
#' @param output_dir Directory used when `save_html = TRUE`.
#' @param save_html Logical. If `TRUE`, save the interactive plot as HTML.
#'   By default, saving is enabled when `output_name` is supplied.
#' @param selfcontained Logical passed to `htmlwidgets::saveWidget()`.
#' @param mirnas Preferred input for selected miRNAs. It may be a character
#'   vector or a data frame or matrix whose first column contains names.
#' @param rpm_matrix Preferred numeric RPM matrix, with miRNAs in rows and
#'   samples in columns.
#'
#' @details
#' Metadata are aligned to the expression matrix using sample identifiers,
#' not row position. Cohort-specific tooltip variables are not assumed;
#' additional fields can be selected through `tooltip_columns`.
#'
#' `miRNAs_DE` and `miRNA_ftd` are retained so calls written for miRPM 0.1.0
#' continue to work. New analyses should use `mirnas` and `rpm_matrix`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `ggplot`: the static `ggplot2` object.
#'   \item `interactive`: the interactive `plotly` object.
#'   \item `data`: the aligned long-format data used for plotting.
#'   \item `selected_mirnas`: miRNAs included in the plot.
#'   \item `missing_mirnas`: requested miRNAs absent from the matrix.
#'   \item `groups`: groups included in the plot.
#'   \item `output_file`: saved HTML path, or `NULL`.
#'   \item `output_folder`: saved output directory, or `NULL`.
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
#'   Condition = c("Disease", "Control", "Disease", "Control"),
#'   Sex = c("F", "M", "M", "F")
#' )
#'
#' output <- miRNA_expression_plot(
#'   mirnas = c("miR-1", "miR-2"),
#'   rpm_matrix = rpm_matrix,
#'   metadata = metadata,
#'   condition_column = "Condition",
#'   groups = c("Control", "Disease"),
#'   tooltip_columns = "Sex",
#'   save_html = FALSE
#' )
#'
#' output$ggplot
miRNA_expression_plot <- function(
    miRNAs_DE = NULL,
    miRNA_ftd = NULL,
    metadata = NULL,
    condition_column = NULL,
    groups = NULL,
    colors = NULL,
    output_name = NULL,
    plot_title = "miRNA Expression Plot",
    sample_column = "Sample",
    tooltip_columns = NULL,
    output_dir = "Interactive_plots",
    save_html = !is.null(output_name),
    selfcontained = TRUE,
    mirnas = NULL,
    rpm_matrix = NULL
) {
  if (
    !is.null(miRNAs_DE) &&
      !is.null(mirnas)
  ) {
    stop(
      "Supply either `mirnas` or `miRNAs_DE`, not both.",
      call. = FALSE
    )
  }

  if (
    !is.null(miRNA_ftd) &&
      !is.null(rpm_matrix)
  ) {
    stop(
      "Supply either `rpm_matrix` or `miRNA_ftd`, not both.",
      call. = FALSE
    )
  }

  selected_mirnas <- if (!is.null(mirnas)) {
    mirnas
  } else {
    miRNAs_DE
  }

  expression_matrix <- if (!is.null(rpm_matrix)) {
    rpm_matrix
  } else {
    miRNA_ftd
  }

  if (is.null(selected_mirnas)) {
    stop(
      "Supply `mirnas`.",
      call. = FALSE
    )
  }

  if (is.null(expression_matrix)) {
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
    !is.character(plot_title) ||
      length(plot_title) != 1L ||
      is.na(plot_title)
  ) {
    stop(
      "`plot_title` must be a single character string.",
      call. = FALSE
    )
  }

  if (
    !is.logical(save_html) ||
      length(save_html) != 1L ||
      is.na(save_html)
  ) {
    stop(
      "`save_html` must be `TRUE` or `FALSE`.",
      call. = FALSE
    )
  }

  if (
    !is.logical(selfcontained) ||
      length(selfcontained) != 1L ||
      is.na(selfcontained)
  ) {
    stop(
      "`selfcontained` must be `TRUE` or `FALSE`.",
      call. = FALSE
    )
  }

  prepared <- .prepare_mirna_plot_data(
    rpm_matrix = expression_matrix,
    metadata = metadata,
    mirnas = selected_mirnas,
    sample_column = sample_column,
    condition_column = condition_column,
    groups = groups
  )

  selected_colors <- .validate_condition_colors(
    colors = colors,
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

  plot_data <- prepared$data

  plot_data$.miRPM_tooltip <- .build_mirna_tooltip(
    data = plot_data,
    sample_column = sample_column,
    condition_column = condition_column,
    tooltip_columns = tooltip_columns
  )

  plot_data$.miRPM_group <- factor(
    as.character(plot_data[[condition_column]]),
    levels = prepared$groups
  )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = miRNA,
      y = RPM,
      color = .miRPM_group,
      text = .miRPM_tooltip
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.7,
      position = ggplot2::position_dodge(
        width = 0.4
      )
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "miRNA",
      y = "RPM per sample",
      color = condition_column
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      )
    )

  if (!is.null(selected_colors)) {
    plot_object <- plot_object +
      ggplot2::scale_color_manual(
        values = selected_colors,
        breaks = prepared$groups,
        drop = FALSE
      )
  }

  interactive_plot <- plotly::ggplotly(
    plot_object,
    tooltip = "text"
  )

  output_file <- NULL
  output_folder <- NULL

  if (isTRUE(save_html)) {
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

    file_name <- .sanitize_plot_filename(
      output_name,
      extension = "html"
    )

    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    output_file <- file.path(
      output_dir,
      file_name
    )

    htmlwidgets::saveWidget(
      interactive_plot,
      file = output_file,
      selfcontained = selfcontained
    )

    output_file <- normalizePath(
      output_file,
      winslash = "/",
      mustWork = TRUE
    )

    output_folder <- dirname(
      output_file
    )
  }

  list(
    ggplot = plot_object,
    interactive = interactive_plot,
    data = plot_data,
    selected_mirnas = prepared$selected_mirnas,
    missing_mirnas = prepared$missing_mirnas,
    groups = prepared$groups,
    output_file = output_file,
    output_folder = output_folder
  )
}
