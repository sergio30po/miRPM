#' Individual miRNA plot for differentially expressed miRNAs
#'
#' This function generates dot plots for differentially expressed miRNAs between specific groups.
#' The plots are saved as PNG files in a folder named after the comparison, inside a main folder called "Individual_plots".
#'
#' @param filtered_results Data frame with the filtered results (FDR < 0.05). The first column should contain miRNA names.
#' @param rpm_matrix Matrix of RPMs (miRNAs in rows, samples in columns).
#' @param metadata Data frame with sample group information. Must have a column for sample names and another for groups.
#' @param condition_column Name of the column in `metadata` that contains the groups.
#' @param sample_column Name of the column in `metadata` that contains the sample names.
#' @param groups_to_include Vector with the names of the groups to include in the comparison.
#' @param condition_colors Named vector of colors for each condition. The names should match the groups in `groups_to_include`.
#' @param adjusted_pvalue_column Name of the column in `filtered_results` that contains the adjusted p-value.
#' @param output_dir Name of the main folder where the plots will be saved (default: "Individual_plots").
#'
#' @return No return value. The plots are saved as PNG files in the specified folder structure.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' individual_miRNA_plot(
#'   filtered_results = filtered_results,
#'   rpm_matrix = rpm_matrix,
#'   metadata = metadata,
#'   condition_column = "Condition",
#'   sample_column = "Samples",
#'   groups_to_include = c("Group1", "Group2"),
#'   condition_colors = c("Group1" = "blue", "Group2" = "red"),
#'   adjusted_pvalue_column = "adj.P.Val"
#' )
#' }

individual_miRNA_plot <- function(
    filtered_results,
    rpm_matrix,
    metadata,
    condition_column,
    sample_column,
    groups_to_include,
    condition_colors,
    adjusted_pvalue_column,
    output_dir = "Individual_plots"
) {
  if (!all(groups_to_include %in% names(condition_colors))) {
    stop(
      "`condition_colors` must contain a color for every selected group.",
      call. = FALSE
    )
  }

  if (!adjusted_pvalue_column %in% colnames(filtered_results)) {
    stop(
      "`adjusted_pvalue_column` was not found in `filtered_results`.",
      call. = FALSE
    )
  }

  required_metadata_columns <- c(sample_column, condition_column)
  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(metadata)
  )

  if (length(missing_metadata_columns) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_metadata_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  miRNAs_DE <- filtered_results[[1]]

  df_long <- as.data.frame(rpm_matrix)
  df_long <- tibble::rownames_to_column(
    df_long,
    var = "miRNA"
  )

  df_long <- dplyr::filter(
    df_long,
    rlang::.data$miRNA %in% miRNAs_DE
  )

  df_long <- tidyr::pivot_longer(
    df_long,
    cols = -1,
    names_to = sample_column,
    values_to = "RPM"
  )

  df_long <- dplyr::left_join(
    df_long,
    metadata,
    by = sample_column
  )

  df_long <- dplyr::filter(
    df_long,
    rlang::.data[[condition_column]] %in% groups_to_include
  )

  comparison_folder <- file.path(
    output_dir,
    paste(groups_to_include, collapse = "_vs_")
  )

  dir.create(
    comparison_folder,
    showWarnings = FALSE,
    recursive = TRUE
  )

  for (current_miRNA in miRNAs_DE) {
    if (!current_miRNA %in% df_long$miRNA) {
      warning(
        "miRNA '",
        current_miRNA,
        "' was not found and will be skipped.",
        call. = FALSE
      )
      next
    }

    adjusted_pvalue <- filtered_results[
      filtered_results[[1]] == current_miRNA,
      adjusted_pvalue_column,
      drop = TRUE
    ]

    df_plot <- dplyr::filter(
      df_long,
      rlang::.data$miRNA == current_miRNA
    )

    p <- ggplot2::ggplot(
      df_plot,
      ggplot2::aes(
        x = rlang::.data[[condition_column]],
        y = rlang::.data$RPM,
        color = rlang::.data[[condition_column]]
      )
    ) +
      ggplot2::geom_point(
        size = 3,
        alpha = 0.7
      ) +
      ggplot2::scale_color_manual(
        values = condition_colors
      ) +
      ggplot2::labs(
        title = current_miRNA,
        x = "Group",
        y = "RPM",
        caption = paste0(
          "Adjusted p-value: ",
          round(adjusted_pvalue[1], 4)
        )
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

    ggplot2::ggsave(
      filename = file.path(
        comparison_folder,
        paste0(current_miRNA, ".png")
      ),
      plot = p,
      width = 6,
      height = 6,
      dpi = 300
    )
  }

  invisible(NULL)
}
