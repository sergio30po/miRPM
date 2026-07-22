#' Generate dot plots for differentially expressed miRNAs
#'
#' This function generates dot plots for differentially expressed miRNAs between specific groups.
#' The plots are saved as interactive HTML files in a folder called "Interactive_plots" with a specified output name.
#' A customizable title can be added to the plot.
#'
#' @param miRNAs_DE File containing differentially expressed miRNAs (miRNAs should be in the first column).
#' @param miRNA_ftd Count matrix with miRNAs as row names and samples as column names.
#' @param metadata Data frame containing sample information. Must have a column with sample names and a condition column.
#' @param condition_column Name of the column in `metadata` that contains the group classifications.
#' @param groups Vector specifying the two groups to compare within the `condition_column`.
#' @param colors Named vector specifying the colors for each group.
#' @param output_name Name for the output HTML file.
#' @param plot_title Title for the plot (default: "miRNA Expression Plot").
#'
#' @return A list containing:
#' \itemize{
#'   \item `ggplot`: A ggplot object of the dot plot.
#'   \item `interactive`: An interactive plotly object of the dot plot.
#'   \item `output_folder`: The path to the "Interactive_plots" folder where the HTML file is saved.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' output <- miRNA_expression_plot(
#'   miRNAs_DE = "differential_miRNAs.xlsx",
#'   miRNA_ftd = miRNA_ftd,
#'   metadata = metadata,
#'   condition_column = "Condition",
#'   groups = c("Control", "Disease"),
#'   colors = c("Control" = "blue", "Disease" = "red"),
#'   output_name = "miRNA_expression",
#'   plot_title = "Differentially Expressed miRNAs in Control vs Disease"
#' )
#'
#' output$ggplot
#' output$interactive
#' output$output_folder
#' }

miRNA_expression_plot <- function(
    miRNAs_DE,
    miRNA_ftd,
    metadata,
    condition_column,
    groups,
    colors,
    output_name,
    plot_title = "miRNA Expression Plot"
) {
  output_folder <- "Interactive_plots"

  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  miRNAs_DE <- miRNAs_DE[[1]]

  df_long <- as.data.frame(miRNA_ftd)
  df_long <- tibble::rownames_to_column(df_long, var = "miRNA")

  df_long <- dplyr::filter(
    df_long,
    rlang::.data$miRNA %in% miRNAs_DE
  )

  df_long <- tidyr::pivot_longer(
    df_long,
    cols = -1,
    names_to = "Sample",
    values_to = "RPM"
  )

  df_long <- dplyr::left_join(
    df_long,
    metadata,
    by = "Sample"
  )

  df_long <- dplyr::filter(
    df_long,
    df_long[[condition_column]] %in% groups
  )

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = rlang::.data$miRNA,
      y = rlang::.data$RPM,
      color = rlang::.data[[condition_column]],
      text = paste(
        "Sample:", rlang::.data$Sample,
        "<br>miRNA:", rlang::.data$miRNA,
        "<br>RPM:", round(rlang::.data$RPM, 2),
        "<br>Condition:", rlang::.data[[condition_column]],
        "<br>Pathology:", rlang::.data$Pathology,
        "<br>Gender:", rlang::.data$Gender,
        "<br>Onset age:", rlang::.data$Onset_age,
        "<br>Age at death:", rlang::.data$Death_age,
        "<br>Braak:", rlang::.data$Braak_stage,
        "<br>APOE:", rlang::.data$APOE,
        "<br>Disease duration:", rlang::.data$Disease_duration,
        "<br>HTT short allele:", rlang::.data$HTT_short_allele,
        "<br>HTT long allele:", rlang::.data$HTT_long_allele
      )
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.7,
      position = ggplot2::position_dodge(width = 0.3)
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "miRNAs",
      y = "RPM per sample",
      color = condition_column
    ) +
    ggplot2::scale_color_manual(values = colors) +
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

  interactive_plot <- plotly::ggplotly(
    p,
    tooltip = "text"
  )

  output_file <- file.path(
    output_folder,
    paste0(output_name, ".html")
  )

  htmlwidgets::saveWidget(
    interactive_plot,
    file = output_file
  )

  list(
    ggplot = p,
    interactive = interactive_plot,
    output_folder = output_folder
  )
}
