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

miRNA_expression_plot <- function(miRNAs_DE, miRNA_ftd, metadata, condition_column, groups, colors, output_name, plot_title = "miRNA Expression Plot") {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(plotly)
  library(readxl)
  library(htmlwidgets)

  # Create "Interactive_plots" folder if it doesn't exist
  output_folder <- "Interactive_plots"
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }

  # Load differentially expressed miRNAs from the dataframe
  miRNAs_DE <- miRNAs_DE %>% pull(1)  # Assuming miRNAs are in the first column

  # Filter count matrix for selected miRNAs
  df_long <- miRNA_ftd %>%
    as.data.frame() %>%
    rownames_to_column(var = "miRNA") %>%
    filter(miRNA %in% miRNAs_DE) %>%
    pivot_longer(cols = -miRNA, names_to = "Sample", values_to = "RPM")

  # Merge with metadata
  df_long <- df_long %>% left_join(metadata, by = "Sample")

  # Ensure only selected groups are included
  df_long <- df_long %>% filter(.data[[condition_column]] %in% groups)

  # Generate ggplot
  p <- ggplot(df_long, aes(x = miRNA, y = RPM, color = .data[[condition_column]],
                           text = paste("Sample:", Sample,
                                        "<br>miRNA:", miRNA,
                                        "<br>RPM:", round(RPM, 2),
                                        "<br>Condition:", .data[[condition_column]],
                                        "<br>Pathology:", Pathology,
                                        "<br>Gender:", Gender,
                                        "<br>Onset age:", Onset_age,
                                        "<br>Age at death:", Death_age,
                                        "<br>Braak:", Braak_stage,
                                        "<br>APOE:", APOE,
                                        "<br>Disease duration:", Disease_duration,
                                        "<br>HTT short allele:", HTT_short_allele,
                                        "<br>HTT long allele:", HTT_long_allele))) +
    geom_point(alpha = 0.7, position = position_dodge(width = 0.3)) +
    labs(
      title = plot_title,  # Add title here
      x = "miRNAs",
      y = "RPM per sample",
      color = condition_column
    ) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  # Customize title appearance
    )

  # Convert to interactive plotly plot
  interactive_plot <- ggplotly(p, tooltip = "text")

  # Save interactive plot in the "Interactive_plots" folder
  output_file <- file.path(output_folder, paste0(output_name, ".html"))
  saveWidget(interactive_plot, file = output_file)

  return(list(ggplot = p, interactive = interactive_plot, output_folder = output_folder))
}
