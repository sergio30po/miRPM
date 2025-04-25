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

individual_miRNA_plot <- function(filtered_results, rpm_matrix, metadata, condition_column, sample_column, groups_to_include, condition_colors, adjusted_pvalue_column, output_dir = "Individual_plots") {
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("tidyverse", quietly = TRUE)) {
    stop("The 'ggplot2' and 'tidyverse' packages are required. Please install them with install.packages('ggplot2') and install.packages('tidyverse').")
  }

  # Load necessary libraries
  library(ggplot2)
  library(tidyverse)

  # Validate condition_colors
  if (!all(groups_to_include %in% names(condition_colors))) {
    stop("The `condition_colors` vector must include colors for all groups in `groups_to_include`.")
  }

  # Validate adjusted_pvalue_column
  if (!adjusted_pvalue_column %in% colnames(filtered_results)) {
    stop("The `adjusted_pvalue_column` does not exist in the `filtered_results` data frame.")
  }

  # 1. Create the main folder if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 2. Extract the names of the miRNAs from the first column of the filtered results
  miRNAs_DE <- filtered_results %>% pull(1)  # Assuming miRNAs are in the first column

  # 3. Filter count matrix for selected miRNAs
  df_long <- rpm_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "miRNA") %>%
    filter(miRNA %in% miRNAs_DE) %>%
    pivot_longer(cols = -miRNA, names_to = sample_column, values_to = "RPM")

  # 4. Merge with metadata
  df_long <- df_long %>% left_join(metadata, by = sample_column)

  # 5. Ensure only selected groups are included
  df_long <- df_long %>% filter(.data[[condition_column]] %in% groups_to_include)

  # 6. Create the comparison folder inside the main folder
  comparison_folder <- file.path(output_dir, paste(groups_to_include, collapse = "_vs_"))
  dir.create(comparison_folder, showWarnings = FALSE, recursive = TRUE)

  # 7. Generate a dot plot for each miRNA in the list
  for (miRNA in miRNAs_DE) {
    # Skip if miRNA is not found in df_long
    if (!miRNA %in% df_long$miRNA) {
      warning(paste("miRNA", miRNA, "not found in the filtered data. Skipping this miRNA."))
      next
    }

    # Extract the adjusted p-value for the current miRNA
    adjusted_pvalue <- filtered_results[filtered_results[[1]] == miRNA, adjusted_pvalue_column]

    # Filter data for the current miRNA
    df_plot <- df_long %>% filter(miRNA == !!miRNA)

    # Create the plot
    p <- ggplot(df_plot, aes(x = .data[[condition_column]], y = RPM, color = .data[[condition_column]])) +
      geom_point(size = 3, alpha = 0.7) +  # Points with custom colors
      scale_color_manual(values = condition_colors) +  # Use custom colors
      labs(
        title = miRNA,
        x = "Group",
        y = "RPM",
        caption = paste0("Adjusted p-value: ", round(adjusted_pvalue, 4))  # Add adjusted p-value
      ) +
      theme_minimal() +  # Minimalist theme
      theme(
        plot.background = element_rect(fill = "white", color = NA),  # White background
        panel.background = element_rect(fill = "white", color = NA),  # White background
        panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Axis lines in black
        axis.text = element_text(color = "black", size = 12),  # Axis text in black and larger
        axis.title = element_text(color = "black", size = 14, face = "bold"),  # Axis titles in bold
        plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),  # Centered and bold title
        plot.caption = element_text(color = "gray50", size = 10, hjust = 1),  # Caption in gray and right-aligned
        legend.position = "none"  # Remove legend
      )

    # Save the plot as a PNG file in the comparison folder
    ggsave(
      filename = file.path(comparison_folder, paste0(miRNA, ".png")),
      plot = p,
      width = 6,
      height = 6,
      dpi = 300
    )
  }
}
