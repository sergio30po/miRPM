#' Perform Statistical Tests on miRNA Expression Data
#'
#' This function performs statistical tests on miRNA expression data based on the number of groups in the specified condition.
#' For two groups, it performs the Mann-Whitney U test. For three or more groups, it performs the Kruskal-Wallis test followed by Dunn's post hoc test.
#' Results are saved in an Excel file inside the "Tests_results" folder with separate sheets for each comparison. Each sheet will also include a column indicating the significance of the results using symbols: `*` for p-value < 0.05, `**` for p-value < 0.01, and `***` for p-value < 0.005.
#'
#' @param miRNA_ftd A count matrix of miRNA expression data (miRNAs in rows, samples in columns).
#' @param metadata A dataframe containing sample information, including the condition to test.
#' @param condition The name of the column in `metadata` that defines the groups to compare.
#' @param output_file The name of the output Excel file (e.g., "Results.xlsx").
#' @return An Excel file with statistical results. For Kruskal-Wallis, it includes a sheet for the global test and sheets for Dunn's post hoc comparisons. For Mann-Whitney, it includes a single sheet with the test results. Each sheet will have a column titled "Significance" showing the corresponding symbol (`*`, `**`, or `***`).
#' @return A dataframe for each comparison.
#' @return The path to the output folder ("Tests_results").
#' @export
#' @examples
#' # Example usage:
#' # miRNA_ftd <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, dimnames = list(c("miRNA1", "miRNA2"), c("Sample1", "Sample2", "Sample3")))
#' # metadata <- data.frame(Sample = c("Sample1", "Sample2", "Sample3"), Condition = c("A", "A", "B"))
#' # perform_statistical_tests(miRNA_ftd, metadata, "Condition", "Test_Results.xlsx")

perform_statistical_tests <- function(miRNA_ftd, metadata, condition, output_file) {
  # Load necessary libraries
  library(dunn.test)
  library(tidyverse)
  library(openxlsx)

  # Ensure "Tests_results" folder exists
  results_dir <- "Tests_results"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }

  # Convert count matrix to long format
  df_long <- miRNA_ftd %>%
    as.data.frame() %>%
    rownames_to_column(var = "miRNA") %>%
    pivot_longer(cols = -miRNA, names_to = "Sample", values_to = "RPM") %>%
    left_join(metadata, by = "Sample")

  # Extract the condition column
  df_long$Condition <- df_long[[condition]]

  # Number of groups
  num_groups <- length(unique(df_long$Condition))

  # Initialize workbook
  wb <- createWorkbook()

  # List to store dataframes for comparisons
  comparison_results_list <- list()

  # Function to add significance symbols
  add_significance_symbol <- function(p_value) {
    if (p_value < 0.005) return("***")
    else if (p_value < 0.01) return("**")
    else if (p_value < 0.05) return("*")
    else return("")
  }

  if (num_groups > 2) {
    # Kruskal-Wallis test
    results_kruskal <- df_long %>%
      group_by(miRNA) %>%
      summarise(p_value = kruskal.test(RPM ~ Condition)$p.value, .groups = "drop")

    results_kruskal$FDR <- p.adjust(results_kruskal$p_value, method = "fdr")
    results_kruskal_sig <- results_kruskal %>% filter(FDR < 0.05)
    results_kruskal_sig$Significance <- sapply(results_kruskal_sig$p_value, add_significance_symbol)

    addWorksheet(wb, "Kruskal_Wallis")
    writeData(wb, "Kruskal_Wallis", results_kruskal_sig)

    # Assign dataframe to environment
    assign("kruskal_wallis_results", results_kruskal_sig, envir = .GlobalEnv)
    comparison_results_list[["Kruskal_Wallis"]] <- results_kruskal_sig

    # Dunn's test
    results_dunn <- data.frame(miRNA = character(), comparison = character(), Z = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

    for (miRNA in results_kruskal_sig$miRNA) {
      data_miRNA <- df_long %>% filter(miRNA == !!miRNA)  # Subset for each miRNA

      dunn_test <- dunn.test(data_miRNA$RPM, data_miRNA$Condition, method = "bh")

      dunn_df <- data.frame(
        miRNA = rep(miRNA, length(dunn_test$comparisons)),  # Ensure correct length
        comparison = dunn_test$comparisons,
        Z = dunn_test$Z,
        p_value = dunn_test$P.adjusted
      )

      dunn_df$Significance <- sapply(dunn_df$p_value, add_significance_symbol)
      results_dunn <- rbind(results_dunn, dunn_df)
    }

    # Ensure p_value is numeric
    results_dunn$p_value <- as.numeric(results_dunn$p_value)

    # Filter for p_value < 0.05
    results_dunn_filtered <- results_dunn[results_dunn$p_value < 0.05, ]

    # Save Dunn's test results to separate sheets
    unique_comparisons <- unique(results_dunn_filtered$comparison)
    for (comp in unique_comparisons) {
      subset_comp <- results_dunn_filtered[results_dunn_filtered$comparison == comp, ]

      if (nrow(subset_comp) > 0) {
        addWorksheet(wb, comp)
        writeData(wb, comp, subset_comp)
        assign(paste0("dunn_test_", gsub(" - ", "_vs_", comp)), subset_comp, envir = .GlobalEnv)
        comparison_results_list[[comp]] <- subset_comp
      }
    }

  } else if (num_groups == 2) {
    # Mann-Whitney U test
    results_mann_whitney <- df_long %>%
      group_by(miRNA) %>%
      summarise(p_value = wilcox.test(RPM ~ Condition)$p.value, .groups = "drop")

    results_mann_whitney$FDR <- p.adjust(results_mann_whitney$p_value, method = "fdr")
    results_mann_whitney_sig <- results_mann_whitney %>% filter(FDR < 0.05)
    results_mann_whitney_sig$Significance <- sapply(results_mann_whitney_sig$p_value, add_significance_symbol)

    addWorksheet(wb, "Mann_Whitney")
    writeData(wb, "Mann_Whitney", results_mann_whitney_sig)

    assign("mann_whitney_results", results_mann_whitney_sig, envir = .GlobalEnv)
    comparison_results_list[["Mann_Whitney"]] <- results_mann_whitney_sig
  } else {
    stop("Error: The condition must have at least 2 groups.")
  }

  # Save workbook in the "Tests_results" folder
  output_path <- file.path(results_dir, output_file)
  saveWorkbook(wb, output_path, overwrite = TRUE)

  # Return list of dataframes and output folder path
  return(list(
    comparison_results = comparison_results_list,
    output_folder = results_dir
  ))
}
