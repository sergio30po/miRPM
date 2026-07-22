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

perform_statistical_tests <- function(
    miRNA_ftd,
    metadata,
    condition,
    output_file
) {
  results_dir <- "Tests_results"

  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  df_long <- as.data.frame(miRNA_ftd)

  df_long <- tibble::rownames_to_column(
    df_long,
    var = "miRNA"
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

  if (!condition %in% colnames(df_long)) {
    stop(
      "The condition column was not found in `metadata`.",
      call. = FALSE
    )
  }

  df_long$Condition <- df_long[[condition]]

  num_groups <- length(unique(df_long$Condition))

  wb <- openxlsx::createWorkbook()
  comparison_results_list <- list()

  add_significance_symbol <- function(p_value) {
    if (p_value < 0.005) {
      "***"
    } else if (p_value < 0.01) {
      "**"
    } else if (p_value < 0.05) {
      "*"
    } else {
      ""
    }
  }

  if (num_groups > 2) {
    results_kruskal <- dplyr::summarise(
      dplyr::group_by(
        df_long,
        rlang::.data$miRNA
      ),
      p_value = stats::kruskal.test(
        rlang::.data$RPM,
        rlang::.data$Condition
      )$p.value,
      .groups = "drop"
    )

    results_kruskal$FDR <- stats::p.adjust(
      results_kruskal$p_value,
      method = "fdr"
    )

    results_kruskal_sig <- results_kruskal[
      results_kruskal$FDR < 0.05,
      ,
      drop = FALSE
    ]

    results_kruskal_sig$Significance <- vapply(
      results_kruskal_sig$p_value,
      add_significance_symbol,
      character(1)
    )

    openxlsx::addWorksheet(
      wb,
      "Kruskal_Wallis"
    )

    openxlsx::writeData(
      wb,
      "Kruskal_Wallis",
      results_kruskal_sig
    )

    assign(
      "kruskal_wallis_results",
      results_kruskal_sig,
      envir = .GlobalEnv
    )

    comparison_results_list[["Kruskal_Wallis"]] <-
      results_kruskal_sig

    results_dunn <- data.frame(
      miRNA = character(),
      comparison = character(),
      Z = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )

    for (current_miRNA in results_kruskal_sig$miRNA) {
      data_miRNA <- df_long[
        df_long$miRNA == current_miRNA,
        ,
        drop = FALSE
      ]

      dunn_result <- dunn.test::dunn.test(
        data_miRNA$RPM,
        data_miRNA$Condition,
        method = "bh"
      )

      dunn_df <- data.frame(
        miRNA = rep(
          current_miRNA,
          length(dunn_result$comparisons)
        ),
        comparison = dunn_result$comparisons,
        Z = dunn_result$Z,
        p_value = dunn_result$P.adjusted,
        stringsAsFactors = FALSE
      )

      dunn_df$Significance <- vapply(
        dunn_df$p_value,
        add_significance_symbol,
        character(1)
      )

      results_dunn <- rbind(
        results_dunn,
        dunn_df
      )
    }

    results_dunn$p_value <- as.numeric(
      results_dunn$p_value
    )

    results_dunn_filtered <- results_dunn[
      results_dunn$p_value < 0.05,
      ,
      drop = FALSE
    ]

    unique_comparisons <- unique(
      results_dunn_filtered$comparison
    )

    for (comparison_name in unique_comparisons) {
      subset_comparison <- results_dunn_filtered[
        results_dunn_filtered$comparison == comparison_name,
        ,
        drop = FALSE
      ]

      if (nrow(subset_comparison) > 0) {
        openxlsx::addWorksheet(
          wb,
          comparison_name
        )

        openxlsx::writeData(
          wb,
          comparison_name,
          subset_comparison
        )

        assign(
          paste0(
            "dunn_test_",
            gsub(
              " - ",
              "_vs_",
              comparison_name
            )
          ),
          subset_comparison,
          envir = .GlobalEnv
        )

        comparison_results_list[[comparison_name]] <-
          subset_comparison
      }
    }
  } else if (num_groups == 2) {
    results_mann_whitney <- dplyr::summarise(
      dplyr::group_by(
        df_long,
        rlang::.data$miRNA
      ),
      p_value = stats::wilcox.test(
        rlang::.data$RPM,
        rlang::.data$Condition
      )$p.value,
      .groups = "drop"
    )

    results_mann_whitney$FDR <- stats::p.adjust(
      results_mann_whitney$p_value,
      method = "fdr"
    )

    results_mann_whitney_sig <- results_mann_whitney[
      results_mann_whitney$FDR < 0.05,
      ,
      drop = FALSE
    ]

    results_mann_whitney_sig$Significance <- vapply(
      results_mann_whitney_sig$p_value,
      add_significance_symbol,
      character(1)
    )

    openxlsx::addWorksheet(
      wb,
      "Mann_Whitney"
    )

    openxlsx::writeData(
      wb,
      "Mann_Whitney",
      results_mann_whitney_sig
    )

    assign(
      "mann_whitney_results",
      results_mann_whitney_sig,
      envir = .GlobalEnv
    )

    comparison_results_list[["Mann_Whitney"]] <-
      results_mann_whitney_sig
  } else {
    stop(
      "The condition must contain at least two groups.",
      call. = FALSE
    )
  }

  output_path <- file.path(
    results_dir,
    output_file
  )

  openxlsx::saveWorkbook(
    wb,
    output_path,
    overwrite = TRUE
  )

  list(
    comparison_results = comparison_results_list,
    output_folder = results_dir
  )
}
