#' Perform statistical tests on miRNA expression data
#'
#' Performs non-parametric statistical tests according to the number of
#' groups defined in the metadata.
#'
#' For two groups, the function performs a Mann-Whitney-Wilcoxon test.
#' For three or more groups, it performs a Kruskal-Wallis test followed
#' by Dunn's post hoc tests for significant miRNAs.
#'
#' Results are adjusted for multiple comparisons using the
#' Benjamini-Hochberg false discovery rate method and are exported to
#' an Excel workbook.
#'
#' @param miRNA_ftd A numeric matrix or data frame with miRNAs in rows
#'   and samples in columns.
#' @param metadata A data frame containing sample information. It must
#'   include a column named `"Sample"`.
#' @param condition Character string naming the metadata column that
#'   defines the comparison groups.
#' @param output_file Character string naming the output Excel file.
#' @param assign_results Logical. If `TRUE`, result objects are also
#'   assigned to the calling environment for backward compatibility.
#'   Default is `TRUE`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `comparison_results`: data frames with statistical results.
#'   \item `output_folder`: directory where the workbook was saved.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' counts <- matrix(
#'   c(10, 20, 30, 40, 50, 60),
#'   nrow = 2,
#'   dimnames = list(
#'     c("miRNA1", "miRNA2"),
#'     c("Sample1", "Sample2", "Sample3")
#'   )
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Sample1", "Sample2", "Sample3"),
#'   Condition = c("Control", "Control", "Disease")
#' )
#'
#' perform_statistical_tests(
#'   miRNA_ftd = counts,
#'   metadata = metadata,
#'   condition = "Condition",
#'   output_file = "Test_results.xlsx"
#' )
#' }
perform_statistical_tests <- function(
    miRNA_ftd,
    metadata,
    condition,
    output_file,
    assign_results = TRUE
) {
  caller_env <- parent.frame()

  if (!is.matrix(miRNA_ftd) && !is.data.frame(miRNA_ftd)) {
    stop(
      "`miRNA_ftd` must be a matrix or data frame.",
      call. = FALSE
    )
  }

  rpm_matrix <- as.matrix(miRNA_ftd)

  if (!is.numeric(rpm_matrix)) {
    stop(
      "`miRNA_ftd` must contain numeric expression values.",
      call. = FALSE
    )
  }

  if (is.null(rownames(rpm_matrix))) {
    stop(
      "`miRNA_ftd` must have miRNA names as row names.",
      call. = FALSE
    )
  }

  if (is.null(colnames(rpm_matrix))) {
    stop(
      "`miRNA_ftd` must have sample names as column names.",
      call. = FALSE
    )
  }

  if (anyNA(rpm_matrix)) {
    stop(
      "`miRNA_ftd` cannot contain missing values.",
      call. = FALSE
    )
  }

  if (!is.data.frame(metadata)) {
    stop(
      "`metadata` must be a data frame.",
      call. = FALSE
    )
  }

  if (!"Sample" %in% colnames(metadata)) {
    stop(
      "`metadata` must contain a column named `Sample`.",
      call. = FALSE
    )
  }

  if (!condition %in% colnames(metadata)) {
    stop(
      "The condition column was not found in `metadata`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(metadata$Sample)) {
    stop(
      "`metadata$Sample` contains duplicated sample identifiers.",
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    colnames(rpm_matrix),
    metadata$Sample
  )

  if (length(missing_samples) > 0) {
    stop(
      "Samples missing from `metadata`: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  results_dir <- "Tests_results"

  if (!dir.exists(results_dir)) {
    dir.create(
      results_dir,
      recursive = TRUE
    )
  }

  df_long <- data.frame(
    miRNA = rep(
      rownames(rpm_matrix),
      times = ncol(rpm_matrix)
    ),
    Sample = rep(
      colnames(rpm_matrix),
      each = nrow(rpm_matrix)
    ),
    RPM = as.vector(rpm_matrix),
    stringsAsFactors = FALSE
  )

  metadata_index <- match(
    df_long$Sample,
    metadata$Sample
  )

  metadata_columns <- setdiff(
    colnames(metadata),
    "Sample"
  )

  df_long <- cbind(
    df_long,
    metadata[
      metadata_index,
      metadata_columns,
      drop = FALSE
    ]
  )

  df_long$Condition <- df_long[[condition]]

  if (anyNA(df_long$Condition)) {
    stop(
      "The condition column contains missing values.",
      call. = FALSE
    )
  }

  groups <- unique(df_long$Condition)
  num_groups <- length(groups)

  if (num_groups < 2) {
    stop(
      "The condition must contain at least two groups.",
      call. = FALSE
    )
  }

  wb <- openxlsx::createWorkbook()
  comparison_results_list <- list()

  add_significance_symbol <- function(p_value) {
    if (is.na(p_value)) {
      return("")
    }

    if (p_value < 0.005) {
      return("***")
    }

    if (p_value < 0.01) {
      return("**")
    }

    if (p_value < 0.05) {
      return("*")
    }

    ""
  }

  miRNA_names <- unique(df_long$miRNA)

  if (num_groups > 2) {
    kruskal_pvalues <- vapply(
      miRNA_names,
      function(current_miRNA) {
        current_data <- df_long[
          df_long$miRNA == current_miRNA,
          ,
          drop = FALSE
        ]

        test_result <- stats::kruskal.test(
          current_data$RPM,
          current_data$Condition
        )

        test_result$p.value
      },
      numeric(1)
    )

    results_kruskal <- data.frame(
      miRNA = miRNA_names,
      p_value = kruskal_pvalues,
      stringsAsFactors = FALSE
    )

    results_kruskal$FDR <- stats::p.adjust(
      results_kruskal$p_value,
      method = "fdr"
    )

    results_kruskal$Significance <- vapply(
      results_kruskal$p_value,
      add_significance_symbol,
      character(1)
    )

    results_kruskal_sig <- results_kruskal[
      !is.na(results_kruskal$FDR) &
        results_kruskal$FDR < 0.05,
      ,
      drop = FALSE
    ]

    openxlsx::addWorksheet(
      wb,
      "Kruskal_Wallis"
    )

    openxlsx::writeData(
      wb,
      "Kruskal_Wallis",
      results_kruskal_sig
    )

    comparison_results_list[["Kruskal_Wallis"]] <-
      results_kruskal_sig

    if (isTRUE(assign_results)) {
      assign(
        "kruskal_wallis_results",
        results_kruskal_sig,
        envir = caller_env
      )
    }

    results_dunn <- data.frame(
      miRNA = character(),
      comparison = character(),
      Z = numeric(),
      p_value = numeric(),
      Significance = character(),
      stringsAsFactors = FALSE
    )

    for (current_miRNA in results_kruskal_sig$miRNA) {
      current_data <- df_long[
        df_long$miRNA == current_miRNA,
        ,
        drop = FALSE
      ]

      dunn_result <- dunn.test::dunn.test(
        current_data$RPM,
        current_data$Condition,
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

    results_dunn_filtered <- results_dunn[
      !is.na(results_dunn$p_value) &
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

      sheet_name <- gsub(
        "[:\\\\/?*\\[\\]]",
        "_",
        comparison_name
      )

      sheet_name <- substr(
        sheet_name,
        1,
        31
      )

      openxlsx::addWorksheet(
        wb,
        sheet_name
      )

      openxlsx::writeData(
        wb,
        sheet_name,
        subset_comparison
      )

      comparison_results_list[[comparison_name]] <-
        subset_comparison

      if (isTRUE(assign_results)) {
        object_name <- paste0(
          "dunn_test_",
          gsub(
            " - ",
            "_vs_",
            comparison_name,
            fixed = TRUE
          )
        )

        assign(
          object_name,
          subset_comparison,
          envir = caller_env
        )
      }
    }
  } else {
    mann_whitney_pvalues <- vapply(
      miRNA_names,
      function(current_miRNA) {
        current_data <- df_long[
          df_long$miRNA == current_miRNA,
          ,
          drop = FALSE
        ]

        test_result <- suppressWarnings(
          stats::wilcox.test(
            RPM ~ Condition,
            data = current_data
          )
        )

        test_result$p.value
      },
      numeric(1)
    )

    results_mann_whitney <- data.frame(
      miRNA = miRNA_names,
      p_value = mann_whitney_pvalues,
      stringsAsFactors = FALSE
    )

    results_mann_whitney$FDR <- stats::p.adjust(
      results_mann_whitney$p_value,
      method = "fdr"
    )

    results_mann_whitney$Significance <- vapply(
      results_mann_whitney$p_value,
      add_significance_symbol,
      character(1)
    )

    results_mann_whitney_sig <- results_mann_whitney[
      !is.na(results_mann_whitney$FDR) &
        results_mann_whitney$FDR < 0.05,
      ,
      drop = FALSE
    ]

    openxlsx::addWorksheet(
      wb,
      "Mann_Whitney"
    )

    openxlsx::writeData(
      wb,
      "Mann_Whitney",
      results_mann_whitney_sig
    )

    comparison_results_list[["Mann_Whitney"]] <-
      results_mann_whitney_sig

    if (isTRUE(assign_results)) {
      assign(
        "mann_whitney_results",
        results_mann_whitney_sig,
        envir = caller_env
      )
    }
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
