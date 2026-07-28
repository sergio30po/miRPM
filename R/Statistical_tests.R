#' Perform non-parametric differential-expression analysis
#'
#' Performs rank-based differential-expression analysis on a normalized
#' miRNA expression matrix. Two-group comparisons use the
#' Mann-Whitney-Wilcoxon rank-sum test. Comparisons involving three or more
#' groups use the Kruskal-Wallis test followed by Dunn post-hoc comparisons.
#'
#' @param miRNA_ftd A numeric matrix or data frame containing normalized
#'   expression values, with miRNAs in rows and samples in columns.
#' @param metadata A data frame containing sample information. It must include
#'   a `Sample` column and the column named by `condition`.
#' @param condition A single character string naming the metadata column that
#'   defines the comparison groups.
#' @param output_file Optional name of the Excel workbook. When `NULL`, no
#'   workbook is written.
#' @param assign_results Logical. If `TRUE`, significant legacy result objects
#'   are also assigned to the calling environment for backward compatibility.
#'   Default is `TRUE`.
#' @param output_dir Directory used when an Excel workbook is written. Default
#'   is `"Tests_results"`.
#' @param alpha FDR threshold used to define significant results. Default is
#'   `0.05`.
#' @param p_adjust_method Multiple-testing correction passed to
#'   `stats::p.adjust()`. Default is `"BH"`.
#' @param detection_threshold Expression value above which a miRNA is
#'   considered detected when prevalence is calculated. Default is `0`.
#' @param pseudocount Non-negative value added to group medians before
#'   calculating the descriptive log2 median ratio. Default is `1`.
#' @param group_order Optional character vector defining group order. For a
#'   two-group comparison, positive rank-biserial correlation indicates
#'   greater values in the first group than in the second group.
#' @param posthoc_only_if_global_significant Logical. For analyses with three
#'   or more groups, perform Dunn comparisons only for miRNAs with a
#'   Kruskal-Wallis FDR below `alpha`. Default is `TRUE`.
#' @param write_excel Logical. If `TRUE`, write an Excel workbook. By default,
#'   this is `TRUE` when `output_file` is supplied and `FALSE` otherwise.
#'
#' @details
#' The function expects an already normalized and, when appropriate, filtered
#' expression matrix. Filtering and normalization should be completed before
#' differential-expression testing.
#'
#' For two groups, the function reports the Wilcoxon statistic, raw p-value,
#' FDR, rank-biserial correlation, probability of superiority, differences
#' between group means and medians, a descriptive log2 median ratio, and
#' group-specific abundance and prevalence summaries.
#'
#' For three or more groups, the function reports the Kruskal-Wallis
#' statistic, FDR and epsilon-squared effect size. Dunn post-hoc p-values are
#' adjusted across miRNAs separately for each pairwise comparison.
#'
#' The Wilcoxon rank-sum test should not automatically be interpreted as a
#' test of medians. Rank-biserial correlation and probability of superiority
#' describe the direction and magnitude of rank separation between groups.
#'
#' The log2 median ratio is descriptive and is not used as a mandatory
#' significance threshold.
#'
#' @return A structured list containing:
#' \itemize{
#'   \item `analysis_type`: `"two_group"` or `"multi_group"`.
#'   \item `groups`: ordered group names.
#'   \item `group_sizes`: sample counts by group.
#'   \item `results`: complete Mann-Whitney or Kruskal-Wallis results.
#'   \item `significant_results`: rows with FDR below `alpha`.
#'   \item `posthoc_results`: complete Dunn results for multigroup analyses.
#'   \item `significant_posthoc_results`: significant Dunn results.
#'   \item `comparison_results`: legacy-compatible significant result tables.
#'   \item `parameters`: analysis settings.
#'   \item `output_file`: written workbook path, or `NULL`.
#'   \item `output_folder`: workbook directory, or `NULL`.
#' }
#'
#' @examples
#' expression_matrix <- matrix(
#'   c(
#'     1, 2, 1, 8, 9, 10,
#'     5, 4, 6, 5, 4, 6
#'   ),
#'   nrow = 2,
#'   byrow = TRUE,
#'   dimnames = list(
#'     c("miR-1", "miR-2"),
#'     paste0("S", 1:6)
#'   )
#' )
#'
#' metadata <- data.frame(
#'   Sample = paste0("S", 1:6),
#'   Condition = rep(c("Control", "Disease"), each = 3)
#' )
#'
#' result <- perform_statistical_tests(
#'   miRNA_ftd = expression_matrix,
#'   metadata = metadata,
#'   condition = "Condition",
#'   output_file = NULL,
#'   assign_results = FALSE
#' )
#'
#' result$results
#'
#' @export
perform_statistical_tests <- function(
    miRNA_ftd,
    metadata,
    condition,
    output_file = NULL,
    assign_results = TRUE,
    output_dir = "Tests_results",
    alpha = 0.05,
    p_adjust_method = "BH",
    detection_threshold = 0,
    pseudocount = 1,
    group_order = NULL,
    posthoc_only_if_global_significant = TRUE,
    write_excel = !is.null(output_file)
) {
  caller_env <- parent.frame()

  validate_flag <- function(x, argument_name) {
    if (
      !is.logical(x) ||
        length(x) != 1L ||
        is.na(x)
    ) {
      stop(
        paste0("`", argument_name, "` must be `TRUE` or `FALSE`."),
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  validate_number <- function(
      x,
      argument_name,
      minimum = -Inf,
      maximum = Inf
  ) {
    if (
      !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x < minimum ||
        x > maximum
    ) {
      stop(
        paste0(
          "`",
          argument_name,
          "` must be a single finite number between ",
          minimum,
          " and ",
          maximum,
          "."
        ),
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  validate_expression_matrix <- function(x) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }

    if (!is.matrix(x)) {
      stop(
        "`miRNA_ftd` must be a numeric matrix or data frame.",
        call. = FALSE
      )
    }

    if (!is.numeric(x)) {
      stop(
        "`miRNA_ftd` must contain numeric expression values.",
        call. = FALSE
      )
    }

    if (nrow(x) == 0L || ncol(x) == 0L) {
      stop(
        "`miRNA_ftd` must contain at least one row and one column.",
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
        "`miRNA_ftd` must have non-empty miRNA names as row names.",
        call. = FALSE
      )
    }

    if (anyDuplicated(feature_names)) {
      stop(
        "`miRNA_ftd` row names must be unique.",
        call. = FALSE
      )
    }

    if (
      is.null(sample_names) ||
        anyNA(sample_names) ||
        any(sample_names == "")
    ) {
      stop(
        "`miRNA_ftd` must have non-empty sample names as column names.",
        call. = FALSE
      )
    }

    if (anyDuplicated(sample_names)) {
      stop(
        "`miRNA_ftd` column names must be unique.",
        call. = FALSE
      )
    }

    if (anyNA(x) || any(!is.finite(x))) {
      stop(
        "`miRNA_ftd` cannot contain missing or infinite values.",
        call. = FALSE
      )
    }

    if (any(x < 0)) {
      stop(
        "`miRNA_ftd` cannot contain negative expression values.",
        call. = FALSE
      )
    }

    x
  }

  validate_metadata <- function(
      metadata,
      sample_names,
      condition,
      group_order
  ) {
    if (!is.data.frame(metadata)) {
      stop(
        "`metadata` must be a data frame.",
        call. = FALSE
      )
    }

    if (!"Sample" %in% names(metadata)) {
      stop(
        "`metadata` must contain a column named `Sample`.",
        call. = FALSE
      )
    }

    if (
      !is.character(condition) ||
        length(condition) != 1L ||
        is.na(condition) ||
        condition == ""
    ) {
      stop(
        "`condition` must be a single non-empty column name.",
        call. = FALSE
      )
    }

    if (!condition %in% names(metadata)) {
      stop(
        "The condition column was not found in `metadata`.",
        call. = FALSE
      )
    }

    metadata_samples <- as.character(metadata$Sample)

    if (
      anyNA(metadata_samples) ||
        any(metadata_samples == "")
    ) {
      stop(
        "`metadata$Sample` cannot contain missing or empty identifiers.",
        call. = FALSE
      )
    }

    if (anyDuplicated(metadata_samples)) {
      stop(
        "`metadata$Sample` contains duplicated sample identifiers.",
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

    raw_groups <- metadata[[condition]]
    group_values <- as.character(raw_groups)

    if (
      anyNA(group_values) ||
        any(group_values == "")
    ) {
      stop(
        "The condition column cannot contain missing or empty values.",
        call. = FALSE
      )
    }

    observed_groups <- unique(group_values)

    if (is.null(group_order)) {
      if (is.factor(raw_groups)) {
        group_levels <- levels(raw_groups)
        group_levels <- group_levels[
          group_levels %in% observed_groups
        ]
      } else {
        group_levels <- observed_groups
      }
    } else {
      if (
        !is.character(group_order) ||
          length(group_order) < 2L ||
          anyNA(group_order) ||
          any(group_order == "") ||
          anyDuplicated(group_order)
      ) {
        stop(
          paste0(
            "`group_order` must be a character vector of unique, ",
            "non-empty group names."
          ),
          call. = FALSE
        )
      }

      if (
        length(group_order) != length(observed_groups) ||
          !setequal(group_order, observed_groups)
      ) {
        stop(
          paste0(
            "`group_order` must contain every observed group exactly once. ",
            "Observed groups: ",
            paste(observed_groups, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      group_levels <- group_order
    }

    if (length(group_levels) < 2L) {
      stop(
        "The condition must contain at least two groups.",
        call. = FALSE
      )
    }

    group_factor <- factor(
      group_values,
      levels = group_levels
    )

    group_sizes <- table(group_factor)

    small_groups <- names(group_sizes)[group_sizes < 2L]

    if (length(small_groups) > 0L) {
      stop(
        paste0(
          "Every group must contain at least two samples. ",
          "Groups with fewer than two samples: ",
          paste(small_groups, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    list(
      metadata = metadata,
      group_factor = group_factor,
      group_levels = group_levels,
      group_sizes = group_sizes
    )
  }

  significance_symbol <- function(p_value) {
    if (is.na(p_value)) {
      return("")
    }

    if (p_value < 0.001) {
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

  probability_effect <- function(x, y) {
    n_x <- length(x)
    n_y <- length(y)

    combined_ranks <- rank(
      c(x, y),
      ties.method = "average"
    )

    u_x <- sum(combined_ranks[seq_len(n_x)]) -
      n_x * (n_x + 1) / 2

    probability_superiority <- u_x / (n_x * n_y)
    rank_biserial <- 2 * probability_superiority - 1

    list(
      mann_whitney_u = unname(u_x),
      probability_superiority = unname(probability_superiority),
      rank_biserial = unname(rank_biserial)
    )
  }

  make_group_summary <- function(
      values,
      group_factor,
      group_levels,
      safe_group_names,
      detection_threshold
  ) {
    result <- list(
      overall_mean = mean(values),
      overall_median = stats::median(values),
      overall_iqr = stats::IQR(values),
      overall_mad = stats::mad(values),
      overall_minimum = min(values),
      overall_maximum = max(values),
      overall_detected_samples = sum(
        values > detection_threshold
      ),
      overall_prevalence = mean(
        values > detection_threshold
      ),
      overall_zero_fraction = mean(values == 0)
    )

    for (current_group in group_levels) {
      safe_group <- safe_group_names[[current_group]]
      current_values <- values[
        group_factor == current_group
      ]

      result[[paste0("n_", safe_group)]] <-
        length(current_values)

      result[[paste0("mean_", safe_group)]] <-
        mean(current_values)

      result[[paste0("median_", safe_group)]] <-
        stats::median(current_values)

      result[[paste0("iqr_", safe_group)]] <-
        stats::IQR(current_values)

      result[[paste0("mad_", safe_group)]] <-
        stats::mad(current_values)

      result[[paste0("minimum_", safe_group)]] <-
        min(current_values)

      result[[paste0("maximum_", safe_group)]] <-
        max(current_values)

      result[[paste0("detected_samples_", safe_group)]] <-
        sum(current_values > detection_threshold)

      result[[paste0("prevalence_", safe_group)]] <-
        mean(current_values > detection_threshold)

      result[[paste0("zero_fraction_", safe_group)]] <-
        mean(current_values == 0)
    }

    result
  }

  safe_sheet_name <- function(x, used_names = character()) {
    x <- gsub(
      "[:\\\\/?*\\[\\]]",
      "_",
      x
    )

    x <- substr(x, 1L, 31L)

    if (x == "") {
      x <- "Results"
    }

    candidate <- x
    suffix <- 1L

    while (candidate %in% used_names) {
      suffix_text <- paste0("_", suffix)
      candidate <- paste0(
        substr(
          x,
          1L,
          31L - nchar(suffix_text)
        ),
        suffix_text
      )
      suffix <- suffix + 1L
    }

    candidate
  }

  create_parameter_table <- function(
      analysis_type,
      group_levels,
      group_sizes,
      alpha,
      p_adjust_method,
      detection_threshold,
      pseudocount,
      posthoc_only_if_global_significant
  ) {
    data.frame(
      parameter = c(
        "analysis_type",
        "groups",
        "group_sizes",
        "alpha",
        "p_adjust_method",
        "detection_threshold",
        "pseudocount",
        "posthoc_only_if_global_significant"
      ),
      value = c(
        analysis_type,
        paste(group_levels, collapse = " | "),
        paste(
          paste0(
            names(group_sizes),
            "=",
            as.integer(group_sizes)
          ),
          collapse = " | "
        ),
        as.character(alpha),
        p_adjust_method,
        as.character(detection_threshold),
        as.character(pseudocount),
        as.character(posthoc_only_if_global_significant)
      ),
      stringsAsFactors = FALSE
    )
  }

  write_results_workbook <- function(
      analysis_type,
      results,
      significant_results,
      posthoc_results,
      significant_posthoc_results,
      parameter_table,
      output_file,
      output_dir
  ) {
    if (
      !is.character(output_file) ||
        length(output_file) != 1L ||
        is.na(output_file) ||
        output_file == ""
    ) {
      stop(
        paste0(
          "`output_file` must be a single non-empty file name when ",
          "`write_excel = TRUE`."
        ),
        call. = FALSE
      )
    }

    if (
      !is.character(output_dir) ||
        length(output_dir) != 1L ||
        is.na(output_dir) ||
        output_dir == ""
    ) {
      stop(
        "`output_dir` must be a single non-empty directory path.",
        call. = FALSE
      )
    }

    if (!grepl("\\.xlsx$", output_file, ignore.case = TRUE)) {
      output_file <- paste0(output_file, ".xlsx")
    }

    if (!dir.exists(output_dir)) {
      dir.create(
        output_dir,
        recursive = TRUE
      )
    }

    output_path <- file.path(
      output_dir,
      output_file
    )

    workbook <- openxlsx::createWorkbook()
    used_sheet_names <- character()

    add_table <- function(sheet_name, data) {
      safe_name <- safe_sheet_name(
        sheet_name,
        used_sheet_names
      )

      openxlsx::addWorksheet(
        workbook,
        safe_name
      )

      openxlsx::writeData(
        workbook,
        safe_name,
        data,
        withFilter = nrow(data) > 0L
      )

      openxlsx::freezePane(
        workbook,
        safe_name,
        firstRow = TRUE
      )

      used_sheet_names <<- c(
        used_sheet_names,
        safe_name
      )

      invisible(safe_name)
    }

    add_table(
      "Analysis_parameters",
      parameter_table
    )

    if (analysis_type == "two_group") {
      add_table(
        "Mann_Whitney_all",
        results
      )

      add_table(
        "Mann_Whitney_significant",
        significant_results
      )
    } else {
      add_table(
        "Kruskal_Wallis_all",
        results
      )

      add_table(
        "Kruskal_Wallis_significant",
        significant_results
      )

      add_table(
        "Dunn_all",
        posthoc_results
      )

      add_table(
        "Dunn_significant",
        significant_posthoc_results
      )
    }

    openxlsx::saveWorkbook(
      workbook,
      output_path,
      overwrite = TRUE
    )

    normalizePath(
      output_path,
      winslash = "/",
      mustWork = FALSE
    )
  }

  validate_flag(
    assign_results,
    "assign_results"
  )

  validate_flag(
    posthoc_only_if_global_significant,
    "posthoc_only_if_global_significant"
  )

  validate_flag(
    write_excel,
    "write_excel"
  )

  validate_number(
    alpha,
    "alpha",
    minimum = 0,
    maximum = 1
  )

  if (alpha == 0) {
    stop(
      "`alpha` must be greater than zero.",
      call. = FALSE
    )
  }

  validate_number(
    detection_threshold,
    "detection_threshold",
    minimum = 0
  )

  validate_number(
    pseudocount,
    "pseudocount",
    minimum = 0
  )

  if (
    !is.character(p_adjust_method) ||
      length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) ||
      !p_adjust_method %in% stats::p.adjust.methods
  ) {
    stop(
      paste0(
        "`p_adjust_method` must be one of: ",
        paste(stats::p.adjust.methods, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  rpm_matrix <- validate_expression_matrix(
    miRNA_ftd
  )

  metadata_info <- validate_metadata(
    metadata = metadata,
    sample_names = colnames(rpm_matrix),
    condition = condition,
    group_order = group_order
  )

  aligned_metadata <- metadata_info$metadata
  group_factor <- metadata_info$group_factor
  group_levels <- metadata_info$group_levels
  group_sizes <- metadata_info$group_sizes
  number_of_groups <- length(group_levels)

  safe_group_names <- make.names(
    group_levels,
    unique = TRUE
  )

  names(safe_group_names) <- group_levels

  miRNA_names <- rownames(rpm_matrix)

  summary_rows <- lapply(
    seq_len(nrow(rpm_matrix)),
    function(index) {
      summary_values <- make_group_summary(
        values = rpm_matrix[index, ],
        group_factor = group_factor,
        group_levels = group_levels,
        safe_group_names = safe_group_names,
        detection_threshold = detection_threshold
      )

      as.data.frame(
        c(
          list(miRNA = miRNA_names[index]),
          summary_values
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  expression_summary <- do.call(
    rbind,
    summary_rows
  )

  rownames(expression_summary) <- NULL

  comparison_results_list <- list()
  posthoc_results <- data.frame()
  significant_posthoc_results <- data.frame()

  if (number_of_groups == 2L) {
    group_1 <- group_levels[1L]
    group_2 <- group_levels[2L]

    test_rows <- lapply(
      seq_len(nrow(rpm_matrix)),
      function(index) {
        values <- rpm_matrix[index, ]

        values_1 <- values[
          group_factor == group_1
        ]

        values_2 <- values[
          group_factor == group_2
        ]

        effect <- probability_effect(
          values_1,
          values_2
        )

        if (length(unique(values)) == 1L) {
          wilcox_statistic <- effect$mann_whitney_u
          p_value <- 1
        } else {
          test_result <- suppressWarnings(
            stats::wilcox.test(
              x = values_1,
              y = values_2,
              alternative = "two.sided",
              exact = FALSE,
              correct = FALSE
            )
          )

          wilcox_statistic <- unname(
            test_result$statistic
          )

          p_value <- test_result$p.value

          if (is.na(p_value) || !is.finite(p_value)) {
            p_value <- 1
          }
        }

        median_1 <- stats::median(values_1)
        median_2 <- stats::median(values_2)
        mean_1 <- mean(values_1)
        mean_2 <- mean(values_2)

        log2_median_ratio <- if (
          median_1 + pseudocount == 0 &&
            median_2 + pseudocount == 0
        ) {
          0
        } else if (
          median_2 + pseudocount == 0
        ) {
          Inf
        } else {
          log2(
            (median_1 + pseudocount) /
              (median_2 + pseudocount)
          )
        }

        data.frame(
          miRNA = miRNA_names[index],
          group_1 = group_1,
          group_2 = group_2,
          wilcox_statistic = wilcox_statistic,
          mann_whitney_u_group_1 =
            effect$mann_whitney_u,
          p_value = p_value,
          rank_biserial_group_1_vs_group_2 =
            effect$rank_biserial,
          probability_group_1_greater_group_2 =
            effect$probability_superiority,
          median_difference_group_1_minus_group_2 =
            median_1 - median_2,
          mean_difference_group_1_minus_group_2 =
            mean_1 - mean_2,
          log2_median_ratio_group_1_vs_group_2 =
            log2_median_ratio,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )

    results <- do.call(
      rbind,
      test_rows
    )

    results$FDR <- stats::p.adjust(
      results$p_value,
      method = p_adjust_method
    )

    results$p_significance <- vapply(
      results$p_value,
      significance_symbol,
      character(1)
    )

    results$FDR_significance <- vapply(
      results$FDR,
      significance_symbol,
      character(1)
    )

    summary_columns <- setdiff(
      names(expression_summary),
      "miRNA"
    )

    results <- cbind(
      results,
      expression_summary[
        match(results$miRNA, expression_summary$miRNA),
        summary_columns,
        drop = FALSE
      ]
    )

    results <- results[
      order(
        results$FDR,
        results$p_value,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    rownames(results) <- NULL

    significant_results <- results[
      !is.na(results$FDR) &
        results$FDR < alpha,
      ,
      drop = FALSE
    ]

    comparison_results_list[["Mann_Whitney"]] <-
      significant_results

    if (isTRUE(assign_results)) {
      assign(
        "mann_whitney_results",
        significant_results,
        envir = caller_env
      )
    }

    analysis_type <- "two_group"
  } else {
    global_rows <- lapply(
      seq_len(nrow(rpm_matrix)),
      function(index) {
        values <- rpm_matrix[index, ]

        if (length(unique(values)) == 1L) {
          statistic <- 0
          degrees_freedom <- number_of_groups - 1L
          p_value <- 1
        } else {
          test_result <- stats::kruskal.test(
            x = values,
            g = group_factor
          )

          statistic <- unname(
            test_result$statistic
          )

          degrees_freedom <- unname(
            test_result$parameter
          )

          p_value <- test_result$p.value

          if (is.na(p_value) || !is.finite(p_value)) {
            p_value <- 1
          }
        }

        denominator <- length(values) -
          number_of_groups

        epsilon_squared <- if (denominator <= 0) {
          NA_real_
        } else {
          max(
            0,
            (
              statistic -
                number_of_groups +
                1
            ) / denominator
          )
        }

        data.frame(
          miRNA = miRNA_names[index],
          kruskal_wallis_statistic = statistic,
          degrees_freedom = degrees_freedom,
          p_value = p_value,
          epsilon_squared = epsilon_squared,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )

    results <- do.call(
      rbind,
      global_rows
    )

    results$FDR <- stats::p.adjust(
      results$p_value,
      method = p_adjust_method
    )

    results$p_significance <- vapply(
      results$p_value,
      significance_symbol,
      character(1)
    )

    results$FDR_significance <- vapply(
      results$FDR,
      significance_symbol,
      character(1)
    )

    summary_columns <- setdiff(
      names(expression_summary),
      "miRNA"
    )

    results <- cbind(
      results,
      expression_summary[
        match(results$miRNA, expression_summary$miRNA),
        summary_columns,
        drop = FALSE
      ]
    )

    results <- results[
      order(
        results$FDR,
        results$p_value,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    rownames(results) <- NULL

    significant_results <- results[
      !is.na(results$FDR) &
        results$FDR < alpha,
      ,
      drop = FALSE
    ]

    comparison_results_list[["Kruskal_Wallis"]] <-
      significant_results

    if (isTRUE(assign_results)) {
      assign(
        "kruskal_wallis_results",
        significant_results,
        envir = caller_env
      )
    }

    posthoc_mirnas <- if (
      isTRUE(posthoc_only_if_global_significant)
    ) {
      significant_results$miRNA
    } else {
      results$miRNA
    }

    group_pairs <- utils::combn(
      group_levels,
      2L,
      simplify = FALSE
    )

    posthoc_rows <- list()
    posthoc_index <- 1L

    for (current_miRNA in posthoc_mirnas) {
      matrix_index <- match(
        current_miRNA,
        miRNA_names
      )

      current_values <- rpm_matrix[
        matrix_index,
      ]

      dunn_result <- NULL

      invisible(
        utils::capture.output(
          dunn_result <- dunn.test::dunn.test(
            x = current_values,
            g = group_factor,
            method = "none",
            kw = FALSE,
            label = FALSE,
            wrap = FALSE,
            table = FALSE,
            list = TRUE,
            rmc = FALSE,
            alpha = alpha,
            altp = TRUE
          )
        )
      )

      if (
        length(dunn_result$Z) != length(group_pairs) ||
          length(dunn_result$altP) != length(group_pairs)
      ) {
        stop(
          paste0(
            "Unexpected Dunn-test output for miRNA `",
            current_miRNA,
            "`."
          ),
          call. = FALSE
        )
      }

      global_row <- results[
        results$miRNA == current_miRNA,
        ,
        drop = FALSE
      ]

      for (pair_index in seq_along(group_pairs)) {
        current_pair <- group_pairs[[pair_index]]
        group_1 <- current_pair[1L]
        group_2 <- current_pair[2L]

        values_1 <- current_values[
          group_factor == group_1
        ]

        values_2 <- current_values[
          group_factor == group_2
        ]

        effect <- probability_effect(
          values_1,
          values_2
        )

        median_1 <- stats::median(values_1)
        median_2 <- stats::median(values_2)
        mean_1 <- mean(values_1)
        mean_2 <- mean(values_2)

        log2_median_ratio <- if (
          median_1 + pseudocount == 0 &&
            median_2 + pseudocount == 0
        ) {
          0
        } else if (
          median_2 + pseudocount == 0
        ) {
          Inf
        } else {
          log2(
            (median_1 + pseudocount) /
              (median_2 + pseudocount)
          )
        }

        posthoc_rows[[posthoc_index]] <- data.frame(
          miRNA = current_miRNA,
          comparison = paste(
            group_1,
            "vs",
            group_2
          ),
          group_1 = group_1,
          group_2 = group_2,
          dunn_z = unname(
            dunn_result$Z[pair_index]
          ),
          p_value = unname(
            dunn_result$altP[pair_index]
          ),
          global_p_value = global_row$p_value,
          global_FDR = global_row$FDR,
          mann_whitney_u_group_1 =
            effect$mann_whitney_u,
          rank_biserial_group_1_vs_group_2 =
            effect$rank_biserial,
          probability_group_1_greater_group_2 =
            effect$probability_superiority,
          median_difference_group_1_minus_group_2 =
            median_1 - median_2,
          mean_difference_group_1_minus_group_2 =
            mean_1 - mean_2,
          log2_median_ratio_group_1_vs_group_2 =
            log2_median_ratio,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        posthoc_index <- posthoc_index + 1L
      }
    }

    if (length(posthoc_rows) > 0L) {
      posthoc_results <- do.call(
        rbind,
        posthoc_rows
      )

      posthoc_results$FDR <- NA_real_

      comparison_names <- unique(
        posthoc_results$comparison
      )

      for (comparison_name in comparison_names) {
        comparison_indices <- which(
          posthoc_results$comparison ==
            comparison_name
        )

        posthoc_results$FDR[
          comparison_indices
        ] <- stats::p.adjust(
          posthoc_results$p_value[
            comparison_indices
          ],
          method = p_adjust_method
        )
      }

      posthoc_results$p_significance <- vapply(
        posthoc_results$p_value,
        significance_symbol,
        character(1)
      )

      posthoc_results$FDR_significance <- vapply(
        posthoc_results$FDR,
        significance_symbol,
        character(1)
      )

      posthoc_results <- posthoc_results[
        order(
          posthoc_results$comparison,
          posthoc_results$FDR,
          posthoc_results$p_value,
          na.last = TRUE
        ),
        ,
        drop = FALSE
      ]

      rownames(posthoc_results) <- NULL

      significant_posthoc_results <- posthoc_results[
        !is.na(posthoc_results$FDR) &
          posthoc_results$FDR < alpha,
        ,
        drop = FALSE
      ]

      for (comparison_name in comparison_names) {
        comparison_significant <- significant_posthoc_results[
          significant_posthoc_results$comparison ==
            comparison_name,
          ,
          drop = FALSE
        ]

        comparison_results_list[[comparison_name]] <-
          comparison_significant

        if (isTRUE(assign_results)) {
          comparison_rows <- posthoc_results[
            posthoc_results$comparison ==
              comparison_name,
            ,
            drop = FALSE
          ]

          group_1 <- comparison_rows$group_1[1L]
          group_2 <- comparison_rows$group_2[1L]

          object_name <- make.names(
            paste0(
              "dunn_test_",
              group_1,
              "_vs_",
              group_2
            )
          )

          assign(
            object_name,
            comparison_significant,
            envir = caller_env
          )
        }
      }
    } else {
      posthoc_results <- data.frame(
        miRNA = character(),
        comparison = character(),
        group_1 = character(),
        group_2 = character(),
        dunn_z = numeric(),
        p_value = numeric(),
        global_p_value = numeric(),
        global_FDR = numeric(),
        mann_whitney_u_group_1 = numeric(),
        rank_biserial_group_1_vs_group_2 = numeric(),
        probability_group_1_greater_group_2 = numeric(),
        median_difference_group_1_minus_group_2 = numeric(),
        mean_difference_group_1_minus_group_2 = numeric(),
        log2_median_ratio_group_1_vs_group_2 = numeric(),
        FDR = numeric(),
        p_significance = character(),
        FDR_significance = character(),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      significant_posthoc_results <- posthoc_results
    }

    analysis_type <- "multi_group"
  }

  parameter_table <- create_parameter_table(
    analysis_type = analysis_type,
    group_levels = group_levels,
    group_sizes = group_sizes,
    alpha = alpha,
    p_adjust_method = p_adjust_method,
    detection_threshold = detection_threshold,
    pseudocount = pseudocount,
    posthoc_only_if_global_significant =
      posthoc_only_if_global_significant
  )

  output_path <- NULL
  output_folder <- NULL

  if (isTRUE(write_excel)) {
    output_path <- write_results_workbook(
      analysis_type = analysis_type,
      results = results,
      significant_results = significant_results,
      posthoc_results = posthoc_results,
      significant_posthoc_results =
        significant_posthoc_results,
      parameter_table = parameter_table,
      output_file = output_file,
      output_dir = output_dir
    )

    output_folder <- dirname(
      output_path
    )
  }

  parameters <- list(
    condition = condition,
    alpha = alpha,
    p_adjust_method = p_adjust_method,
    detection_threshold = detection_threshold,
    pseudocount = pseudocount,
    group_order = group_levels,
    posthoc_only_if_global_significant =
      posthoc_only_if_global_significant,
    assign_results = assign_results,
    write_excel = write_excel
  )

  list(
    analysis_type = analysis_type,
    groups = group_levels,
    group_sizes = stats::setNames(
      as.integer(group_sizes),
      names(group_sizes)
    ),
    results = results,
    significant_results = significant_results,
    posthoc_results = posthoc_results,
    significant_posthoc_results =
      significant_posthoc_results,
    comparison_results = comparison_results_list,
    parameters = parameters,
    output_file = output_path,
    output_folder = output_folder,
    aligned_metadata = aligned_metadata
  )
}
