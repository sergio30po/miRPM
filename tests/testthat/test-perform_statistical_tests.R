make_two_group_de_data <- function() {
  expression_matrix <- matrix(
    c(
      1, 2, 1, 2, 1, 10, 11, 9, 12, 10,
      5, 4, 6, 5, 4, 5, 4, 6, 5, 4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      paste0("S", 1:10)
    )
  )

  metadata <- data.frame(
    Sample = paste0(
      "S",
      c(6:10, 1:5)
    ),
    Condition = c(
      rep("Disease", 5),
      rep("Control", 5)
    ),
    stringsAsFactors = FALSE
  )

  list(
    expression_matrix = expression_matrix,
    metadata = metadata
  )
}

make_multi_group_de_data <- function() {
  expression_matrix <- matrix(
    c(
      1, 1, 2, 2,
      5, 6, 5, 6,
      10, 11, 10, 11,
      5, 5, 5, 5,
      5, 5, 5, 5,
      5, 5, 5, 5
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      paste0("S", 1:12)
    )
  )

  metadata <- data.frame(
    Sample = paste0(
      "S",
      c(9:12, 1:4, 5:8)
    ),
    Group = c(
      rep("C", 4),
      rep("A", 4),
      rep("B", 4)
    ),
    stringsAsFactors = FALSE
  )

  list(
    expression_matrix = expression_matrix,
    metadata = metadata
  )
}

test_that("two-group analysis returns complete structured results", {
  test_data <- make_two_group_de_data()

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "Control",
      "Disease"
    )
  )

  expect_identical(
    result$analysis_type,
    "two_group"
  )

  expect_identical(
    result$groups,
    c(
      "Control",
      "Disease"
    )
  )

  expect_equal(
    result$group_sizes,
    c(
      Control = 5L,
      Disease = 5L
    )
  )

  expect_equal(
    nrow(result$results),
    2L
  )

  expect_setequal(
    result$results$miRNA,
    c(
      "miR-1",
      "miR-2"
    )
  )

  expected_columns <- c(
    "miRNA",
    "group_1",
    "group_2",
    "wilcox_statistic",
    "mann_whitney_u_group_1",
    "p_value",
    "FDR",
    "rank_biserial_group_1_vs_group_2",
    "probability_group_1_greater_group_2",
    "median_difference_group_1_minus_group_2",
    "mean_difference_group_1_minus_group_2",
    "log2_median_ratio_group_1_vs_group_2",
    "median_Control",
    "median_Disease",
    "prevalence_Control",
    "prevalence_Disease"
  )

  expect_true(
    all(
      expected_columns %in%
        names(result$results)
    )
  )

  expect_null(
    result$output_file
  )

  expect_null(
    result$output_folder
  )

  expect_false(
    result$parameters$write_excel
  )
})

test_that("two-group effect sizes have the expected direction", {
  test_data <- make_two_group_de_data()

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "Control",
      "Disease"
    )
  )

  mir_1 <- result$results[
    result$results$miRNA == "miR-1",
    ,
    drop = FALSE
  ]

  expect_equal(
    mir_1$mann_whitney_u_group_1,
    0
  )

  expect_equal(
    mir_1$probability_group_1_greater_group_2,
    0
  )

  expect_equal(
    mir_1$rank_biserial_group_1_vs_group_2,
    -1
  )

  expect_lt(
    mir_1$median_difference_group_1_minus_group_2,
    0
  )

  expect_lt(
    mir_1$mean_difference_group_1_minus_group_2,
    0
  )

  expect_lt(
    mir_1$log2_median_ratio_group_1_vs_group_2,
    0
  )

  expect_lt(
    mir_1$p_value,
    0.05
  )

  expect_lt(
    mir_1$FDR,
    0.05
  )

  expect_equal(
    mir_1$prevalence_Control,
    1
  )

  expect_equal(
    mir_1$prevalence_Disease,
    1
  )
})

test_that("equivalent two-group distributions have null effect size", {
  test_data <- make_two_group_de_data()

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "Control",
      "Disease"
    )
  )

  mir_2 <- result$results[
    result$results$miRNA == "miR-2",
    ,
    drop = FALSE
  ]

  expect_equal(
    mir_2$probability_group_1_greater_group_2,
    0.5
  )

  expect_equal(
    mir_2$rank_biserial_group_1_vs_group_2,
    0
  )

  expect_equal(
    mir_2$median_difference_group_1_minus_group_2,
    0
  )

  expect_equal(
    mir_2$mean_difference_group_1_minus_group_2,
    0
  )

  expect_equal(
    mir_2$log2_median_ratio_group_1_vs_group_2,
    0
  )

  expect_equal(
    mir_2$p_value,
    1
  )
})

test_that("metadata are aligned by sample identifier", {
  test_data <- make_two_group_de_data()

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "Control",
      "Disease"
    )
  )

  expect_identical(
    as.character(
      result$aligned_metadata$Sample
    ),
    colnames(
      test_data$expression_matrix
    )
  )

  expect_identical(
    as.character(
      result$aligned_metadata$Condition
    ),
    c(
      rep("Control", 5),
      rep("Disease", 5)
    )
  )
})

test_that("group order controls effect-size direction", {
  test_data <- make_two_group_de_data()

  control_first <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "Control",
      "Disease"
    )
  )

  disease_first <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "Disease",
      "Control"
    )
  )

  control_effect <- control_first$results[
    control_first$results$miRNA == "miR-1",
    "rank_biserial_group_1_vs_group_2"
  ]

  disease_effect <- disease_first$results[
    disease_first$results$miRNA == "miR-1",
    "rank_biserial_group_1_vs_group_2"
  ]

  expect_equal(
    control_effect,
    -1
  )

  expect_equal(
    disease_effect,
    1
  )

  expect_equal(
    control_effect,
    -disease_effect
  )
})

test_that("detection threshold controls prevalence summaries", {
  expression_matrix <- matrix(
    c(
      0, 5, 6, 10,
      0, 0, 0, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2", "S3", "S4")
    )
  )

  metadata <- data.frame(
    Sample = c(
      "S1",
      "S2",
      "S3",
      "S4"
    ),
    Condition = c(
      "A",
      "A",
      "B",
      "B"
    ),
    stringsAsFactors = FALSE
  )

  result <- perform_statistical_tests(
    miRNA_ftd = expression_matrix,
    metadata = metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c("A", "B"),
    detection_threshold = 5
  )

  mir_1 <- result$results[
    result$results$miRNA == "miR-1",
    ,
    drop = FALSE
  ]

  expect_equal(
    mir_1$detected_samples_A,
    0
  )

  expect_equal(
    mir_1$prevalence_A,
    0
  )

  expect_equal(
    mir_1$detected_samples_B,
    2
  )

  expect_equal(
    mir_1$prevalence_B,
    1
  )

  expect_equal(
    mir_1$overall_detected_samples,
    2
  )

  expect_equal(
    mir_1$overall_prevalence,
    0.5
  )
})

test_that("multigroup analysis returns Kruskal-Wallis results", {
  test_data <- make_multi_group_de_data()

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Group",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "A",
      "B",
      "C"
    )
  )

  expect_identical(
    result$analysis_type,
    "multi_group"
  )

  expect_identical(
    result$groups,
    c(
      "A",
      "B",
      "C"
    )
  )

  expect_equal(
    result$group_sizes,
    c(
      A = 4L,
      B = 4L,
      C = 4L
    )
  )

  expect_equal(
    nrow(result$results),
    2L
  )

  expected_columns <- c(
    "miRNA",
    "kruskal_wallis_statistic",
    "degrees_freedom",
    "p_value",
    "FDR",
    "epsilon_squared",
    "median_A",
    "median_B",
    "median_C",
    "prevalence_A",
    "prevalence_B",
    "prevalence_C"
  )

  expect_true(
    all(
      expected_columns %in%
        names(result$results)
    )
  )

  mir_1 <- result$results[
    result$results$miRNA == "miR-1",
    ,
    drop = FALSE
  ]

  expect_gt(
    mir_1$kruskal_wallis_statistic,
    0
  )

  expect_equal(
    mir_1$degrees_freedom,
    2
  )

  expect_lt(
    mir_1$p_value,
    0.05
  )

  expect_lt(
    mir_1$FDR,
    0.05
  )

  expect_gt(
    mir_1$epsilon_squared,
    0
  )
})

test_that("multigroup analysis produces Dunn comparisons", {
  test_data <- make_multi_group_de_data()

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Group",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "A",
      "B",
      "C"
    ),
    posthoc_only_if_global_significant = TRUE
  )

  expect_equal(
    nrow(result$posthoc_results),
    3L
  )

  expect_identical(
    unique(result$posthoc_results$miRNA),
    "miR-1"
  )

  expect_setequal(
    result$posthoc_results$comparison,
    c(
      "A vs B",
      "A vs C",
      "B vs C"
    )
  )

  expected_columns <- c(
    "miRNA",
    "comparison",
    "group_1",
    "group_2",
    "dunn_z",
    "p_value",
    "FDR",
    "global_p_value",
    "global_FDR",
    "rank_biserial_group_1_vs_group_2",
    "probability_group_1_greater_group_2",
    "median_difference_group_1_minus_group_2"
  )

  expect_true(
    all(
      expected_columns %in%
        names(result$posthoc_results)
    )
  )

  expect_true(
    all(
      result$posthoc_results$
        rank_biserial_group_1_vs_group_2 < 0
    )
  )

  expect_true(
    all(
      result$posthoc_results$
        median_difference_group_1_minus_group_2 < 0
    )
  )

  expect_true(
    all(
      result$posthoc_results$global_FDR < 0.05
    )
  )
})

test_that("constant multigroup expression returns a null result", {
  expression_matrix <- matrix(
    rep(5, 12),
    nrow = 1,
    dimnames = list(
      "miR-constant",
      paste0("S", 1:12)
    )
  )

  metadata <- data.frame(
    Sample = paste0("S", 1:12),
    Group = rep(
      c("A", "B", "C"),
      each = 4
    ),
    stringsAsFactors = FALSE
  )

  result <- perform_statistical_tests(
    miRNA_ftd = expression_matrix,
    metadata = metadata,
    condition = "Group",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "A",
      "B",
      "C"
    )
  )

  expect_equal(
    result$results$kruskal_wallis_statistic,
    0
  )

  expect_equal(
    result$results$p_value,
    1
  )

  expect_equal(
    result$results$FDR,
    1
  )

  expect_equal(
    result$results$epsilon_squared,
    0
  )

  expect_equal(
    nrow(result$significant_results),
    0L
  )

  expect_equal(
    nrow(result$posthoc_results),
    0L
  )
})

test_that("complete and significant results are returned separately", {
  test_data <- make_two_group_de_data()

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = NULL,
    assign_results = FALSE,
    group_order = c(
      "Control",
      "Disease"
    )
  )

  expect_equal(
    nrow(result$results),
    2L
  )

  expect_true(
    "miR-1" %in%
      result$significant_results$miRNA
  )

  expect_false(
    "miR-2" %in%
      result$significant_results$miRNA
  )

  expect_identical(
    result$comparison_results$Mann_Whitney,
    result$significant_results
  )
})

test_that("legacy result assignment remains available", {
  test_data <- make_two_group_de_data()

  assignment_environment <- new.env(
    parent = environment()
  )

  assignment_environment$expression_matrix <-
    test_data$expression_matrix

  assignment_environment$metadata <-
    test_data$metadata

  evalq(
    result <- perform_statistical_tests(
      miRNA_ftd = expression_matrix,
      metadata = metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = TRUE,
      group_order = c(
        "Control",
        "Disease"
      )
    ),
    envir = assignment_environment
  )

  expect_true(
    exists(
      "mann_whitney_results",
      envir = assignment_environment,
      inherits = FALSE
    )
  )

  assigned_results <- get(
    "mann_whitney_results",
    envir = assignment_environment,
    inherits = FALSE
  )

  expect_identical(
    assigned_results,
    assignment_environment$
      result$significant_results
  )
})

test_that("Excel output is optional and can be written", {
  test_data <- make_two_group_de_data()

  output_directory <- tempfile(
    pattern = "miRPM-statistical-tests-"
  )

  dir.create(
    output_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      output_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  result <- perform_statistical_tests(
    miRNA_ftd = test_data$expression_matrix,
    metadata = test_data$metadata,
    condition = "Condition",
    output_file = "two_group_results",
    output_dir = output_directory,
    assign_results = FALSE,
    group_order = c(
      "Control",
      "Disease"
    ),
    write_excel = TRUE
  )

  expect_true(
    file.exists(result$output_file)
  )

  expect_true(
    grepl(
      "\\.xlsx$",
      result$output_file,
      ignore.case = TRUE
    )
  )

  expect_identical(
    result$output_folder,
    dirname(result$output_file)
  )

  sheet_names <- openxlsx::getSheetNames(
    result$output_file
  )

  expect_setequal(
    sheet_names,
    c(
      "Analysis_parameters",
      "Mann_Whitney_all",
      "Mann_Whitney_significant"
    )
  )
})

test_that("the historical positional interface still works", {
  test_data <- make_two_group_de_data()

  output_directory <- tempfile(
    pattern = "miRPM-legacy-tests-"
  )

  dir.create(
    output_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      output_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  result <- perform_statistical_tests(
    test_data$expression_matrix,
    test_data$metadata,
    "Condition",
    "legacy_results.xlsx",
    FALSE,
    output_dir = output_directory,
    group_order = c(
      "Control",
      "Disease"
    )
  )

  expect_true(
    file.exists(result$output_file)
  )

  expect_identical(
    result$analysis_type,
    "two_group"
  )

  expect_false(
    result$parameters$assign_results
  )
})

test_that("invalid expression matrices are rejected", {
  metadata <- data.frame(
    Sample = c(
      "S1",
      "S2",
      "S3",
      "S4"
    ),
    Condition = c(
      "A",
      "A",
      "B",
      "B"
    ),
    stringsAsFactors = FALSE
  )

  negative_matrix <- matrix(
    c(
      1, -1, 2, 3
    ),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2", "S3", "S4")
    )
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = negative_matrix,
      metadata = metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE
    ),
    "negative"
  )

  missing_matrix <- matrix(
    c(
      1, NA, 2, 3
    ),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2", "S3", "S4")
    )
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = missing_matrix,
      metadata = metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE
    ),
    "missing or infinite"
  )

  unnamed_matrix <- matrix(
    c(
      1, 2, 3, 4
    ),
    nrow = 1
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = unnamed_matrix,
      metadata = metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE
    ),
    "row names"
  )
})

test_that("invalid metadata are rejected", {
  test_data <- make_two_group_de_data()

  missing_sample_metadata <- test_data$metadata[
    test_data$metadata$Sample != "S10",
    ,
    drop = FALSE
  ]

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = missing_sample_metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE
    ),
    "S10"
  )

  duplicated_metadata <- rbind(
    test_data$metadata,
    test_data$metadata[1, , drop = FALSE]
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = duplicated_metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE
    ),
    "duplicated"
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "MissingCondition",
      output_file = NULL,
      assign_results = FALSE
    ),
    "condition column"
  )
})

test_that("groups require valid sample sizes and ordering", {
  expression_matrix <- matrix(
    c(
      1, 2, 3
    ),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2", "S3")
    )
  )

  metadata_small_group <- data.frame(
    Sample = c(
      "S1",
      "S2",
      "S3"
    ),
    Condition = c(
      "A",
      "A",
      "B"
    ),
    stringsAsFactors = FALSE
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = expression_matrix,
      metadata = metadata_small_group,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE
    ),
    "fewer than two"
  )

  test_data <- make_two_group_de_data()

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE,
      group_order = c(
        "Control",
        "Unknown"
      )
    ),
    "every observed group"
  )
})

test_that("analysis parameters are validated", {
  test_data <- make_two_group_de_data()

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE,
      alpha = 0
    ),
    "greater than zero"
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE,
      alpha = 1.1
    ),
    "alpha"
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE,
      detection_threshold = -1
    ),
    "detection_threshold"
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE,
      pseudocount = -1
    ),
    "pseudocount"
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE,
      p_adjust_method = "invalid_method"
    ),
    "p_adjust_method"
  )

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = NA
    ),
    "assign_results"
  )
})

test_that("Excel writing requires a valid output file", {
  test_data <- make_two_group_de_data()

  expect_error(
    perform_statistical_tests(
      miRNA_ftd = test_data$expression_matrix,
      metadata = test_data$metadata,
      condition = "Condition",
      output_file = NULL,
      assign_results = FALSE,
      write_excel = TRUE
    ),
    "output_file"
  )
})
