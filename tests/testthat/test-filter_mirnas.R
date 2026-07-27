make_filter_test_data <- function() {
  rpm_matrix <- matrix(
    c(
      5, 6, 0, 0,
      0, 0, 10, 0,
      0, 0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2", "miR-3"),
      c("S1", "S2", "S3", "S4")
    )
  )

  metadata <- data.frame(
    Sample = c("S4", "S2", "S1", "S3"),
    Group = c(
      "Disease",
      "Control",
      "Control",
      "Disease"
    ),
    stringsAsFactors = FALSE
  )

  list(
    rpm_matrix = rpm_matrix,
    metadata = metadata
  )
}

test_that("filter_mirnas retains miRNAs expressed in any group", {
  test_data <- make_filter_test_data()

  result <- filter_mirnas(
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    min_rpm = 5,
    min_count = NULL,
    min_prevalence = 0.5,
    min_samples = NULL,
    group_column = "Group",
    prevalence_scope = "any_group",
    min_groups = 1,
    return_diagnostics = TRUE
  )

  expect_identical(
    rownames(result$filtered_matrix),
    c("miR-1", "miR-2")
  )

  expect_equal(
    result$diagnostics$passes_filter,
    c(TRUE, TRUE, FALSE)
  )

  expect_equal(
    result$diagnostics$prevalence_Control,
    c(1, 0, 0)
  )

  expect_equal(
    result$diagnostics$prevalence_Disease,
    c(0, 0.5, 0)
  )

  expect_equal(
    result$summary$input_mirnas,
    3L
  )

  expect_equal(
    result$summary$retained_mirnas,
    2L
  )

  expect_equal(
    result$summary$removed_mirnas,
    1L
  )

  expect_equal(
    result$summary$retained_fraction,
    2 / 3
  )

  expect_false(
    result$summary$legacy_interface
  )
})

test_that("filter_mirnas supports overall prevalence", {
  rpm_matrix <- matrix(
    c(
      5, 5, 0, 0,
      5, 0, 0, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2", "S3", "S4")
    )
  )

  result <- filter_mirnas(
    rpm_matrix = rpm_matrix,
    min_rpm = 5,
    min_prevalence = 0.5,
    prevalence_scope = "overall",
    return_diagnostics = TRUE
  )

  expect_identical(
    rownames(result$filtered_matrix),
    "miR-1"
  )

  expect_equal(
    result$diagnostics$overall_prevalence,
    c(0.5, 0.25)
  )

  expect_equal(
    result$diagnostics$required_samples_overall,
    c(2, 2)
  )

  expect_equal(
    result$diagnostics$passes_filter,
    c(TRUE, FALSE)
  )
})

test_that("filter_mirnas supports all-groups filtering", {
  rpm_matrix <- matrix(
    c(
      5, 0, 5, 0,
      5, 0, 0, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2", "S3", "S4")
    )
  )

  metadata <- data.frame(
    Sample = c("S1", "S2", "S3", "S4"),
    Group = c(
      "Control",
      "Control",
      "Disease",
      "Disease"
    ),
    stringsAsFactors = FALSE
  )

  result <- filter_mirnas(
    rpm_matrix = rpm_matrix,
    metadata = metadata,
    min_rpm = 5,
    min_prevalence = 0.5,
    group_column = "Group",
    prevalence_scope = "all_groups",
    return_diagnostics = TRUE
  )

  expect_identical(
    rownames(result$filtered_matrix),
    "miR-1"
  )

  expect_equal(
    result$diagnostics$groups_passing,
    c(2, 1)
  )

  expect_equal(
    result$diagnostics$passes_filter,
    c(TRUE, FALSE)
  )
})

test_that("filter_mirnas supports a minimum number of groups", {
  rpm_matrix <- matrix(
    c(
      5, 0, 5, 0,
      5, 0, 0, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2", "S3", "S4")
    )
  )

  metadata <- data.frame(
    Sample = c("S1", "S2", "S3", "S4"),
    Group = c(
      "Control",
      "Control",
      "Disease",
      "Disease"
    ),
    stringsAsFactors = FALSE
  )

  result <- filter_mirnas(
    rpm_matrix = rpm_matrix,
    metadata = metadata,
    min_rpm = 5,
    min_prevalence = 0.5,
    group_column = "Group",
    prevalence_scope = "any_group",
    min_groups = 2,
    return_diagnostics = TRUE
  )

  expect_identical(
    rownames(result$filtered_matrix),
    "miR-1"
  )

  expect_equal(
    result$diagnostics$groups_passing,
    c(2, 1)
  )
})

test_that("filter_mirnas combines prevalence and minimum samples", {
  rpm_matrix <- matrix(
    c(
      5, 5, 0, 0,
      5, 0, 0, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2", "S3", "S4")
    )
  )

  result <- filter_mirnas(
    rpm_matrix = rpm_matrix,
    min_rpm = 5,
    min_prevalence = 0.25,
    min_samples = 2,
    prevalence_scope = "overall",
    return_diagnostics = TRUE
  )

  expect_identical(
    rownames(result$filtered_matrix),
    "miR-1"
  )

  expect_equal(
    result$diagnostics$required_samples_overall,
    c(2, 2)
  )

  expect_equal(
    result$diagnostics$overall_detected_samples,
    c(2, 1)
  )
})

test_that("filter_mirnas applies an explicit raw-count threshold", {
  rpm_matrix <- matrix(
    c(
      5, 6,
      5, 6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2")
    )
  )

  raw_counts <- matrix(
    c(
      20, 15,
      14, 15
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-2", "miR-1"),
      c("S2", "S1")
    )
  )

  result <- filter_mirnas(
    rpm_matrix = rpm_matrix,
    raw_counts = raw_counts,
    min_rpm = 5,
    min_count = 15,
    min_prevalence = 1,
    prevalence_scope = "overall",
    return_diagnostics = TRUE
  )

  expect_identical(
    rownames(result$filtered_matrix),
    "miR-2"
  )

  expect_equal(
    result$diagnostics$overall_detected_samples,
    c(1, 2)
  )

  expect_equal(
    result$summary$count_requirement_source,
    "raw_counts"
  )
})

test_that("filter_mirnas infers counts from normalization metadata", {
  raw_counts <- matrix(
    c(
      10, 10,
      20, 20
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2")
    )
  )

  rpm_matrix <- normalize_rpm(
    count_matrix = raw_counts,
    library_sizes = c(
      S1 = 1e6,
      S2 = 1e6
    )
  )

  result <- filter_mirnas(
    rpm_matrix = rpm_matrix,
    min_rpm = 5,
    min_count = 15,
    min_prevalence = 1,
    prevalence_scope = "overall",
    return_diagnostics = TRUE
  )

  expect_identical(
    rownames(result$filtered_matrix),
    "miR-2"
  )

  expect_equal(
    result$summary$count_requirement_source,
    "inferred_from_normalization"
  )

  expect_equal(
    result$summary$rpm_threshold_read_equivalents,
    c(S1 = 5, S2 = 5)
  )

  normalization_info <- attr(
    result$filtered_matrix,
    "miRPM_normalization",
    exact = TRUE
  )

  expect_equal(
    normalization_info$library_sizes,
    c(S1 = 1e6, S2 = 1e6)
  )
})

test_that("filter_mirnas stores diagnostics as matrix attributes", {
  rpm_matrix <- matrix(
    c(
      5, 5,
      0, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2")
    )
  )

  result <- filter_mirnas(
    rpm_matrix = rpm_matrix,
    min_rpm = 5,
    min_prevalence = 1,
    prevalence_scope = "overall"
  )

  expect_true(
    is.matrix(result)
  )

  expect_identical(
    rownames(result),
    "miR-1"
  )

  expect_true(
    is.data.frame(
      attr(
        result,
        "filter_diagnostics",
        exact = TRUE
      )
    )
  )

  expect_true(
    is.list(
      attr(
        result,
        "filter_summary",
        exact = TRUE
      )
    )
  )

  expect_true(
    is.list(
      attr(
        result,
        "filter_parameters",
        exact = TRUE
      )
    )
  )
})

test_that("filter_mirnas retains the version 0.1.0 interface", {
  rpm_matrix <- matrix(
    c(
      1000, 1100, 0,
      0, 0, 2000,
      1000, 0, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2", "miR-3"),
      c("S1", "S2", "S3")
    )
  )

  metadata <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Condition = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- filter_mirnas(
      rpm_matrix,
      metadata,
      threshold = 1,
      min_reads = 1000,
      threshold_comparison = ">=",
      read_comparison = ">="
    ),
    "deprecated"
  )

  expect_identical(
    rownames(result),
    c("miR-1", "miR-2")
  )

  filter_summary <- attr(
    result,
    "filter_summary",
    exact = TRUE
  )

  expect_true(
    filter_summary$legacy_interface
  )
})

test_that("filter_mirnas validates grouped filtering inputs", {
  test_data <- make_filter_test_data()

  expect_error(
    filter_mirnas(
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      min_rpm = 5,
      min_prevalence = 0.5,
      prevalence_scope = "any_group"
    ),
    "group_column"
  )

  incomplete_metadata <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = test_data$rpm_matrix,
      metadata = incomplete_metadata,
      min_rpm = 5,
      min_prevalence = 0.5,
      group_column = "Group",
      prevalence_scope = "any_group"
    ),
    "S4"
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      min_rpm = 5,
      min_prevalence = 0.5,
      group_column = "Group",
      prevalence_scope = "any_group",
      min_groups = 3
    ),
    "cannot exceed"
  )
})

test_that("filter_mirnas validates abundance and prevalence arguments", {
  rpm_matrix <- matrix(
    c(
      5, 5,
      0, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2")
    )
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = rpm_matrix,
      min_rpm = -1,
      prevalence_scope = "overall"
    ),
    "min_rpm"
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = rpm_matrix,
      min_prevalence = 1.1,
      prevalence_scope = "overall"
    ),
    "min_prevalence"
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = rpm_matrix,
      min_samples = 1.5,
      prevalence_scope = "overall"
    ),
    "positive whole number"
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = rpm_matrix,
      min_count = 15,
      prevalence_scope = "overall"
    ),
    "requires `raw_counts`"
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = rpm_matrix,
      raw_counts = rpm_matrix,
      prevalence_scope = "overall"
    ),
    "only used when"
  )
})

test_that("filter_mirnas rejects invalid expression matrices", {
  negative_matrix <- matrix(
    c(
      5, -1
    ),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2")
    )
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = negative_matrix,
      prevalence_scope = "overall"
    ),
    "negative"
  )

  missing_names <- matrix(
    c(
      5, 5
    ),
    nrow = 1
  )

  expect_error(
    filter_mirnas(
      rpm_matrix = missing_names,
      prevalence_scope = "overall"
    ),
    "row names"
  )
})
