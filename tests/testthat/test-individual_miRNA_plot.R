make_individual_plot_data <- function() {
  rpm_matrix <- matrix(
    c(
      5, 6, 10, 11,
      2, 3, 8, 9
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
      "S3",
      "S1",
      "S4",
      "S2"
    ),
    Condition = c(
      "Disease",
      "Control",
      "Disease",
      "Control"
    ),
    Sex = c(
      "F",
      "M",
      "M",
      "F"
    ),
    Age = c(
      72,
      65,
      70,
      68
    ),
    stringsAsFactors = FALSE
  )

  filtered_results <- data.frame(
    miRNA = c(
      "miR-1",
      "miR-2"
    ),
    FDR = c(
      0.01,
      0.04
    ),
    stringsAsFactors = FALSE
  )

  list(
    rpm_matrix = rpm_matrix,
    metadata = metadata,
    filtered_results = filtered_results
  )
}

test_that("individual plots are returned without writing files", {
  test_data <- make_individual_plot_data()

  result <- individual_miRNA_plot(
    filtered_results = test_data$filtered_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups_to_include = c(
      "Control",
      "Disease"
    ),
    condition_colors = c(
      Control = "blue",
      Disease = "red"
    ),
    adjusted_pvalue_column = "FDR",
    save_png = FALSE
  )

  expect_identical(
    names(result$plots),
    c(
      "miR-1",
      "miR-2"
    )
  )

  expect_s3_class(
    result$plots[["miR-1"]],
    "ggplot"
  )

  expect_s3_class(
    result$plots[["miR-2"]],
    "ggplot"
  )

  expect_length(
    result$files,
    0L
  )

  expect_identical(
    result$selected_mirnas,
    c(
      "miR-1",
      "miR-2"
    )
  )

  expect_identical(
    result$missing_mirnas,
    character()
  )

  expect_identical(
    result$groups,
    c(
      "Control",
      "Disease"
    )
  )

  expect_null(
    result$output_folder
  )
})

test_that("metadata are aligned by sample identifier", {
  test_data <- make_individual_plot_data()

  result <- individual_miRNA_plot(
    filtered_results = test_data$filtered_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups_to_include = c(
      "Control",
      "Disease"
    ),
    save_png = FALSE
  )

  expect_identical(
    as.character(result$data$Sample),
    rep(
      c(
        "S1",
        "S2",
        "S3",
        "S4"
      ),
      each = 2
    )
  )

  expect_identical(
    as.character(result$data$Condition),
    rep(
      c(
        "Control",
        "Control",
        "Disease",
        "Disease"
      ),
      each = 2
    )
  )

  expect_equal(
    result$data$RPM,
    c(
      5, 2,
      6, 3,
      10, 8,
      11, 9
    )
  )

  expect_identical(
    as.character(result$data$Sex),
    rep(
      c(
        "M",
        "F",
        "F",
        "M"
      ),
      each = 2
    )
  )

  expect_equal(
    result$data$Age,
    rep(
      c(
        65,
        68,
        72,
        70
      ),
      each = 2
    )
  )
})

test_that("group and miRNA order are preserved", {
  test_data <- make_individual_plot_data()

  reversed_results <- test_data$filtered_results[
    2:1,
    ,
    drop = FALSE
  ]

  result <- individual_miRNA_plot(
    filtered_results = reversed_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups_to_include = c(
      "Disease",
      "Control"
    ),
    save_png = FALSE
  )

  expect_identical(
    result$selected_mirnas,
    c(
      "miR-2",
      "miR-1"
    )
  )

  expect_identical(
    names(result$plots),
    c(
      "miR-2",
      "miR-1"
    )
  )

  expect_identical(
    result$groups,
    c(
      "Disease",
      "Control"
    )
  )

  expect_identical(
    levels(
      result$plots[["miR-2"]]$
        data$.miRPM_group
    ),
    c(
      "Disease",
      "Control"
    )
  )
})

test_that("adjusted p-values are preserved and shown in captions", {
  test_data <- make_individual_plot_data()

  result <- individual_miRNA_plot(
    filtered_results = test_data$filtered_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    adjusted_pvalue_column = "FDR",
    save_png = FALSE
  )

  expect_equal(
    result$adjusted_pvalues,
    c(
      "miR-1" = 0.01,
      "miR-2" = 0.04
    )
  )

  expect_identical(
    result$plots[["miR-1"]]$
      labels$caption,
    "FDR: 0.01"
  )

  expect_identical(
    result$plots[["miR-2"]]$
      labels$caption,
    "FDR: 0.04"
  )
})

test_that("direct miRNA input does not require statistical results", {
  test_data <- make_individual_plot_data()

  result <- individual_miRNA_plot(
    mirnas = c(
      "miR-2",
      "miR-1"
    ),
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    save_png = FALSE
  )

  expect_identical(
    result$selected_mirnas,
    c(
      "miR-2",
      "miR-1"
    )
  )

  expect_length(
    result$adjusted_pvalues,
    0L
  )

  expect_null(
    result$plots[["miR-1"]]$
      labels$caption
  )
})

test_that("custom and legacy miRNA columns are supported", {
  test_data <- make_individual_plot_data()

  custom_results <- data.frame(
    feature = c(
      "miR-2",
      "miR-1"
    ),
    q_value = c(
      0.02,
      0.03
    ),
    stringsAsFactors = FALSE
  )

  custom <- individual_miRNA_plot(
    filtered_results = custom_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    mirna_column = "feature",
    adjusted_pvalue_column = "q_value",
    save_png = FALSE
  )

  expect_identical(
    custom$selected_mirnas,
    c(
      "miR-2",
      "miR-1"
    )
  )

  expect_equal(
    custom$adjusted_pvalues,
    c(
      "miR-2" = 0.02,
      "miR-1" = 0.03
    )
  )

  legacy_first_column <- individual_miRNA_plot(
    filtered_results = custom_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    mirna_column = "missing_column",
    adjusted_pvalue_column = "q_value",
    save_png = FALSE
  )

  expect_identical(
    legacy_first_column$selected_mirnas,
    c(
      "miR-2",
      "miR-1"
    )
  )
})

test_that("the default FDR column is optional", {
  test_data <- make_individual_plot_data()

  results_without_fdr <- data.frame(
    miRNA = c(
      "miR-1",
      "miR-2"
    ),
    stringsAsFactors = FALSE
  )

  result <- individual_miRNA_plot(
    filtered_results = results_without_fdr,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    save_png = FALSE
  )

  expect_length(
    result$adjusted_pvalues,
    0L
  )

  expect_null(
    result$plots[["miR-1"]]$
      labels$caption
  )
})

test_that("missing miRNAs are reported and omitted", {
  test_data <- make_individual_plot_data()

  missing_results <- data.frame(
    miRNA = c(
      "miR-1",
      "miR-missing"
    ),
    FDR = c(
      0.01,
      0.02
    ),
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- individual_miRNA_plot(
      filtered_results = missing_results,
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "miR-missing"
  )

  expect_identical(
    result$selected_mirnas,
    "miR-1"
  )

  expect_identical(
    result$missing_mirnas,
    "miR-missing"
  )

  expect_identical(
    names(result$plots),
    "miR-1"
  )

  expect_equal(
    result$adjusted_pvalues,
    c(
      "miR-1" = 0.01
    )
  )
})

test_that("named colors are applied to all selected groups", {
  test_data <- make_individual_plot_data()

  result <- individual_miRNA_plot(
    filtered_results = test_data$filtered_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups_to_include = c(
      "Control",
      "Disease"
    ),
    condition_colors = c(
      Disease = "red",
      Control = "blue"
    ),
    save_png = FALSE
  )

  color_scale <- result$plots[["miR-1"]]$
    scales$get_scales("colour")

  expect_false(
    is.null(color_scale)
  )

  expect_identical(
    color_scale$breaks,
    c(
      "Control",
      "Disease"
    )
  )
})

test_that("PNG files are optional and use safe names", {
  test_data <- make_individual_plot_data()

  rownames(test_data$rpm_matrix)[1] <-
    "miR:1/alpha?"

  special_results <- data.frame(
    miRNA = "miR:1/alpha?",
    FDR = 0.01,
    stringsAsFactors = FALSE
  )

  output_directory <- tempfile(
    pattern = "miRPM-individual-plots-"
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

  result <- individual_miRNA_plot(
    filtered_results = special_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups_to_include = c(
      "Control",
      "Disease"
    ),
    output_dir = output_directory,
    save_png = TRUE,
    width = 2,
    height = 2,
    dpi = 72
  )

  expect_identical(
    names(result$files),
    "miR:1/alpha?"
  )

  expect_true(
    file.exists(
      result$files[["miR:1/alpha?"]]
    )
  )

  expect_identical(
    basename(
      result$files[["miR:1/alpha?"]]
    ),
    "miR_1_alpha_.png"
  )

  expect_identical(
    basename(result$output_folder),
    "Control_vs_Disease"
  )
})

test_that("the historical positional interface remains supported", {
  test_data <- make_individual_plot_data()

  result <- individual_miRNA_plot(
    test_data$filtered_results,
    test_data$rpm_matrix,
    test_data$metadata,
    "Condition",
    "Sample",
    c(
      "Control",
      "Disease"
    ),
    c(
      Control = "blue",
      Disease = "red"
    ),
    "FDR",
    "unused-output-directory",
    FALSE
  )

  expect_s3_class(
    result$plots[["miR-1"]],
    "ggplot"
  )

  expect_length(
    result$files,
    0L
  )

  expect_null(
    result$output_folder
  )

  expect_equal(
    result$adjusted_pvalues,
    c(
      "miR-1" = 0.01,
      "miR-2" = 0.04
    )
  )
})

test_that("input source arguments are validated", {
  test_data <- make_individual_plot_data()

  expect_error(
    individual_miRNA_plot(
      filtered_results =
        test_data$filtered_results,
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "not both"
  )

  expect_error(
    individual_miRNA_plot(
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "mirnas.*filtered_results"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "rpm_matrix"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "metadata"
  )
})

test_that("statistical result tables are validated", {
  test_data <- make_individual_plot_data()

  duplicated_results <- data.frame(
    miRNA = c(
      "miR-1",
      "miR-1"
    ),
    FDR = c(
      0.01,
      0.02
    ),
    stringsAsFactors = FALSE
  )

  expect_error(
    individual_miRNA_plot(
      filtered_results = duplicated_results,
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "duplicated miRNA"
  )

  invalid_fdr <- test_data$filtered_results
  invalid_fdr$FDR[1] <- 1.2

  expect_error(
    individual_miRNA_plot(
      filtered_results = invalid_fdr,
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "between 0 and 1"
  )

  expect_error(
    individual_miRNA_plot(
      filtered_results =
        test_data$filtered_results,
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      adjusted_pvalue_column =
        "missing_adjusted_pvalue",
      save_png = FALSE
    ),
    "was not found"
  )
})

test_that("matrix, metadata, groups and colors are validated", {
  test_data <- make_individual_plot_data()

  negative_matrix <- test_data$rpm_matrix
  negative_matrix[1, 1] <- -1

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = negative_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "negative"
  )

  incomplete_metadata <- test_data$metadata[
    test_data$metadata$Sample != "S4",
    ,
    drop = FALSE
  ]

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = incomplete_metadata,
      condition_column = "Condition",
      save_png = FALSE
    ),
    "S4"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      groups_to_include = "Unknown",
      save_png = FALSE
    ),
    "Unknown"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      groups_to_include = c(
        "Control",
        "Disease"
      ),
      condition_colors = c(
        Control = "blue"
      ),
      save_png = FALSE
    ),
    "Disease"
  )
})

test_that("plot and output parameters are validated", {
  test_data <- make_individual_plot_data()

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_png = NA
    ),
    "save_png"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      point_size = 0,
      save_png = FALSE
    ),
    "point_size"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      point_alpha = 1.1,
      save_png = FALSE
    ),
    "point_alpha"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      jitter_width = -0.1,
      save_png = FALSE
    ),
    "jitter_width"
  )

  expect_error(
    individual_miRNA_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      output_dir = "",
      save_png = TRUE
    ),
    "output_dir"
  )
})
