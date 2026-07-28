make_expression_plot_data <- function() {
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

  list(
    rpm_matrix = rpm_matrix,
    metadata = metadata
  )
}

test_that("miRNA_expression_plot returns static and interactive plots", {
  test_data <- make_expression_plot_data()

  result <- miRNA_expression_plot(
    mirnas = c(
      "miR-1",
      "miR-2"
    ),
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups = c(
      "Control",
      "Disease"
    ),
    tooltip_columns = c(
      "Sex",
      "Age"
    ),
    save_html = FALSE
  )

  expect_s3_class(
    result$ggplot,
    "ggplot"
  )

  expect_s3_class(
    result$interactive,
    "plotly"
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
    result$output_file
  )

  expect_null(
    result$output_folder
  )
})

test_that("metadata are aligned by sample identifier", {
  test_data <- make_expression_plot_data()

  result <- miRNA_expression_plot(
    mirnas = c(
      "miR-1",
      "miR-2"
    ),
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups = c(
      "Control",
      "Disease"
    ),
    save_html = FALSE
  )

  expect_identical(
    as.character(
      result$data$Sample
    ),
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
    as.character(
      result$data$Condition
    ),
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
    as.character(
      result$data$Sex
    ),
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

test_that("miRNA and group order are preserved", {
  test_data <- make_expression_plot_data()

  result <- miRNA_expression_plot(
    mirnas = c(
      "miR-2",
      "miR-1"
    ),
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups = c(
      "Disease",
      "Control"
    ),
    save_html = FALSE
  )

  expect_identical(
    result$selected_mirnas,
    c(
      "miR-2",
      "miR-1"
    )
  )

  expect_identical(
    levels(result$data$miRNA),
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
    levels(result$data$Condition),
    c(
      "Disease",
      "Control"
    )
  )

  expect_identical(
    levels(result$data$.miRPM_group),
    c(
      "Disease",
      "Control"
    )
  )
})

test_that("additional tooltip columns are configurable", {
  test_data <- make_expression_plot_data()

  result <- miRNA_expression_plot(
    mirnas = "miR-1",
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    tooltip_columns = c(
      "Sex",
      "Age"
    ),
    save_html = FALSE
  )

  expect_true(
    all(
      grepl(
        "<b>Sex:</b>",
        result$data$.miRPM_tooltip,
        fixed = TRUE
      )
    )
  )

  expect_true(
    all(
      grepl(
        "<b>Age:</b>",
        result$data$.miRPM_tooltip,
        fixed = TRUE
      )
    )
  )

  expect_true(
    all(
      grepl(
        "<b>Sample:</b>",
        result$data$.miRPM_tooltip,
        fixed = TRUE
      )
    )
  )

  expect_true(
    all(
      grepl(
        "<b>RPM:</b>",
        result$data$.miRPM_tooltip,
        fixed = TRUE
      )
    )
  )
})

test_that("missing miRNAs are reported and omitted", {
  test_data <- make_expression_plot_data()

  expect_warning(
    result <- miRNA_expression_plot(
      mirnas = c(
        "miR-1",
        "miR-missing"
      ),
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_html = FALSE
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

  expect_true(
    all(
      as.character(result$data$miRNA) ==
        "miR-1"
    )
  )
})

test_that("table input for miRNAs remains supported", {
  test_data <- make_expression_plot_data()

  mirna_results <- data.frame(
    miRNA = c(
      "miR-2",
      "miR-1"
    ),
    FDR = c(
      0.01,
      0.02
    ),
    stringsAsFactors = FALSE
  )

  result <- miRNA_expression_plot(
    mirnas = mirna_results,
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    save_html = FALSE
  )

  expect_identical(
    result$selected_mirnas,
    c(
      "miR-2",
      "miR-1"
    )
  )
})

test_that("named colors are applied to all selected groups", {
  test_data <- make_expression_plot_data()

  result <- miRNA_expression_plot(
    mirnas = "miR-1",
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    groups = c(
      "Control",
      "Disease"
    ),
    colors = c(
      Disease = "red",
      Control = "blue"
    ),
    save_html = FALSE
  )

  color_scale <- result$ggplot$scales$get_scales(
    "colour"
  )

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

test_that("HTML output is optional and uses a safe file name", {
  test_data <- make_expression_plot_data()

  output_directory <- tempfile(
    pattern = "miRPM-expression-plot-"
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

  result <- miRNA_expression_plot(
    mirnas = "miR-1",
    rpm_matrix = test_data$rpm_matrix,
    metadata = test_data$metadata,
    condition_column = "Condition",
    output_name = "miRNA:plot?test.html",
    output_dir = output_directory,
    save_html = TRUE,
    selfcontained = FALSE
  )

  expect_true(
    file.exists(result$output_file)
  )

  expect_identical(
    basename(result$output_file),
    "miRNA_plot_test.html"
  )

  expect_identical(
    result$output_folder,
    dirname(result$output_file)
  )
})

test_that("the historical positional interface still works", {
  test_data <- make_expression_plot_data()

  mirna_results <- data.frame(
    miRNA = c(
      "miR-1",
      "miR-2"
    ),
    stringsAsFactors = FALSE
  )

  result <- miRNA_expression_plot(
    mirna_results,
    test_data$rpm_matrix,
    test_data$metadata,
    "Condition",
    c(
      "Control",
      "Disease"
    ),
    c(
      Control = "blue",
      Disease = "red"
    ),
    NULL,
    "Legacy plot",
    save_html = FALSE
  )

  expect_s3_class(
    result$ggplot,
    "ggplot"
  )

  expect_s3_class(
    result$interactive,
    "plotly"
  )

  expect_identical(
    result$selected_mirnas,
    c(
      "miR-1",
      "miR-2"
    )
  )
})

test_that("expression plot validates matrix inputs", {
  test_data <- make_expression_plot_data()

  negative_matrix <- test_data$rpm_matrix
  negative_matrix[1, 1] <- -1

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = negative_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_html = FALSE
    ),
    "negative"
  )

  unnamed_matrix <- unname(
    test_data$rpm_matrix
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = unnamed_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_html = FALSE
    ),
    "row names"
  )
})

test_that("expression plot validates metadata and groups", {
  test_data <- make_expression_plot_data()

  incomplete_metadata <- test_data$metadata[
    test_data$metadata$Sample != "S4",
    ,
    drop = FALSE
  ]

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = incomplete_metadata,
      condition_column = "Condition",
      save_html = FALSE
    ),
    "S4"
  )

  duplicated_metadata <- rbind(
    test_data$metadata,
    test_data$metadata[1, , drop = FALSE]
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = duplicated_metadata,
      condition_column = "Condition",
      save_html = FALSE
    ),
    "duplicated"
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "MissingCondition",
      save_html = FALSE
    ),
    "Missing metadata columns"
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      groups = "Unknown",
      save_html = FALSE
    ),
    "Unknown"
  )
})

test_that("expression plot validates colors and tooltips", {
  test_data <- make_expression_plot_data()

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      groups = c(
        "Control",
        "Disease"
      ),
      colors = c(
        Control = "blue"
      ),
      save_html = FALSE
    ),
    "Disease"
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      tooltip_columns = "MissingTooltip",
      save_html = FALSE
    ),
    "MissingTooltip"
  )
})

test_that("expression plot validates aliases and output arguments", {
  test_data <- make_expression_plot_data()

  expect_error(
    miRNA_expression_plot(
      miRNAs_DE = "miR-1",
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_html = FALSE
    ),
    "not both"
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      miRNA_ftd = test_data$rpm_matrix,
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_html = FALSE
    ),
    "not both"
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_html = TRUE
    ),
    "output_name"
  )

  expect_error(
    miRNA_expression_plot(
      mirnas = "miR-1",
      rpm_matrix = test_data$rpm_matrix,
      metadata = test_data$metadata,
      condition_column = "Condition",
      save_html = NA
    ),
    "save_html"
  )
})
