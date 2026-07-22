expect_rpm_matrix_equal <- function(
    result,
    expected,
    tolerance = 1e-8
) {
  expect_equal(dim(result), dim(expected))
  expect_equal(dimnames(result), dimnames(expected))

  expect_equal(
    as.numeric(result),
    as.numeric(expected),
    tolerance = tolerance
  )
}

test_that("normalize_rpm calculates RPM from the complete count matrix", {
  count_matrix <- matrix(
    c(
      100, 200,
      300, 400
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2")
    )
  )

  result <- normalize_rpm(count_matrix)

  expected <- matrix(
    c(
      250000, 333333.333333333,
      750000, 666666.666666667
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = dimnames(count_matrix)
  )

  expect_rpm_matrix_equal(
    result,
    expected,
    tolerance = 1e-8
  )

  expect_equal(
    colSums(result),
    c(S1 = 1e6, S2 = 1e6)
  )

  normalization_info <- attr(
    result,
    "miRPM_normalization",
    exact = TRUE
  )

  expect_equal(
    normalization_info$method,
    "RPM"
  )

  expect_equal(
    normalization_info$scale,
    1e6
  )

  expect_equal(
    normalization_info$library_sizes,
    c(S1 = 400, S2 = 600)
  )

  expect_equal(
    normalization_info$reads_per_rpm,
    c(S1 = 400, S2 = 600) / 1e6
  )

  expect_equal(
    normalization_info$library_size_source,
    "matrix_column_sums"
  )

  expect_equal(
    normalization_info$input_features,
    2L
  )

  expect_equal(
    normalization_info$input_samples,
    2L
  )
})

test_that("normalize_rpm accepts numeric data frames", {
  count_data <- data.frame(
    S1 = c(100, 300),
    S2 = c(200, 400),
    row.names = c("miR-1", "miR-2")
  )

  result <- normalize_rpm(count_data)

  expect_true(is.matrix(result))

  expect_equal(
    dimnames(result),
    dimnames(as.matrix(count_data))
  )

  expect_equal(
    colSums(result),
    c(S1 = 1e6, S2 = 1e6)
  )

  normalization_info <- attr(
    result,
    "miRPM_normalization",
    exact = TRUE
  )

  expect_equal(
    normalization_info$library_size_source,
    "matrix_column_sums"
  )
})

test_that("normalize_rpm aligns named external library sizes", {
  count_matrix <- matrix(
    c(
      10, 20,
      30, 40
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2")
    )
  )

  library_sizes <- c(
    S2 = 200,
    S1 = 100
  )

  result <- normalize_rpm(
    count_matrix,
    library_sizes = library_sizes
  )

  expected <- sweep(
    count_matrix,
    MARGIN = 2,
    STATS = c(S1 = 100, S2 = 200),
    FUN = "/"
  ) * 1e6

  expect_rpm_matrix_equal(
    result,
    expected
  )

  normalization_info <- attr(
    result,
    "miRPM_normalization",
    exact = TRUE
  )

  expect_equal(
    normalization_info$library_sizes,
    c(S1 = 100, S2 = 200)
  )

  expect_equal(
    normalization_info$reads_per_rpm,
    c(S1 = 100, S2 = 200) / 1e6
  )

  expect_equal(
    normalization_info$library_size_source,
    "external_vector"
  )
})

test_that("normalize_rpm rejects invalid count values", {
  negative_matrix <- matrix(
    c(10, -1),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2")
    )
  )

  missing_matrix <- matrix(
    c(10, NA_real_),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2")
    )
  )

  infinite_matrix <- matrix(
    c(10, Inf),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2")
    )
  )

  expect_error(
    normalize_rpm(negative_matrix),
    "negative"
  )

  expect_error(
    normalize_rpm(missing_matrix),
    "missing or infinite"
  )

  expect_error(
    normalize_rpm(infinite_matrix),
    "missing or infinite"
  )
})

test_that("normalize_rpm rejects samples with zero library size", {
  count_matrix <- matrix(
    c(
      10, 0,
      20, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("miR-1", "miR-2"),
      c("S1", "S2")
    )
  )

  expect_error(
    normalize_rpm(count_matrix),
    "S2"
  )
})

test_that("normalize_rpm detects missing external library sizes", {
  count_matrix <- matrix(
    c(10, 20),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2")
    )
  )

  expect_error(
    normalize_rpm(
      count_matrix,
      library_sizes = c(S1 = 100)
    ),
    "missing samples: S2"
  )
})

test_that("normalize_rpm retains the deprecated version 0.1.0 interface", {
  count_matrix <- matrix(
    c(10, 20),
    nrow = 1,
    dimnames = list(
      "miR-1",
      c("S1", "S2")
    )
  )

  metrics <- data.frame(
    Sample = c("S2", "S1"),
    TotalReads = c(200, 100)
  )

  expect_warning(
    result <- normalize_rpm(
      count_matrix = count_matrix,
      metrics = metrics,
      reads_column = "TotalReads"
    ),
    "deprecated"
  )

  expected <- matrix(
    c(100000, 100000),
    nrow = 1,
    dimnames = dimnames(count_matrix)
  )

  expect_rpm_matrix_equal(
    result,
    expected
  )

  normalization_info <- attr(
    result,
    "miRPM_normalization",
    exact = TRUE
  )

  expect_equal(
    normalization_info$library_size_source,
    "legacy_metrics"
  )

  expect_equal(
    normalization_info$library_sizes,
    c(S1 = 100, S2 = 200)
  )

  expect_equal(
    normalization_info$reads_per_rpm,
    c(S1 = 100, S2 = 200) / 1e6
  )
})
