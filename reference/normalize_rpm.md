# Normalize a count matrix to RPM without scientific notation

This function normalizes a count matrix to Reads Per Million (RPM) using
the total reads per sample. The resulting matrix is rounded to 2 decimal
places and avoids scientific notation.

## Usage

``` r
normalize_rpm(count_matrix, metrics, reads_column)
```

## Arguments

- count_matrix:

  A count matrix (genes in rows, samples in columns).

- metrics:

  A dataframe containing sample information, including total reads.

- reads_column:

  The name of the column in metrics that contains the total reads.

## Value

A normalized matrix in RPM without scientific notation.

## Examples

``` r
# Example usage:
# count_matrix <- matrix(c(100, 200, 300, 400), nrow = 2, dimnames = list(c("Gene1", "Gene2"), c("Sample1", "Sample2")))
# metrics <- data.frame(Sample = c("Sample1", "Sample2"), TotalReads = c(1000, 2000))
# rpm_matrix <- normalize_rpm(count_matrix, metrics, "TotalReads")
```
