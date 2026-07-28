# Normalize a count matrix to reads per million

Converts a complete count matrix to reads per million (RPM) using the
library size of each sample.

## Usage

``` r
normalize_rpm(
  count_matrix,
  metrics = NULL,
  reads_column = NULL,
  library_sizes = NULL
)
```

## Arguments

- count_matrix:

  A numeric matrix or data frame with features in rows and samples in
  columns.

- metrics:

  Deprecated. A data frame containing a \`Sample\` column and the
  library-size column indicated by \`reads_column\`.

- reads_column:

  Deprecated. Name of the library-size column in \`metrics\`.

- library_sizes:

  Optional named numeric vector containing one positive library size per
  sample. Names must match the column names of \`count_matrix\`.

## Value

A numeric matrix with the same dimensions and dimnames as
\`count_matrix\`. The \`"miRPM_normalization"\` attribute records the
normalization method, library sizes, their source and reads per RPM.

## Details

When \`library_sizes\` is \`NULL\`, library sizes are calculated from
the complete input matrix using \`colSums(count_matrix)\`.

Normalization should therefore be performed before filtering miRNAs.
Applying the function to an already filtered matrix changes the
denominator and the biological meaning of the resulting RPM values.

The deprecated \`metrics\` and \`reads_column\` arguments are retained
temporarily for compatibility with miRPM 0.1.0 pipelines.

The returned matrix stores normalization information in the
\`"miRPM_normalization"\` attribute. This includes the library sizes,
their source and the approximate number of reads represented by one RPM
in each sample.

## Examples

``` r
count_matrix <- matrix(
  c(100, 200, 300, 400),
  nrow = 2,
  dimnames = list(
    c("miR-1", "miR-2"),
    c("Sample1", "Sample2")
  )
)

rpm_matrix <- normalize_rpm(count_matrix)

colSums(rpm_matrix)
#> Sample1 Sample2 
#>   1e+06   1e+06 

attr(rpm_matrix, "miRPM_normalization")
#> $method
#> [1] "RPM"
#> 
#> $scale
#> [1] 1e+06
#> 
#> $library_sizes
#> Sample1 Sample2 
#>     300     700 
#> 
#> $reads_per_rpm
#> Sample1 Sample2 
#>   3e-04   7e-04 
#> 
#> $library_size_source
#> [1] "matrix_column_sums"
#> 
#> $input_features
#> [1] 2
#> 
#> $input_samples
#> [1] 2
#> 
```
