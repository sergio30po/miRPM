# Filter miRNAs based on expression in a minimum percentage of samples in any group

This function filters miRNAs that have reads in at least a specified
percentage of samples in any group. The comparison operator for both the
threshold and the read count can be customized (\>, \>=, \<, \<=).

## Usage

``` r
filter_mirnas(
  count_matrix,
  metadata,
  threshold = 0.5,
  min_reads = 1000,
  threshold_comparison = ">=",
  read_comparison = ">"
)
```

## Arguments

- count_matrix:

  A count matrix (miRNAs in rows, samples in columns).

- metadata:

  A data frame containing sample information, including a 'Condition'
  column with the group for each sample.

- threshold:

  Minimum percentage of samples with reads (value between 0 and 1,
  default is 0.5).

- min_reads:

  Minimum number of reads to consider a miRNA expressed in a sample
  (default is 1000).

- threshold_comparison:

  Comparison operator for the threshold. Options are "\>", "\>=", "\<",
  or "\<=" (default is "\>=").

- read_comparison:

  Comparison operator for the read count. Options are "\>", "\>=", "\<",
  or "\<=" (default is "\>").

## Value

A filtered count matrix containing miRNAs that meet the criteria in at
least one group.

## Examples

``` r
# Example usage:
# count_matrix <- matrix(c(1000, 0, 2000, 1500, 0, 3000), nrow = 2, dimnames = list(c("miRNA1", "miRNA2"), c("Sample1", "Sample2", "Sample3")))
# metadata <- data.frame(Sample = c("Sample1", "Sample2", "Sample3"), Condition = c("A", "A", "B"))
# filtered_matrix <- filter_mirnas(count_matrix, metadata, threshold = 0.5, min_reads = 1000, threshold_comparison = ">=", read_comparison = ">")
```
