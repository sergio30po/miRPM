# Perform non-parametric differential-expression analysis

Performs rank-based differential-expression analysis on a normalized
miRNA expression matrix. Two-group comparisons use the
Mann-Whitney-Wilcoxon rank-sum test. Comparisons involving three or more
groups use the Kruskal-Wallis test followed by Dunn post-hoc
comparisons.

## Usage

``` r
perform_statistical_tests(
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
)
```

## Arguments

- miRNA_ftd:

  A numeric matrix or data frame containing normalized expression
  values, with miRNAs in rows and samples in columns.

- metadata:

  A data frame containing sample information. It must include a
  \`Sample\` column and the column named by \`condition\`.

- condition:

  A single character string naming the metadata column that defines the
  comparison groups.

- output_file:

  Optional name of the Excel workbook. When \`NULL\`, no workbook is
  written.

- assign_results:

  Logical. If \`TRUE\`, significant legacy result objects are also
  assigned to the calling environment for backward compatibility.
  Default is \`TRUE\`.

- output_dir:

  Directory used when an Excel workbook is written. Default is
  \`"Tests_results"\`.

- alpha:

  FDR threshold used to define significant results. Default is \`0.05\`.

- p_adjust_method:

  Multiple-testing correction passed to \`stats::p.adjust()\`. Default
  is \`"BH"\`.

- detection_threshold:

  Expression value above which a miRNA is considered detected when
  prevalence is calculated. Default is \`0\`.

- pseudocount:

  Non-negative value added to group medians before calculating the
  descriptive log2 median ratio. Default is \`1\`.

- group_order:

  Optional character vector defining group order. For a two-group
  comparison, positive rank-biserial correlation indicates greater
  values in the first group than in the second group.

- posthoc_only_if_global_significant:

  Logical. For analyses with three or more groups, perform Dunn
  comparisons only for miRNAs with a Kruskal-Wallis FDR below \`alpha\`.
  Default is \`TRUE\`.

- write_excel:

  Logical. If \`TRUE\`, write an Excel workbook. By default, this is
  \`TRUE\` when \`output_file\` is supplied and \`FALSE\` otherwise.

## Value

A structured list containing:

- \`analysis_type\`: \`"two_group"\` or \`"multi_group"\`.

- \`groups\`: ordered group names.

- \`group_sizes\`: sample counts by group.

- \`results\`: complete Mann-Whitney or Kruskal-Wallis results.

- \`significant_results\`: rows with FDR below \`alpha\`.

- \`posthoc_results\`: complete Dunn results for multigroup analyses.

- \`significant_posthoc_results\`: significant Dunn results.

- \`comparison_results\`: legacy-compatible significant result tables.

- \`parameters\`: analysis settings.

- \`output_file\`: written workbook path, or \`NULL\`.

- \`output_folder\`: workbook directory, or \`NULL\`.

## Details

The function expects an already normalized and, when appropriate,
filtered expression matrix. Filtering and normalization should be
completed before differential-expression testing.

For two groups, the function reports the Wilcoxon statistic, raw
p-value, FDR, rank-biserial correlation, probability of superiority,
differences between group means and medians, a descriptive log2 median
ratio, and group-specific abundance and prevalence summaries.

For three or more groups, the function reports the Kruskal-Wallis
statistic, FDR and epsilon-squared effect size. Dunn post-hoc p-values
are adjusted across miRNAs separately for each pairwise comparison.

The Wilcoxon rank-sum test should not automatically be interpreted as a
test of medians. Rank-biserial correlation and probability of
superiority describe the direction and magnitude of rank separation
between groups.

The log2 median ratio is descriptive and is not used as a mandatory
significance threshold.

## Examples

``` r
expression_matrix <- matrix(
  c(
    1, 2, 1, 8, 9, 10,
    5, 4, 6, 5, 4, 6
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("miR-1", "miR-2"),
    paste0("S", 1:6)
  )
)

metadata <- data.frame(
  Sample = paste0("S", 1:6),
  Condition = rep(c("Control", "Disease"), each = 3)
)

result <- perform_statistical_tests(
  miRNA_ftd = expression_matrix,
  metadata = metadata,
  condition = "Condition",
  output_file = NULL,
  assign_results = FALSE
)

result$results
#>   miRNA group_1 group_2 wilcox_statistic mann_whitney_u_group_1    p_value
#> 1 miR-1 Control Disease              0.0                    0.0 0.04630159
#> 2 miR-2 Control Disease              4.5                    4.5 1.00000000
#>   rank_biserial_group_1_vs_group_2 probability_group_1_greater_group_2
#> 1                               -1                                 0.0
#> 2                                0                                 0.5
#>   median_difference_group_1_minus_group_2 mean_difference_group_1_minus_group_2
#> 1                                      -8                             -7.666667
#> 2                                       0                              0.000000
#>   log2_median_ratio_group_1_vs_group_2        FDR p_significance
#> 1                            -2.321928 0.09260319              *
#> 2                             0.000000 1.00000000               
#>   FDR_significance overall_mean overall_median overall_iqr overall_mad
#> 1                      5.166667              5         7.5      5.9304
#> 2                      5.000000              5         1.5      1.4826
#>   overall_minimum overall_maximum overall_detected_samples overall_prevalence
#> 1               1              10                        6                  1
#> 2               4               6                        6                  1
#>   overall_zero_fraction n_Control mean_Control median_Control iqr_Control
#> 1                     0         3     1.333333              1         0.5
#> 2                     0         3     5.000000              5         1.0
#>   mad_Control minimum_Control maximum_Control detected_samples_Control
#> 1      0.0000               1               2                        3
#> 2      1.4826               4               6                        3
#>   prevalence_Control zero_fraction_Control n_Disease mean_Disease
#> 1                  1                     0         3            9
#> 2                  1                     0         3            5
#>   median_Disease iqr_Disease mad_Disease minimum_Disease maximum_Disease
#> 1              9           1      1.4826               8              10
#> 2              5           1      1.4826               4               6
#>   detected_samples_Disease prevalence_Disease zero_fraction_Disease
#> 1                        3                  1                     0
#> 2                        3                  1                     0
```
