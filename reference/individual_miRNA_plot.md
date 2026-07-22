# Individual miRNA plot for differentially expressed miRNAs

This function generates dot plots for differentially expressed miRNAs
between specific groups. The plots are saved as PNG files in a folder
named after the comparison, inside a main folder called
"Individual_plots".

## Usage

``` r
individual_miRNA_plot(
  filtered_results,
  rpm_matrix,
  metadata,
  condition_column,
  sample_column,
  groups_to_include,
  condition_colors,
  adjusted_pvalue_column,
  output_dir = "Individual_plots"
)
```

## Arguments

- filtered_results:

  Data frame with the filtered results (FDR \< 0.05). The first column
  should contain miRNA names.

- rpm_matrix:

  Matrix of RPMs (miRNAs in rows, samples in columns).

- metadata:

  Data frame with sample group information. Must have a column for
  sample names and another for groups.

- condition_column:

  Name of the column in \`metadata\` that contains the groups.

- sample_column:

  Name of the column in \`metadata\` that contains the sample names.

- groups_to_include:

  Vector with the names of the groups to include in the comparison.

- condition_colors:

  Named vector of colors for each condition. The names should match the
  groups in \`groups_to_include\`.

- adjusted_pvalue_column:

  Name of the column in \`filtered_results\` that contains the adjusted
  p-value.

- output_dir:

  Name of the main folder where the plots will be saved (default:
  "Individual_plots").

## Value

No return value. The plots are saved as PNG files in the specified
folder structure.

## Examples

``` r
if (FALSE) { # \dontrun{
individual_miRNA_plot(
  filtered_results = filtered_results,
  rpm_matrix = rpm_matrix,
  metadata = metadata,
  condition_column = "Condition",
  sample_column = "Samples",
  groups_to_include = c("Group1", "Group2"),
  condition_colors = c("Group1" = "blue", "Group2" = "red"),
  adjusted_pvalue_column = "adj.P.Val"
)
} # }
```
