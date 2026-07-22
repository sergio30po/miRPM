# Generate dot plots for differentially expressed miRNAs

This function generates dot plots for differentially expressed miRNAs
between specific groups. The plots are saved as interactive HTML files
in a folder called "Interactive_plots" with a specified output name. A
customizable title can be added to the plot.

## Usage

``` r
miRNA_expression_plot(
  miRNAs_DE,
  miRNA_ftd,
  metadata,
  condition_column,
  groups,
  colors,
  output_name,
  plot_title = "miRNA Expression Plot"
)
```

## Arguments

- miRNAs_DE:

  File containing differentially expressed miRNAs (miRNAs should be in
  the first column).

- miRNA_ftd:

  Count matrix with miRNAs as row names and samples as column names.

- metadata:

  Data frame containing sample information. Must have a column with
  sample names and a condition column.

- condition_column:

  Name of the column in \`metadata\` that contains the group
  classifications.

- groups:

  Vector specifying the two groups to compare within the
  \`condition_column\`.

- colors:

  Named vector specifying the colors for each group.

- output_name:

  Name for the output HTML file.

- plot_title:

  Title for the plot (default: "miRNA Expression Plot").

## Value

A list containing:

- \`ggplot\`: A ggplot object of the dot plot.

- \`interactive\`: An interactive plotly object of the dot plot.

- \`output_folder\`: The path to the "Interactive_plots" folder where
  the HTML file is saved.

## Examples

``` r
if (FALSE) { # \dontrun{
output <- miRNA_expression_plot(
  miRNAs_DE = "differential_miRNAs.xlsx",
  miRNA_ftd = miRNA_ftd,
  metadata = metadata,
  condition_column = "Condition",
  groups = c("Control", "Disease"),
  colors = c("Control" = "blue", "Disease" = "red"),
  output_name = "miRNA_expression",
  plot_title = "Differentially Expressed miRNAs in Control vs Disease"
)

output$ggplot
output$interactive
output$output_folder
} # }
```
