# Perform Statistical Tests on miRNA Expression Data

This function performs statistical tests on miRNA expression data based
on the number of groups in the specified condition. For two groups, it
performs the Mann-Whitney U test. For three or more groups, it performs
the Kruskal-Wallis test followed by Dunn's post hoc test. Results are
saved in an Excel file inside the "Tests_results" folder with separate
sheets for each comparison. Each sheet will also include a column
indicating the significance of the results using symbols: \`\*\` for
p-value \< 0.05, \`\*\*\` for p-value \< 0.01, and \`\*\*\*\` for
p-value \< 0.005.

## Usage

``` r
perform_statistical_tests(miRNA_ftd, metadata, condition, output_file)
```

## Arguments

- miRNA_ftd:

  A count matrix of miRNA expression data (miRNAs in rows, samples in
  columns).

- metadata:

  A dataframe containing sample information, including the condition to
  test.

- condition:

  The name of the column in \`metadata\` that defines the groups to
  compare.

- output_file:

  The name of the output Excel file (e.g., "Results.xlsx").

## Value

An Excel file with statistical results. For Kruskal-Wallis, it includes
a sheet for the global test and sheets for Dunn's post hoc comparisons.
For Mann-Whitney, it includes a single sheet with the test results. Each
sheet will have a column titled "Significance" showing the corresponding
symbol (\`\*\`, \`\*\*\`, or \`\*\*\*\`).

A dataframe for each comparison.

The path to the output folder ("Tests_results").

## Examples

``` r
# Example usage:
# miRNA_ftd <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, dimnames = list(c("miRNA1", "miRNA2"), c("Sample1", "Sample2", "Sample3")))
# metadata <- data.frame(Sample = c("Sample1", "Sample2", "Sample3"), Condition = c("A", "A", "B"))
# perform_statistical_tests(miRNA_ftd, metadata, "Condition", "Test_Results.xlsx")
```
