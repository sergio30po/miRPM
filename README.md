# miRPM

<!-- badges: start -->
[![R-CMD-check](https://github.com/sergio30po/miRPM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sergio30po/miRPM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/sergio30po/miRPM/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/sergio30po/miRPM/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- badges: end -->

`miRPM` is an R package for transparent and reproducible analysis of
microRNA sequencing count data using reads per million (RPM)
normalization.

The package provides a modular workflow for:

- RPM normalization with explicit library-size handling.
- Abundance and prevalence filtering.
- Non-parametric analysis of group differences.
- Multiple-testing correction and effect-size estimation.
- Static and interactive expression plots.
- Optional Excel, HTML and PNG output.

## Scope

RPM normalization rescales each sample by its sequencing depth:

$$
\mathrm{RPM}_{ij} =
\frac{\mathrm{count}_{ij}}
{\mathrm{library\ size}_{j}}
\times 10^6
$$

where $i$ identifies a miRNA and $j$ identifies a sample.

RPM is a transparent scaling method, not a count-based statistical
model. It is useful for exploratory analysis, visualization and
workflows in which normalized abundance has a direct interpretation.

For confirmatory differential-expression analyses, results should be
benchmarked against established count-based methods such as DESeq2,
edgeR or limma-voom, particularly when composition bias, complex
designs or strong mean-variance relationships are expected.

## Installation

Install the development version from GitHub with `pak`:

```r
# install.packages("pak")
pak::pak("sergio30po/miRPM")
```

Alternatively, use `remotes`:

```r
# install.packages("remotes")
remotes::install_github("sergio30po/miRPM")
```

Load the package:

```r
library(miRPM)
```

## Core functions

| Function | Purpose |
|---|---|
| `normalize_rpm()` | Normalize a complete count matrix using column sums or externally supplied library sizes. |
| `filter_mirnas()` | Filter miRNAs by RPM abundance, raw-count support, prevalence and group-specific criteria. |
| `perform_statistical_tests()` | Run Mann-Whitney-Wilcoxon or Kruskal-Wallis/Dunn analyses with FDR correction and effect sizes. |
| `miRNA_expression_plot()` | Create a joint static and interactive expression plot for selected miRNAs. |
| `individual_miRNA_plot()` | Create one static expression plot per selected miRNA, with optional PNG export. |

## Minimal workflow

### Example data

`miRPM` includes a small synthetic dataset for reproducible examples. It
contains six miRNAs measured in eight samples divided into two groups.

```r
data("miRPM_example", package = "miRPM")

counts <- miRPM_example$counts
metadata <- miRPM_example$metadata

dim(counts)
metadata

stopifnot(
  identical(colnames(counts), metadata$Sample)
)
```

The count matrix must contain miRNAs in rows and samples in columns.
Sample identifiers in `colnames(counts)` must correspond to the values
in `metadata$Sample`.

### 1. Normalize to RPM

By default, `normalize_rpm()` uses the complete column sums of the input
matrix as library sizes.

```r
rpm <- normalize_rpm(counts)
rpm
```

Externally calculated library sizes can also be supplied as a named
vector:

```r
library_sizes <- setNames(
  c(
    2.1e6, 2.3e6, 2.0e6, 2.2e6,
    1.9e6, 2.4e6, 2.1e6, 2.0e6
  ),
  colnames(counts)
)

rpm_external <- normalize_rpm(
  count_matrix = counts,
  library_sizes = library_sizes
)
```

Normalization metadata are stored with the returned matrix:

```r
attr(rpm, "miRPM_normalization")
```

### 2. Filter miRNAs

The example below retains miRNAs reaching at least 5 RPM in at least
50% of samples from any group. A minimum raw-count criterion is also
applied.

```r
filtered <- filter_mirnas(
  rpm_matrix = rpm,
  metadata = metadata,
  raw_counts = counts,
  min_rpm = 5,
  min_count = 10,
  min_prevalence = 0.5,
  min_samples = 2,
  group_column = "Condition",
  prevalence_scope = "any_group",
  return_diagnostics = TRUE
)

filtered$filtered_matrix
filtered$summary
filtered$diagnostics
```

Filtering thresholds are study-dependent and should be selected before
testing, justified biologically and evaluated through sensitivity
analyses.

### 3. Analyze group differences

For two independent groups, `perform_statistical_tests()` uses the
Mann-Whitney-Wilcoxon rank-sum test.

```r
de_result <- perform_statistical_tests(
  miRNA_ftd = filtered$filtered_matrix,
  metadata = metadata,
  condition = "Condition",
  group_order = c("Control", "Disease"),
  output_file = NULL,
  assign_results = FALSE
)

de_result$results
de_result$significant_results
```

The returned results include:

- Raw and adjusted p-values.
- Group-specific descriptive statistics.
- Detection prevalence.
- Rank-biserial correlation.
- Probability of superiority.
- Mean and median differences.
- A descriptive log2 median ratio.

For three or more groups, the same function uses Kruskal-Wallis tests
and can generate Dunn post-hoc comparisons.

### 4. Create a joint interactive plot

```r
joint_plot <- miRNA_expression_plot(
  mirnas = rownames(filtered$filtered_matrix),
  rpm_matrix = filtered$filtered_matrix,
  metadata = metadata,
  condition_column = "Condition",
  sample_column = "Sample",
  groups = c("Control", "Disease"),
  save_html = FALSE
)

joint_plot$ggplot
joint_plot$interactive
```

To save the interactive plot:

```r
joint_plot_saved <- miRNA_expression_plot(
  mirnas = rownames(filtered$filtered_matrix),
  rpm_matrix = filtered$filtered_matrix,
  metadata = metadata,
  condition_column = "Condition",
  groups = c("Control", "Disease"),
  output_name = "miRNA-expression",
  output_dir = "Interactive_plots",
  save_html = TRUE
)
```

### 5. Create individual miRNA plots

```r
individual_plots <- individual_miRNA_plot(
  filtered_results = de_result$results,
  rpm_matrix = filtered$filtered_matrix,
  metadata = metadata,
  condition_column = "Condition",
  sample_column = "Sample",
  groups_to_include = c("Control", "Disease"),
  adjusted_pvalue_column = "FDR",
  save_png = FALSE
)

individual_plots$plots
```

Individual PNG files can be generated with `save_png = TRUE`.

## Statistical interpretation

A small adjusted p-value alone should not define biological relevance.
Interpretation should also consider:

- Effect-size magnitude and direction.
- Expression abundance.
- Detection prevalence.
- Consistency across samples.
- Group imbalance and zero inflation.
- Sensitivity to filtering and normalization choices.
- Replication in independent datasets.

The Wilcoxon rank-sum test evaluates differences in rank distributions.
It should not be described simply as a test of medians unless additional
distributional assumptions are satisfied.

## Reproducibility

The package currently includes:

- Automated tests for all exported functions.
- Local `R CMD check` validation.
- Multi-platform GitHub Actions checks on Linux, Windows and macOS.
- A pkgdown workflow for automated website deployment.

Record the package version and R session used for each analysis:

```r
packageVersion("miRPM")
sessionInfo()
```

## Documentation

Full function documentation and examples are available at:

<https://sergio30po.github.io/miRPM/>

## Development status

The current development version is preparing `miRPM` 0.2.0. The API
retains compatibility with the original 0.1.0 interfaces where
practical, while new code should use the named arguments documented in
the current function references.

## License

`miRPM` is released under the MIT License.
