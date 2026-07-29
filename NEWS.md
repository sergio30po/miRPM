# miRPM 0.2.0

## Package structure and reproducibility

- Reorganized `miRPM` as a standard R package with documented public
  functions, automated tests, continuous integration, and a pkgdown website.
- Added the synthetic `miRPM_example` dataset for reproducible examples,
  package testing, and documentation.
- Added multi-platform R CMD check workflows for released, development, and
  previous R versions on Linux, macOS, and Windows.

## RPM normalization

- Updated `normalize_rpm()` to normalize complete miRNA count matrices using
  sample-specific library sizes.
- By default, library sizes are calculated from the complete column sums of
  the input count matrix.
- Added support for externally calculated, named library sizes.
- Added validation of count matrices, sample names, library sizes, and
  normalization inputs.
- Added normalization metadata to the returned RPM matrix.

## miRNA filtering

- Redesigned `filter_mirnas()` to combine RPM abundance, raw-count,
  prevalence, and minimum-sample criteria.
- Added filtering across all samples, any group, all groups, or a specified
  minimum number of groups.
- Added detailed filtering diagnostics and summary information.
- Retained selected compatibility arguments from the original interface.

## Statistical analysis

- Redesigned `perform_statistical_tests()` for nonparametric group
  comparisons of filtered RPM data.
- Added Wilcoxon rank-sum testing for two independent groups.
- Added Kruskal-Wallis testing and Dunn post-hoc comparisons for analyses
  involving more than two groups.
- Added multiple-testing correction, including Benjamini-Hochberg false
  discovery rate adjustment.
- Added descriptive statistics, effect-size information, and structured
  result tables.
- Added optional Excel export while allowing analyses to run without writing
  files.

## Visualization

- Redesigned `miRNA_expression_plot()` to generate both static `ggplot2`
  graphics and interactive `plotly` visualizations.
- Redesigned `individual_miRNA_plot()` to generate one expression plot per
  selected miRNA.
- Added sample alignment, group selection, color validation, tooltip
  information, and optional file export.
- Functions now return structured objects containing plots, data, selected
  miRNAs, missing miRNAs, groups, and output paths.

## Documentation

- Added a pkgdown website with Home, Reference, and License sections.
- Expanded the README with a complete normalization, filtering, statistical
  analysis, and visualization workflow.
- Documented the orientation and alignment requirements for count matrices
  and sample metadata.
- Clarified that RPM is a transparent abundance-scaling method and is not a
  replacement for count-based differential-expression models.
