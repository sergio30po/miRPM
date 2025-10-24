miRPM: miRNA RPM Normalization
================

# miRPM

------------------------------------------------------------------------

## Overview

**miRPM** is a lightweight R package for **normalization and analysis of
microRNA-seq (miRNA-seq) count data** using *Reads Per Million (RPM)*.  
This approach provides a straightforward and biologically interpretable
way to compare miRNA expression levels across samples with different
sequencing depths.

miRPM is designed for reproducible analysis of small RNA-seq datasets,
providing core functionalities for **RPM normalization**, **filtering**,
**visualization**, and **basic statistical testing**.

------------------------------------------------------------------------

## Concept

Normalization is a critical step in miRNA-seq analysis to correct for
differences in library size and sequencing depth between samples.  
While advanced methods (e.g., DESeq2, TMM, RUV) handle composition
biases, **RPM normalization** remains an essential and transparent
approach for exploratory analyses and publication-ready summaries.

RPM normalization rescales each sample by its total read count,
expressed in millions, using:

$$
RPM_{ij} = \frac{counts_{ij}}{\sum_i counts_{ij}} \times 10^6
$$

where $counts_{ij}$ is the raw count for miRNA *i* in sample *j*.

------------------------------------------------------------------------

## Functionality

miRPM provides simple, modular functions for different stages of the
normalization and exploration workflow:

| Function | Description |
|----|----|
| `normalize_rpm()` | Performs RPM normalization on a count matrix. |
| `filter_mirnas()` | Removes low-abundance miRNAs based on count thresholds. |
| `miRNA_expression_plot()` | Visualizes normalized miRNA expression distributions. |
| `individual_miRNA_plot()` | Displays the expression pattern of a single miRNA across conditions. |
| `perform_statistical_tests()` | Conducts non-parametric statistical tests (e.g., Kruskal–Wallis, Dunn). |

------------------------------------------------------------------------

## Metrics and Interpretation

RPM normalization is not a statistical model—it is a **scaling
method**.  
Its value lies in interpretability and consistency, especially for small
to medium sample sizes or visualization purposes.

In practice, you can evaluate the normalization by: - Inspecting overall
expression distributions before and after normalization.  
- Assessing whether sample-specific sequencing depth bias has been
removed.  
- Comparing normalized values across experimental groups.

------------------------------------------------------------------------

## Installation and Usage

You can install the development version directly from GitHub using
**remotes**:

``` r
# install.packages("remotes")
remotes::install_github("sergio30po/miRPM")
```

Please refer to the [documentation page](https://sergio30po.github.io/miRPM)
for detailed function references.
