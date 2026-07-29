# Synthetic miRNA-Seq Example Data

A small synthetic dataset designed to demonstrate the normalization,
filtering, statistical-analysis, and plotting functions provided by
\`miRPM\`.

## Usage

``` r
miRPM_example
```

## Format

A list with two components:

- counts:

  An integer matrix containing six miRNAs in rows and eight samples in
  columns.

- metadata:

  A data frame containing one row per sample and the variables
  \`Sample\` and \`Condition\`.

## Source

Synthetic data generated for the \`miRPM\` package.

## Details

The dataset is intended exclusively for package examples and testing. It
does not represent measurements from biological samples and must not be
used for scientific benchmarking or biological interpretation.
