# Filter miRNAs by abundance and prevalence

Filters a normalized RPM matrix using configurable abundance,
prevalence, sample-count and group-specific criteria.

## Usage

``` r
filter_mirnas(
  rpm_matrix = NULL,
  metadata = NULL,
  threshold = NULL,
  min_reads = NULL,
  threshold_comparison = NULL,
  read_comparison = NULL,
  count_matrix = NULL,
  raw_counts = NULL,
  min_rpm = 5,
  min_count = NULL,
  min_prevalence = 0.5,
  min_samples = NULL,
  group_column = NULL,
  prevalence_scope = c("overall", "any_group", "all_groups"),
  min_groups = 1,
  return_diagnostics = FALSE
)
```

## Arguments

- rpm_matrix:

  A numeric matrix or data frame containing RPM values, with miRNAs in
  rows and samples in columns.

- metadata:

  Optional data frame containing sample information. Grouped filtering
  requires a \`Sample\` column and the column named by \`group_column\`.

- threshold:

  Deprecated. Minimum prevalence used by the miRPM 0.1.0 interface.

- min_reads:

  Deprecated. Minimum expression value used by the miRPM 0.1.0
  interface. In the historical workflow this was normally applied to an
  RPM-normalized matrix despite the argument name.

- threshold_comparison:

  Deprecated comparison operator for \`threshold\`.

- read_comparison:

  Deprecated comparison operator for \`min_reads\`.

- count_matrix:

  Deprecated named alias for the normalized matrix used by miRPM 0.1.0.

- raw_counts:

  Optional raw count matrix corresponding to \`rpm_matrix\`. When
  supplied with \`min_count\`, a sample must satisfy both the RPM and
  raw-count thresholds.

- min_rpm:

  Minimum RPM required for a miRNA to be considered detected in a
  sample. The default is 5.

- min_count:

  Optional minimum raw count required in addition to \`min_rpm\`. When
  \`raw_counts\` is omitted, counts are inferred from the normalization
  metadata created by \`normalize_rpm()\`.

- min_prevalence:

  Minimum proportion of samples in which a miRNA must satisfy the
  abundance criterion. Must be between 0 and 1.

- min_samples:

  Optional minimum absolute number of samples that must satisfy the
  abundance criterion. When both \`min_prevalence\` and \`min_samples\`
  are supplied, both requirements are enforced.

- group_column:

  Optional name of the metadata column defining groups.

- prevalence_scope:

  Filtering scope. One of \`"overall"\`, \`"any_group"\` or
  \`"all_groups"\`.

- min_groups:

  Minimum number of groups in which the criterion must be satisfied when
  \`prevalence_scope = "any_group"\`.

- return_diagnostics:

  If \`TRUE\`, return a list containing the filtered matrix,
  diagnostics, summary and parameters. If \`FALSE\`, return only the
  filtered matrix and store the same information as attributes.

## Value

A filtered RPM matrix, or a list containing the filtered matrix,
diagnostics, summary and parameters when \`return_diagnostics = TRUE\`.

## Details

The recommended workflow is to normalize the complete count matrix with
\`normalize_rpm()\` before filtering. Filtering thresholds are
configurable because their interpretation depends on sequencing depth,
sample type, library complexity and study objectives.

A provisional starting profile for differential-expression analysis is
\`min_rpm = 5\`, \`min_count = 15\`, \`min_prevalence = 0.5\`,
\`min_samples = 3\` and \`prevalence_scope = "any_group"\`. These values
are recommendations rather than universal requirements.

With \`prevalence_scope = "overall"\`, the criterion is evaluated across
all samples. With \`"any_group"\`, it must be satisfied in at least
\`min_groups\` groups. With \`"all_groups"\`, it must be satisfied in
every group.

The historical miRPM 0.1.0 interface is retained temporarily for
backward compatibility.

## Examples

``` r
rpm_matrix <- matrix(
  c(
    5, 6, 0, 0,
    0, 0, 10, 0,
    0, 0, 0, 0
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    c("miR-1", "miR-2", "miR-3"),
    c("S1", "S2", "S3", "S4")
  )
)

metadata <- data.frame(
  Sample = c("S1", "S2", "S3", "S4"),
  Group = c("Control", "Control", "Disease", "Disease")
)

filtered <- filter_mirnas(
  rpm_matrix = rpm_matrix,
  metadata = metadata,
  min_rpm = 5,
  min_prevalence = 0.5,
  group_column = "Group",
  prevalence_scope = "any_group"
)
```
