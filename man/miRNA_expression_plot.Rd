% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dotplot.R
\name{miRNA_expression_plot}
\alias{miRNA_expression_plot}
\title{Generate dot plots for differentially expressed miRNAs}
\usage{
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
}
\arguments{
\item{miRNAs_DE}{File containing differentially expressed miRNAs (miRNAs should be in the first column).}

\item{miRNA_ftd}{Count matrix with miRNAs as row names and samples as column names.}

\item{metadata}{Data frame containing sample information. Must have a column with sample names and a condition column.}

\item{condition_column}{Name of the column in `metadata` that contains the group classifications.}

\item{groups}{Vector specifying the two groups to compare within the `condition_column`.}

\item{colors}{Named vector specifying the colors for each group.}

\item{output_name}{Name for the output HTML file.}

\item{plot_title}{Title for the plot (default: "miRNA Expression Plot").}
}
\value{
A list containing:
\itemize{
  \item `ggplot`: A ggplot object of the dot plot.
  \item `interactive`: An interactive plotly object of the dot plot.
  \item `output_folder`: The path to the "Interactive_plots" folder where the HTML file is saved.
}
}
\description{
This function generates dot plots for differentially expressed miRNAs between specific groups.
The plots are saved as interactive HTML files in a folder called "Interactive_plots" with a specified output name.
A customizable title can be added to the plot.
}
\examples{
\dontrun{
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
}
}
