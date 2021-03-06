% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{view_segments}
\alias{view_segments}
\title{Produce ploidetect-style plots for a copy number profile at the given chromosome and positions}
\usage{
view_segments(cna, chrom, plot_pos, plot_end)
}
\arguments{
\item{cna}{a data.table containing the bin-level copy number calls}

\item{chrom}{the chromosome of interest}

\item{plot_pos}{the start position of the plot}

\item{plot_end}{the end position of the plot}
}
\value{
a ggplot2 plot
}
\description{
\code{view_segments} is a wrapper for plot_ploidetect to make filtering for specific regions simple
}
