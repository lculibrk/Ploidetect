% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{curate_segments}
\alias{curate_segments}
\title{Merge noisy segments with their locally dominant segment}
\usage{
curate_segments(cna, regions)
}
\arguments{
\item{cna}{a data.table containing the bin-level copy number calls}

\item{regions}{the output from map_noisy_segments}
}
\value{
a \code{cna} data.table with the modified segments.
}
\description{
\code{curate_segments} takes in \code{cna} of the bin-level copy number calls and
\code{regions} from \code{map_noisy_segments}. Multiple outputs from \code{map_noisy_segments}
can be combined with \code{c()}
}
