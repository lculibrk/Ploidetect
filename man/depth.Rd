% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc_functions.R
\name{depth}
\alias{depth}
\title{Compute predicted depths based on coverage characteristics}
\usage{
depth(maxpeak, d, P, n)
}
\arguments{
\item{maxpeak}{the depth of the most common copy number of the genome 
(eg. where CN = Ploidy)}

\item{d}{the depth difference per copy number. Can be calculated from 
\code{get_coverage_characteristics}}

\item{P}{the ploidy (most common copy number)}

\item{n}{a vector of all copy numbers to compute depths for}
}
\value{
a numeric vector of predicted depths for each n, with names equal to 
n
}
\description{
\code{depth} returns the depth of of each specified copy number given the 
input parameters
}
