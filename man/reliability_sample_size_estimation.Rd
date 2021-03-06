% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reliability_sample_size_estimation.R
\name{reliability_sample_size_estimation}
\alias{reliability_sample_size_estimation}
\title{Title: Sample Size Estimations For Reliability Studies}
\usage{
reliability_sample_size_estimation(w, n, z, p_h, type)
}
\arguments{
\item{w}{numeric value or vector of values for width of the confidence interval}

\item{n}{numeric value for number of replicates}

\item{z}{numeric value indicating the z score value for a given confidence interval (e.g., 1.96 for 95\% CI)}

\item{p_h}{numeric value or vector of values for expected/planned ICC}

\item{type}{string input taking form of "single" or "multiple".}
}
\value{
dataframe consisting of two values (type = "single") for k and K1, or large dataframe with k and k1 values for a range of w and p_h values.
}
\description{
Function that takes a given width of a confidence interval (w), number of replicates (n), z score for confidence interval required (z), and planned ICC value (p_h)
Formulas given by Shoukri et al., 2004 (equation 7)
Shoukri MM, Asyali MH & Donner A (2004). Sample size requirements for the design of reliability study: review and new results. Statistical Methods in Medical Research 13, 251-271.
}
\details{
For: w = 0.2, n = 3, z = 1.96, p_h = 0.8, type = "single, should equal k = 35, k1 = 36 to fit with Table 3 in Shoukri et al., 2004
}
\keyword{coefficient}
\keyword{correlation}
\keyword{intra-class}
\keyword{reliability,}
