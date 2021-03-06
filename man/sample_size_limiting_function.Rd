% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_size_limiting_function.R
\name{sample_size_limiting_function}
\alias{sample_size_limiting_function}
\title{Title: Limit the sample sizes required to only the MDC values}
\usage{
sample_size_limiting_function(stages_info, MDC_main_table, sample_size_table)
}
\arguments{
\item{stages_info}{vector of strings denoting which levels ("BL", "LH", "MH", and "HH")}

\item{MDC_main_table}{dataframe of the MDC estimates as standardised effect sizes}

\item{sample_size_table}{dataframe of the large sample size dataframe}
}
\value{
dataframe of sample sizes required for each MDC estimate
}
\description{
NOTES: Used for sample size tab in the supplementary material document
}
\keyword{coefficient}
\keyword{coefficient,}
\keyword{correlation}
\keyword{intra-class}
\keyword{of}
\keyword{reliability,}
\keyword{variation}
