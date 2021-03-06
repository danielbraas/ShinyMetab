% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_data.R
\name{medium_RelA}
\alias{medium_RelA}
\title{medium_RelA calculates statistical metrics from all replicates (including fresh medium blanks).}
\usage{
medium_RelA(dat, Title = "")
}
\arguments{
\item{dat}{A formated data frame produced by prepare_data().}

\item{Title}{Project title.}
}
\value{
A data frame with integrated metabolomics data as well as
bar graphs of all measured metabolites.
}
\description{
medium_RelA calculates statistical metrics from all replicates (including fresh medium blanks).
}
\author{
Daniel Braas
}
