% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_data.R
\name{make_FC}
\alias{make_FC}
\title{This function produces a data frame with fractional contribution data.}
\usage{
make_FC(DF, Title = "", folder = getwd())
}
\arguments{
\item{DF}{The input data. This should be a data frame that contains MID data generated
with the 'make_MID' function.}

\item{Title}{The title to be used.}

\item{folder}{The folder where the data is supposed to be saved. Default is the current directory.}
}
\value{
A data frame with fractional contribution data.
}
\description{
This function produces a data frame with fractional contribution data.
}
\author{
Daniel Braas
}
