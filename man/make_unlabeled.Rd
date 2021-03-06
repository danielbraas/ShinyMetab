% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_data.R
\name{make_unlabeled}
\alias{make_unlabeled}
\title{This function produces a data frame with fractional contribution data.}
\usage{
make_unlabeled(DF, Title = "", max.iso = 50, folder = getwd())
}
\arguments{
\item{DF}{The input data. This should be a data frame that contains the MID envelope.
The data frame should have the following columns: Used_ID, Iso and each sample separated
into individual columns.}

\item{Title}{The title to be used.}

\item{max.iso}{The maximal number of isotopologues in the data (default is 50).}

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
