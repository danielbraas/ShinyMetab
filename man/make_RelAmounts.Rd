% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_data.R
\name{make_RelAmounts}
\alias{make_RelAmounts}
\title{This function produces a data frame with relative metabolite amounts.}
\usage{
make_RelAmounts(DF, Title = "", folder = getwd())
}
\arguments{
\item{DF}{The input data. This should be a data frame with Name, Iso,
Condition, Exp, Value, Used_ID, KEGG.ID, Nr.C, Rt and Formula columns.}

\item{Title}{The title to be used.}

\item{folder}{The folder for the data to be saved. Default is the current directory.}
}
\value{
A data frame with relative amount information both averaged as well as averaged and normalized to the first condition.
}
\description{
This function produces a data frame with relative metabolite amounts.
}
\author{
Daniel Braas
}
