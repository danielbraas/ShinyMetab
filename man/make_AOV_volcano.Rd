% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_AOV_volcano.R
\name{make_AOV_volcano}
\alias{make_AOV_volcano}
\title{This function takes in a data and produces a volcano plot relative to a specified condition level.}
\usage{
make_AOV_volcano(DF, Ctrl, sig_val = 1.5, p_val = 0.05)
}
\arguments{
\item{DF}{The data frame to be used.}

\item{Ctrl}{The condition that everything should be measured against.}

\item{sig_val}{The ratio cutoff (default is 1.5).}

\item{p_val}{The ANOVA p-value. Can also be a regular p-value from a t-test (default is 0.05).}
}
\value{
A pdf file with a volcano plot.
}
\description{
This function takes in a data and produces a volcano plot relative to a specified condition level.
}
\author{
Daniel Braas
}
