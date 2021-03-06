% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_data.R
\name{gather_data}
\alias{gather_data}
\title{This function gathers the individual integration files produced by TraceFinder.}
\usage{
gather_data(
  Title,
  info = info,
  save_raw = TRUE,
  save_dropbox = TRUE,
  drop_vol = drop_vol
)
}
\arguments{
\item{Title}{Project title}

\item{info}{A data frame with information about the metabolomics experiment perfromed.}

\item{save_raw}{(boolean)Should the raw data be saved.}

\item{save_dropbox}{(boolean) Should the raw data be saved on your Dropbox.
This parameter expects that the path to the Dropbox folder is known.}

\item{drop_vol}{The path to the Dropbox project folder to be used.}
}
\value{
A data frame with integrated metabolomics data.
}
\description{
This function gathers the individual integration files produced by TraceFinder.
}
\author{
Daniel Braas
}
