% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queryOcc.R
\name{queryOcc}
\alias{queryOcc}
\title{Query Species Occurrences}
\usage{
queryOcc(
  species_name = NULL,
  data_sources = c("gbif", "bien"),
  download_dir = "",
  output_dir = "",
  run_name = "Test",
  mc_cores = NULL
)
}
\arguments{
\item{species_name}{Either a vector of species names or a csv file where species names are stored in the first column.}

\item{data_sources}{A (vector of) character string of the name of the database to search for occurrence data. Must be either `gbif` or `bien` or both for now.}

\item{download_dir}{A directory where to save species occurrence records data downloaded from online databases.}

\item{output_dir}{A directory where to save species occurrence records data after formatting.}

\item{mc_core}{A numeric integer specifying the number of cores to be used for parallel computing. Should be >1 for parallel processing.}
}
\value{
None
}
\description{
Queries species occurrences from different databases online and offline
}
