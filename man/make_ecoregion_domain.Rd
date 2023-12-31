% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_geographic_background.R
\name{make_ecoregion_domain}
\alias{make_ecoregion_domain}
\title{Estimating the potential geographic distribution area of a species}
\usage{
make_ecoregion_domain(
  loc_dat,
  coordHeaders = NULL,
  output_dir = NULL,
  output_name = NULL,
  do.alpha_hull = TRUE,
  dissolve = FALSE,
  path_to_alpha_hull = NULL,
  crs_proj = "+init=epsg:4326",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{loc_dat}{A two-column matrix or data.frame or data.table or a path to a csv file with longitude and latitude coordinates.}

\item{coordHeaders}{A character string vector of length two giving the names of the coordinates in \code{loc_dat}}

\item{output_dir}{An optional directory where to save the domain defined.}

\item{output_name}{An optional character string specifying the name to be given to the domain.}

\item{do.alpha_hull}{A logical. Should an alpha hull be built from the set of points ? Default is \code{TRUE}.}

\item{dissolve}{A logical. Should the borders between adjacent polygons be dissolved ? Default is \code{FALSE}.}

\item{path_to_alpha_hull}{A character string specifying the path to the alpha-hull shapefile. Ignored if \code{do.alpha_hull=TRUE}.}

\item{crs_proj}{A character string specifying the Coordinate Reference System of the locations. Default is set to 4326 and should not be changed.}

\item{...}{Parameters to be passed to \code{make_alpha_hulls} function}
}
\value{
An sf object
}
\description{
Works exactly as \code{make_geographic_domain} but is only applied at the ecoregion level (instead of biome + ecoregion level)
}
