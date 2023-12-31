% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping_functions.R
\name{add_legend}
\alias{add_legend}
\title{Color bar legend}
\usage{
add_legend(
  zlevels,
  xoffset = 0.8,
  yoffset = -0.3,
  title = "Probability\\nof occurrence",
  color.pal = function(n) grDevices::hcl.colors(n, palette = "Spectral", rev = T),
  ...
)
}
\arguments{
\item{zlevels}{A vector of raster values levels.}

\item{xoffset}{A numeric betwen 0 and 1 specifying the position of the bar legend along the x axis relative to the plot region.}

\item{yoffset}{A numeric betwen 0 and 1 specifying the position of the bar legend along the y axis relative to the plot region.}

\item{title}{A character string specifying the title to be given to the legend}

\item{color.pal}{A function generating a color palette.}

\item{...}{additional graphic parameters}
}
\description{
Add a color bar legend to an existing plot
}
