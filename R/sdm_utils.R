#' Remove temporary raster files
#' @export
removeTMPFiles <- function(h=24) {

  # remove files in the temp folder that are > h hours old
  warnopt <- getOption('warn')
  on.exit(options('warn'= warnopt))

  tmpdir <- raster::tmpDir(create=FALSE)
  if (!is.na(tmpdir)) {

    d <- raster:::.removeTrailingSlash(tmpdir)
    f <- list.files(path=d, pattern='r_tmp*', full.names=TRUE, include.dirs=TRUE)
    #		f <- list.files(path=d, pattern='[.]gr[di]', full.names=TRUE, include.dirs=TRUE)
    fin <- file.info(f)
    dif <- Sys.time() - fin$mtime
    dif <- as.numeric(dif, units="hours")

    f <- f[which(dif > h)]
    unlink(f, recursive=TRUE)
  }
  options('warn'=warnopt)
}

#'              From meters to decimal degrees
#'
#' Converts distance in meters to distance in decimal degrees given a geometry
#'
#' @param x An sf object representative of a location from where the distance in degree decimal will be calculated. It can be points, polygons,...etc
#' @param dist A numeric value specifying the distance in meters to be converted in degrees
#' @export
meters2degrees <- function(x, dist, ...){

  if(missing(x))
    stop("x is missing.")
  if(missing(dist))
    stop("dist is missing with no default value.")

  if(!inherits(x,"sf"))
    stop("x must be an sf object.")

  if(!is.numeric(dist))
    stop("dist must be numeric")

  if(dist<=0)
    stop("dist must be > 0")

  cent <- suppressWarnings(sf::st_coordinates(sf::st_centroid(st_union(x))))
  degreeLat <- geosphere::destPoint(p=cent, b=0, d=dist,...)[,"lat"] - cent[,"Y"]
  degreeLon <- geosphere::destPoint(p=cent, b=90, d=dist,...)[,"lon"] - cent[,"X"]

  return(list(degreeLat = degreeLat, degreeLon=degreeLon))
}



#'              Most frequent value(s)
#'
#' Find the most frequent n value(s) from a vector
#'
#' @param x A vector
#' @param n A numeric integer specifying the number of values to return. Default is 1, i.e. the most frequent value.
#' @param ... Additional parameters to be passed to the function \code{table}.
#' @seealso \code{table}
#' @export
modal <- function(x, n=1, ...){
  sort(table(x,...),decreasing=TRUE)[1:n]
}

#'              Remove hole(s)
#' from https://stackoverflow.com/questions/52654701/removing-holes-from-polygons-in-r-sf
#' Remove holes from the geometry of a polygon
#'
#' @param x A sf or sfc object.
#' @export
st_remove_holes <- function(x){

  if(!inherits(x,c("sf","sfc")))
    stop("x must be a sf or sfc object.")

  if(any(sf::st_is(x,"MULTIPOLYGON"))){
    x %<>%
      sf::st_cast("POLYGON",warn=FALSE)
  }

  if(inherits(x,"sf")){
    # remove holes
    xgeom <- 1:length(sf::st_geometry(x)) %>%
      lapply(function(p) sf::st_multipolygon(lapply(`[`(sf::st_geometry(x),p), function(y) y[1]))) %>%
      st_as_sfc()
    # set crs
    xgeom %<>%
      sf::st_set_crs(sf::st_crs(x))
    # set new geometry
    x %<>%
      sf::st_set_geometry(xgeom)
  }else{
    if(length(sf::st_geometry(x))>1){
      # remove holes
      xgeom <- 1:length(sf::st_geometry(x)) %>%
        lapply(function(p) sf::st_multipolygon(lapply(`[`(sf::st_geometry(x),p), function(p) y[1]))) %>%
        st_as_sfc()

      #set crs
      xgeom %<>%
        sf::st_set_crs(sf::st_crs(x))

      x<- xgeom
    }else{
      # remove holes
      x %<>%
        purrr::map(function(x) x[1]) %>%
        sf::st_multipolygon()
    }
  }

  return(x)
}

#' Crop/mask raster with a spatial object while accounting for cell coverage
#'
#'
#' crop then mask a raster object with a spatial shape while making sure that boundary cells in contact with the shape are not masked out
#'
#' @param x, A `Raster*` object
#' @param y, A `Spatial*` object
#' @return A `Raster*` object
#' @export
maskCover <- function(x, y, do.crop=TRUE, ...){

  # check inputs
  stopifnot(inherits(x,c("RasterLayer","RasterStack")))
  stopifnot(inherits(y,c("SpatialPolygons","SpatialPolygonsDataFrame","SpatialLines","SpatialPoints")))

  # first match the extent of the vector shape
  if(do.crop)
    stack_reference <- raster::crop(x, y)
  else
    stack_reference <- x

  # then identify cells that at least partially overlap with the shape
  mask_raster <- raster::rasterize(y, stack_reference[[1]], getCover=TRUE)

  # finally mask out cells not overlapping with the shape
  mask_raster[mask_raster==0] <- NA
  raster::stack(raster::mask(stack_reference,mask=mask_raster,...))
}
