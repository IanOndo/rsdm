
#' utility function: read multiple raster layers from a directory.
#' @param path A character string specifying the path to the directory where the raster layers are stored
#' @param nlayers A numeric integer specifying the number of layers to read. Default is all layers.
#' @export
read_layers <- function(path, nlayers=NULL){

  if(!dir.exists(path))
    stop("Unable to find directory ",path,".")

  # Regular expression to check for raster layers with different formats
  layer_formats <- "(*.grd$)|(*.asc$)|(*.bil$)|(*.sdat$)|(*.rst$)|(*.tif$)|(*.envi$)|(*.img$)|(*hdr.adf$)"
  list_layers <- list.files(path, pattern = layer_formats, full.names = TRUE, recursive=TRUE)

  if(length(list_layers)==0L)
    stop("No raster layers found in directory ",path)

  if(is.null(nlayers))
    nlayers = length(list_layers)

  # try with stars
  layers <- try(stars::read_stars(head(list_layers, nlayers)),silent=TRUE)

  if(inherits(layers,"try-error")){

    layers <- lapply(head(list_layers, nlayers), raster::raster) %>%
      raster::stack()

    if(inherits(layers,"try-error")){
      cat(layers)
      stop("Unable to read raster layers. See error message.")
    }
  }

  return(layers)
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


#' From Lon/Lat to UTM zone
#'
#' Converts coordinates from longitude/latitude to UTM zone
#'
#' @param longitude A numeric vector of longitude coordinates in decimal degrees
#' @param latitude A numeric vector of latitude coordinates in decimal degrees
#' @return A numeric vector specifying the UTM zone the coordinates belong to.
#' @export
LL2UTMZone <- function(longitude, latitude, add_latitude_band=FALSE) {

  if(add_latitude_band){

    #  band letters
    band <- LETTERS[3:24]
    band <- band[!band %in% c("I","O")] # letters "I" and "0" are skipped

    # horizontal bands spanning 8 degrees of latitude
    brks <- seq(from=-80, to=84, by=8)
    brks[length(brks)] <- 84 # band 'X' spans 12 degree

    latitude_band <- as.character(cut(latitude, breaks=brks, labels=band))
  }

  # Special zones for Svalbard and Norway
  if(latitude >= 72.0 && latitude < 84.0 )
    if (longitude >= 0.0  && longitude <  9.0)
      return(if(add_latitude_band) paste0(31,latitude_band) else 31);
  if (longitude >= 9.0  && longitude < 21.0)
    return(if(add_latitude_band) paste0(33,latitude_band) else 33)
  if (longitude >= 21.0 && longitude < 33.0)
    return(if(add_latitude_band) paste0(35,latitude_band) else 35)
  if (longitude >= 33.0 && longitude < 42.0)
    return(if(add_latitude_band) paste0(37,latitude_band) else 37)

  if(add_latitude_band) paste0( (floor((longitude + 180) / 6) %% 60) + 1, latitude_band) else (floor((longitude + 180) / 6) %% 60) + 1
}
