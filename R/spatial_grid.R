#'                  Convert a set of point to a grid
#' adapted from https://stackoverflow.com/questions/25547826/generate-regularly-spaced-points-in-polygon
#' @param loc_dat A two-column matrix,a data.frame, a data.table with longitude and latitude coordinates.
#' @param domain A sf object specifying the extent of the grid to build. Default is the bounding box of the location data.
#' @param dst A numeric value specifying the distance in meters between each grid points.Default is 20000m, i.e. 20km.
#' @export
pointsTogrid <- function(loc_dat, domain, dst=20000){
  if(missing(loc_dat))
    stop("loc_dat is missing.")
  if(missing(domain))
    stop("domain is missing.")

  # get bounding box
  bbx <- sf::st_bbox(ext)

  # calculate cell centers
  grid_points <- grid_pixels(loc_dat, ext=bbx, dst=dst)

  # add a z value of ones
  z <- rep(1,nrow(grid_points))

  # return a raster
  raster::rasterFromXYZ(cbind(grid_points,z))
}

#'    Build a regular grid from a set of point
#' modified from https://stackoverflow.com/questions/25547826/generate-regularly-spaced-points-in-polygon
#' @param loc_dat A two-column matrix,a data.frame, a data.table with longitude and latitude coordinates.
#' @param ext A bounding box object specifying the extent of the grid to build. Default is the bounding box of the location data.
#' @param dst A numeric value specifying the distance in meters between each grid points.Default is 20000m, i.e. 20km.
#' @export
grid_pixels <- function(loc_dat, ext=NULL, dst=20000){

  if(!Require("magrittr"))
    stop("Package magrittr needs to be installed to use this function.")

  if(missing(loc_dat))
    stop("loc_dat is missing.")

  if(nrow(loc_dat)< 2)
    stop("loc_dat must have > 2 points.")

  if(ncol(loc_dat)> 2)
    stop("loc_dat has > 2 columns.")

  sf_loc_data <- loc_dat %>%
    sf::st_as_sf(., coords=1:2, crs=sf::st_crs('+proj=longlat +datum=WGS84'))

  if(is.null(ext)){
    # Extract bounding box from
    ext <- sf::st_bbox(sf_loc_data)
  }else{
    if(!inherits(ext,'bbox'))
      stop('ext must be a bbox object.')
  }

  #=====================================
  #=1. Find one of the two closest point
  #=====================================
  distm <- sf::st_distance(sf_loc_data) # calculate the distance matrix
  distm[lower.tri(distm,diag=TRUE)] <- NA # set the diagonal and lower part of the matrix to NA
  ind.orig <- which(distm == min(distm, na.rm=TRUE), arr.ind =TRUE)[,'row'][1]
  pts.orig <- sf::st_coordinates(sf_loc_data)[ind.orig,]

  #=================================================
  #=2. Calculate coordinates along the vertical axis
  #=================================================
  pts <- data.frame()
  #-----------------------------------------
  #=a. number of points on the vertical axis
  #-----------------------------------------
  toN <- pts.orig[2] != ext[['ymax']]
  toS <- pts.orig[2] != ext[['ymin']]
  ymn <- ext[['ymin']]
  ymx <- ext[['ymax']]

  if(toN){ # extends to north (0°)
    ny <- ceiling(geosphere::distGeo(p1 = pts.orig,
                                     p2 = c(pts.orig[1], ext[['ymax']])) / dst)
    pts <- pts %>%
      rbind(geosphere::destPoint(p = pts.orig,
                                 b = 0, d = dst * (1:ny - 1)))
    ymx <- tail(pts,1)[,2]
  }

  if(toS){ # extends to south (180°)
    ny <- ceiling(geosphere::distGeo(p1 = pts.orig,
                                     p2 = c(ext[['xmin']], pts.orig[2])) / dst)
    pts <- pts %>%
      rbind(geosphere::destPoint(p = pts.orig,
                                 b = 180, d = dst * (1:ny - 1)))
    ymn <- tail(pts,1)[,2]
  }
  Ny <- nrow(pts)

  #-------------------------------------------
  #=b. number of points on the horizontal axis
  #-------------------------------------------
  toW <- pts.orig[1] != ext[['xmin']]
  toE <- pts.orig[1] != ext[['xmax']]
  xmn <- ext[['xmin']]
  xmx <- ext[['xmax']]
  # This needs to be calculated for the lowermost and uppermost horizontal lines
  # as the distance between latitudinal lines varies when the longitude changes.
  Nx <- 0
  if(toW){ # extends to west (270°)
    # number of points at the upper limit of the bounding box
    nxLL <- ceiling(geosphere::distGeo(p1 = c(pts.orig[1], ymn),
                                       p2 = c(xmn, ymn)) / dst)
    # number of points at the upper limit of the bounding box
    nxUL <- ceiling(geosphere::distGeo(p1 = c(pts.orig[1], ymx),
                                       p2 = c(xmn, ymx)) / dst)
    nx <- max(nxLL, nxUL)
    # Calculate coordinates.
    for(j in 1:Ny)
      pts <- pts %>%
      rbind(geosphere::destPoint(.[j,], b = 270, dst * 1:(nx - 1)))

    Nx <- Nx + nx
  }

  if(toE){ # extends to east (90°)
    nxLL <- ceiling(geosphere::distGeo(p1 = c(pts.orig[1], ymn),
                                       p2 = c(xmx, ymn)) / dst)
    nxUL <- ceiling(geosphere::distGeo(p1 = c(pts.orig[1], ymx),
                                       p2 = c(xmx, ymn)) / dst)
    nx <- max(nxLL, nxUL)
    # Calculate coordinates.
    for(j in 1:Ny)
      pts <- pts %>%
      rbind(geosphere::destPoint(.[j,], b = 90, dst * 1:(nx - 1)))

    Nx <- Nx + nx
  }

  return(pts)
}

