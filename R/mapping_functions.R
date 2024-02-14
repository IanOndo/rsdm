#'                           Map the species distribution
#'
#'    Plot the raster map of the probability of occurrence of a species
#'
#' @param x A RasterLayer or SpatRaster object, or a path to the raster layer file.
#' @param y An optional sf object describing spatial features to be displayed on the map: points, polygons,...etc.
#' @param inset An optional inset map to be added to the plot. Default is `NULL`.
#' @param outputdir A character string specifying the output directory where the map created will be stored.
#' @param main_title A character string specifying the title of the map.
#' @param legend_title A character string specifying the title of the legend.
#' @param output_format A character string specifying in which format the map should be exported. Default is `pdf` format.
#' @param inset_shape A character string specifying the shape of the thumbnail for the inset map: `rectangular` (default) or `circular`
#' @param logz A logical. Should the raster values be log-transformed ? Default is `FALSE`.
#' @param reproject A logical. Should the raster layer be reprojected. Ignored if argument `projection` and the raster projection are the same.
#' @param projection A character string specifying the new projection in which the raster must be reprojected. Ignored if \code{reproject=FALSE}.
#' @param ... Additional graphic parameters
#' @export
map_raster <- function(x,
                       y=NULL,
                       inset=NULL,
                       outputdir,
                       main_title=NULL,
                       legend_title=NULL,
                       output_format="pdf",
                       inset_shape="rectangular",
                       logz=FALSE,
                       reproject=FALSE,
                       projection="+proj=wintri",
                       verbose=TRUE, ...){

  if(missing(x))
    stop("raster is missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  file_flag = tryCatch(file.exists(x), error=function(err) FALSE) && !tryCatch(dir.exists(x), error=function(err) FALSE)
  data_flag = !file_flag & any(inherits(x, c("RasterLayer","SpatRaster")))

  if(!file_flag & !data_flag)
    stop('Unable to read input raster. Please provide valid data')

  if(!dir.exists(outputdir))
    stop("The directory ",outputdir," does not exist or is not accessible.")

  if(file_flag){
    layer_formats <- "(*.grd$)|(*.asc$)|(*.bil$)|(*.sdat$)|(*.rst$)|(*.tif$)|(*.envi$)|(*.img$)|(*hdr.adf$)"
    if(!grepl(layer_formats, basename(x)))
      stop("Incompatible raster layer format.")
    x <- terra::rast(x)
  }else{
    if(inherits(x,"RasterLayer"))
      x %<>% terra::rast()
  }

  if(is.na(sf::st_is_longlat(x)))
     stop("x has no valid Coordinate Reference System.")

  feature_type = 0
  if(!is.null(y)){
    if(!inherits(y,"sf"))
      stop("y must be an sf object")

    if(all(sf::st_is(y,c("POINT","MULTIPOINT"))))
      feature_type = 1

    if(all(sf::st_is(y,c("POLYGON","MULTIPOLYGON"))))
      feature_type = 2
  }

  add_inset = FALSE
  if(!is.null(inset)){
    if(!inherits(inset,c("sf","sfc")))
      stop("inset must be an sf or sfc object")
    add_inset=TRUE
  }

  bboxx <- suppressMessages(sf::st_bbox(x))
  if(abs(diff(c(bboxx["xmin"],bboxx["xmax"]))) >=100 | abs(diff(c(bboxx["ymin"],bboxx["ymax"]))) >=100)
    reproject=TRUE


  if(reproject){

    if(sf::st_is_longlat(x) & !tryCatch(sf::st_is_longlat(projection),
                                        error=function(err) FALSE)){
      if(verbose) cat(">...[reprojecting raster]...\n")

      # save raster name
      ras_name <- names(x)

      if(pmatch("+proj=wintri", projection, nomatch=0)>0){

        if(inherits(try(x %<>% lwgeom::st_transform_proj(projection),silent=TRUE),"try-error"))
          x %<>%
          terra::project(projection) %>%
          `crs<-`(NA)

          if(!is.null(y))
            y %<>% lwgeom::st_transform_proj(projection)
          if(!is.null(inset))
            inset %<>% lwgeom::st_transform_proj(projection)

        worldbackground_reproj <- RGeodata::terrestrial_lands %>%
          lwgeom::st_transform_proj(projection)
      }else{
        x %<>%
          terra::project(projection)
        if(!is.null(y))
          y %<>% sf::st_transform(projection)
        if(!is.null(inset))
          inset %<>% sf::st_transform(projection)

        worldbackground_reproj <- RGeodata::terrestrial_lands %>%
          lwgeom::st_transform_proj(projection)
      }
    }
    else if(!sf::st_is_longlat(x) & !tryCatch(sf::st_is_longlat(projection),
                                              error=function(err) FALSE)){
      x %<>%
        terra::project(projection)
      if(!is.null(y))
        y %<>% sf::st_transform(projection)
      if(!is.null(inset))
        inset %<>% sf::st_transform(projection)

      worldbackground_reproj <- RGeodata::terrestrial_lands %>%
        lwgeom::st_transform_proj(projection)
    }

    names(x) <- ras_name
  }

  # get species name
  species_name <- if(is.na(strip_extension(names(x)))) unlist(stringr::str_split(names(x),"_")) else unlist(stringr::str_split(strip_extension(names(x)),"_"))

  # define color palette function
  if(verbose) cat(">...[setting color palette function]...\n")
  hcl_settings <- list()
  dot.args <- list(...)
  if(length(dot.args)>0L){
    hcl_extra_settings <- dot.args[names(dot.args) %in% c("n","palette", "alpha", "rev", "fixup")]
    palette_name <- hcl_extra_settings[["palette"]]
    ncolors <- hcl_extra_settings[["n"]]
    if(!is.null(palette_name)){
      if(!palette_name %in% grDevices::hcl.pals())
        stop("Palette ", palette_name, "is not available.")
      if(palette_name %in% rownames(RColorBrewer::brewer.pal.info)){
        # Get maximun number of colors from palette name
        maxColors = RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
        hcl_extra_settings[["n"]] <- min(ncolors, maxColors)
      }
    }else{
      hcl_extra_settings[["palette"]] <- "Spectral"
    }
    if(!is.null(ncolors)){
      if(ncolors<=1){
        warnings("ncolors must be >=1")
        hcl_extra_settings[["n"]] <- 2
      }
    }else{
      hcl_extra_settings[["n"]] <- 3
    }
    hcl_settings <- append(hcl_settings, hcl_extra_settings)
  }else{
    hcl_settings <- list(n=11, palette="Spectral", rev=TRUE)
  }
  color.pal <- function(n){
    hcl_settings[["n"]] <- n
    call <- rlang::call2("hcl.colors",!!!hcl_settings, .ns="grDevices")
    eval(call)
  }

  # set exporting function
  export <- switch(output_format,
                   pdf = grDevices::pdf,
                   png = grDevices::png,
                   jpg = grDevices::jpeg
  )
  if(is.null(export))
    stop("format '",output_format,"' not available")

  if(verbose) cat(">...[setting raster values format]...\n")
  # get raster values in matrix format
  #z <- terra::as.matrix(terra::flip(x), wide=TRUE)
  rotateRaster <- function(x) terra::flip(terra::trans(x), direction="vertical")
  # Very very UGLY !! TODO: Change this ASAP !
  z <- matrix(terra::values(rotateRaster(x)), nrow(x), ncol(x)) %>% t() %>% apply(1,rev) %>% t() %>% apply(2,rev)
  zlevels <- pretty(range(z,finite=TRUE), n=40)

  if(verbose) cat(">...[setting map area]...\n")
  # set up plot area
  bboxx <- suppressMessages(sf::st_bbox(x))
  offset <- if(reproject) 110000 else 1
  brd_bboxx <- sf::st_bbox(c(bboxx["xmin"]-offset, bboxx["xmax"]+offset, bboxx["ymin"]-offset, bboxx["ymax"]+offset))
  if(reproject)
    brd <- sf::st_intersection(worldbackground_reproj,
                               st_rectangle(brd_bboxx, tolerance=0) %>% sf::st_set_crs(sf::st_crs(worldbackground_reproj))) %>%
    suppressMessages() %>%
    suppressWarnings()
  else brd <- sf::st_intersection(RGeodata::terrestrial_lands,
                                  st_rectangle(brd_bboxx, tolerance=0) %>% sf::st_set_crs(sf::st_crs(RGeodata::terrestrial_lands))) %>%
    suppressMessages() %>%
    suppressWarnings()

  if(verbose) cat(">...[mapping]...\n")
  if(is.null(main_title))
    main_title = bquote(italic(.(species_name[1]))~.(species_name[2]))
  if(is.null(legend_title))
    legend_title = 'Probability\nof occurrence'

  # plot
  export(file=file.path(outputdir, paste0(paste(species_name,collapse="_"),'.',output_format)))
  on.exit(grDevices::dev.off())

  par(mar = c(4, 4, 3.5, 6), oma=c(3,3,0,1), cex.main=2, cex.lab=1.2, cex.axis=if(reproject) 0.6 else 1, family='serif')
  plot(sf::st_geometry(brd), reset=FALSE, bgc="#bee8ff",
       axes=TRUE, las=1,
       xaxs="i", yaxs="i", asp=1,
       xlim = c(bboxx["xmin"], bboxx["xmax"]), ylim=c(bboxx["ymin"], bboxx["ymax"]),
       main=title(main=main_title,
                  xlab = if(reproject) expression('Easting') else expression('Longitude'),
                  ylab = if(reproject) expression('Northing') else expression('Latitude')),
       col = "#dddddd",
       border = "#888888",
       lwd = 0.5)
  .filled.contour(x=seq(bboxx["xmin"], bboxx["xmax"], length.out=nrow(z)),
                  y=seq(bboxx["ymin"], bboxx["ymax"], length.out=ncol(z)),
                  z=z,
                  levels= zlevels,
                  col=color.pal(length(zlevels) - 1))
  if(feature_type){
    switch(feature_type,
           {plot(sf::st_geometry(y), pch=3, cex=0.8, add=TRUE)},
           {plot(sf::st_geometry(y), lwd=2, add=TRUE)})
  }
  suppressWarnings(add_legend(zlevels, title=legend_title, color.pal=color.pal))
  if(add_inset){
    if(verbose) cat(">...[creating inset map]...\n")
    add_inset_map(inset, shape=inset_shape)
  }

  if(verbose) cat(">...[completed]...\n\n")
}

#'                           Color bar legend
#'
#'    Add a color bar legend to an existing plot
#'
#' @param zlevels A vector of raster values levels.
#' @param xoffset A numeric betwen 0 and 1 specifying the position of the bar legend along the x axis relative to the plot region.
#' @param yoffset A numeric betwen 0 and 1 specifying the position of the bar legend along the y axis relative to the plot region.
#' @param title A character string specifying the title to be given to the legend
#' @param color.pal A function generating a color palette.
#' @param ... additional graphic parameters
#' @export
add_legend <- function(zlevels, xoffset=0.8, yoffset=-0.3, title='Probability\nof occurrence',color.pal=function(n) grDevices::hcl.colors(n, palette="Spectral", rev=T), ...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')

  zvals <- seq(0, 1, by=0.02)

  xleft=xoffset;
  ybottom=zvals[-length(zvals)]+yoffset;
  xright=xoffset+xoffset/10;
  ytop=zvals[-1L]+yoffset

  # draw legend ticks
  is_valid <- dplyr::between(pretty(zlevels), min(zlevels), max(zlevels))
  ticks_labels <- pretty(zlevels)[is_valid]
  ticks_labels <- replace(ticks_labels, decimalnumcount(ticks_labels)>4, as.numeric(formatC(ticks_labels[decimalnumcount(ticks_labels)>4], format="e", digits=0)))
  at <- (ticks_labels-min(ticks_labels))/(max(ticks_labels)-min(ticks_labels))
  segments(xleft, at+yoffset, xright+1.e-02, at+yoffset)

  # draw legend bar
  rect(xleft, ybottom, xright, ytop, col=color.pal(length(zvals) - 1), border = NA, ...)

  # draw ticks labels
  text(xright+1.e-01, at+yoffset, labels=ticks_labels, adj=c(0.5,0))

  # draw legend title
  text(xright+1.e-02/2, max(ytop)+0.1, labels=bquote(bold(.(title))), adj=c(0.5,0))
}

#'                           Inset map
#'
#'    Add an inset map to an existing plot
#'
#' @param inset A sf object of the map to plot as thumbnail.
#' @param where A character string specifying `where` the inset should be plotted.
#' @param width A numeric value specifying the width of the map.
#' @param heigth A numeric value specifying the height of the map.
#' @export
add_inset_map <- function(inset, where="bottomright", shape="rectangular", width=0.22, height=0.25){

  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  switch(where,
         bottomright={
           xx= 0.77
           yy= 0.21
         },
         bottomleft={
           xx= 0.38
           yy= 0.35
         },
         topright={
           xx= 0.6
           yy= 0.6
         },
         topleftt={
           xx= 0.4
           yy= 0.6
         })
  # set viewport position and size
  grid::pushViewport(grid::viewport(x=xx, y=yy, width = width, height = height))

  suppressWarnings(suppressMessages(inset %<>%
                                      sf::st_buffer(0) %>%
                                      sf::st_union()
  ))
  projected = !sf::st_is_longlat(inset) | is.na(sf::st_is_longlat(inset))
  offset = if(projected) 111000*10 else 10
  inset_bboxx <- suppressMessages(sf::st_bbox(inset))
  inset_brd_bbox <- sf::st_bbox(c(inset_bboxx["xmin"]-offset, inset_bboxx["xmax"]+offset, inset_bboxx["ymin"]-offset, inset_bboxx["ymax"]+offset))

  if(!shape %in% c("circular","rectangular")){
    warning("The shape '", shape,"' is not available. A 'rectangular' shape will be set by default.")
    shape = "rectangular"
  }

  switch(shape,

         circular = {
           radius = max(inset_brd_bbox["xmax"]-inset_brd_bbox["xmin"], inset_brd_bbox["ymax"]-inset_brd_bbox["ymin"])/2
           cent = sf::st_point(c(x = c(inset_brd_bbox["xmax"]+inset_brd_bbox["xmin"])/2, y = c(inset_brd_bbox["ymax"]+inset_brd_bbox["ymin"])/2))
           buf = sf::st_buffer(cent, radius)[[1]] %>%
             list() %>%
             st_polygon() %>%
             sf::st_sfc(crs=sf::st_crs(inset))

           if(projected){
             projection <- if(is.na(sf::st_crs(inset))) "+proj=wintri" else sf::st_crs(inset)
             worldbackground_reproj <- RGeodata::terrestrial_lands %>%
               lwgeom::st_transform_proj(projection)
             inset_brd <- suppressWarnings(
               suppressMessages(sf::st_intersection(worldbackground_reproj, buf))
             )
           }else{
             inset_brd <- suppressWarnings(
               suppressMessages(sf::st_intersection(RGeodata::terrestrial_lands, buf))
             )
           }

         },

         rectangular = {
           if(projected){
             projection <- if(is.na(sf::st_crs(inset))) "+proj=wintri" else sf::st_crs(inset)
             worldbackground_reproj <- RGeodata::terrestrial_lands %>%
               lwgeom::st_transform_proj(projection)
             inset_brd <- suppressWarnings(
               suppressMessages(sf::st_intersection(worldbackground_reproj, st_rectangle(inset_brd_bbox)))
             )
           }else{
             inset_brd <- suppressWarnings(
               suppressMessages(sf::st_intersection(RGeodata::terrestrial_lands, st_rectangle(inset_brd_bbox)))
             )
           }
         }
  )

  # set viewport map
  grid::pushViewport(sf::st_viewport(inset_brd))

  switch(shape,

         circular = {
           suppressWarnings(suppressMessages(buf %>%
                                               sf::st_buffer(2) %>%
                                               sf::st_geometry() %>%
                                               sf::st_as_grob(gp=grid::gpar(fill="#ffffff", lwd=0.1)) %>%
                                               grid::grid.draw())
           )
         },

         rectangular = {
           suppressWarnings(suppressMessages(
             st_rectangle(inset_brd_bbox, tolerance = 0) %>%
               sf::st_buffer(2) %>%
               sf::st_geometry() %>%
               sf::st_as_grob(gp=grid::gpar(fill="#ffffff", lwd=0.1)) %>%
               grid::grid.draw()
           ))
         }
  )

  sf::st_geometry(inset_brd) %>%
    sf::st_as_grob(gp=grid::gpar(fill = "#dddddd",
                                 lwd = 0.1)) %>%
    grid::grid.draw()

  sf::st_geometry(inset) %>%
    sf::st_as_grob(gp=grid::gpar(fill = "#74b374",
                                 lwd=0.1)) %>%
    grid::grid.draw()

}



