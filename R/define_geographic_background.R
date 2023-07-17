#'                                 Estimating the potential geographic distribution area of a species
#'
#' Delineate a geographic distribution area of a species by combining a set occurrence records locations with biomes and ecoregions
#'
#' @param loc_dat A two-column matrix or data.frame or data.table or a path to a csv file with longitude and latitude coordinates.
#' @param coordHeaders A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param output_dir An optional directory where to save the domain defined.
#' @param output_name An optional character string specifying the name to be given to the domain.
#' @param do.alpha_hull A logical. Should an alpha hull be built from the set of points ? Default is \code{TRUE}.
#' @param dissolve A logical. Should the borders between adjacent polygons be dissolved ? Default is \code{FALSE}.
#' @param path_to_alpha_hull A character string specifying the path to the alpha-hull shapefile. Ignored if \code{do.alpha_hull=TRUE}.
#' @param min_area_cover A numeric value specifying the minimum area cover expressed in percent the point distribution should have relative to the broad background estimate.
#' Calculated as the ratio of the area of the alpha hull to the area of the broad estimate of the distribution area.
#' If this ratio is < `min_area_cover`, then the background is estimated from ecoregions with the function \code{make_ecoregion_domain}.
#' @param min_occ_number A numeric value specifying the minimum number of occurrence records required that must fall in a biome to be selected.
#' @param ... Parameters to be passed to \code{make_alpha_hulls} function
#' @return An sf object
#' @export
make_geographic_domain <- function(loc_dat, coordHeaders=NULL, output_dir=NULL, output_name=NULL, do.alpha_hull=TRUE, dissolve=FALSE, path_to_alpha_hull=NULL, min_area_cover=0.1, min_occ_number=10,verbose=TRUE, ...) {

  if(verbose){
    cat('#==================\n')
    cat('#= 1. Check inputs\n')
    cat('#==================\n\n')
  }
  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  data_flag = !file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!file_flag & !data_flag)
    stop('Unable to read input data. Please provide valid input data')

  if(file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(output_name)){
    if(file_flag){
      output_name=gsub(".csv","",basename(loc_dat))
      if(nchar(output_name)==0L)
        output_name="Unknown"
    }
    else
      output_name="Unknown"
  }

  save_bg <- !is.null(output_dir)
  # ensure that the output directory provided exists
  if(!is.null(output_dir) && !dir.exists(output_dir)){
    warning(paste("Output directory:", output_dir,"does not exist or could not be found."));flush.console()
    save_bg = FALSE
  }

  do.convex_hull = FALSE
  if(!is.null(path_to_alpha_hull)){
    if(!file.exists(path_to_alpha_hull)){
      warning(paste("Path:", path_to_alpha_hull,"does not exist or could not be found."));flush.console()
      do.alpha_hull =TRUE
    }
  }else if (!do.alpha_hull){
    do.convex_hull = TRUE
  }

  if(verbose){
    cat('\n')
    cat('#==================\n')
    cat('#= 2. Make domain \n')
    cat('#==================\n\n')
  }

  if(do.alpha_hull){
    if(verbose) cat('> [...building alpha-hull...]\n')
    dot.args <- list(...)
    if(length(dot.args)>0L){
      dot.args <-dot.args[!names(dot.args) %in% c("land_file","fraction","partCount","initialAlpha","alphaIncrement","default_buffer","maxIter","do.parallel","save.outputs")]
    }
    # try to guess species column
    alpha_hull <- suppressWarnings(try(do.call(AlphaHullRangeModeller::make_alpha_hulls,
                          append(list(loc_data=loc_dat, save.outputs=FALSE, verbose=FALSE), dot.args)),silent=TRUE))
    if(is.null(alpha_hull) | inherits(alpha_hull,"try-error") ){
      cat(alpha_hull)
      stop("Unable to build an alpha hull from this set of points. See error(s).")
    }
  }
  else{
    if(verbose) cat('> [...read alpha-hull from file...]\n')
    alpha_hull <- try(suppressMessages(sf::st_read(path_to_alpha_hull)), silent=TRUE)
    if(is.null(alpha_hull) | inherits(alpha_hull,'try-error') ){
      cat(alpha_hull)
      stop("Unable to read alpha hull from this file specified.")
    }
  }

  if(verbose){
    cat('> [...Done...]')
    cat('\n\n')
    cat('#----------------------------------------------\n')
    cat('#= 2.a. Step 1: Coarse scale domain estimation \n')
    cat('#----------------------------------------------\n')
  }

  if(file_flag){
    loc_dat <- try(read.csv(loc_dat,h=TRUE),silent=TRUE)
    if(inherits(loc_dat,"try-error"))
      stop("Unable to read the csv file of the input dataset.")
    sf_loc_data <- loc_dat %>%
      as.data.frame(row.names=NULL)
  }else{
    sf_loc_data <- loc_dat %>%
      as.data.frame(row.names=NULL)
  }

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|^[Xx]$",x = names(sf_loc_data), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|^[Yy]$",x = names(sf_loc_data), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }
  sf_loc_data <- sf_loc_data %>%
      dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
      dplyr::select(c(X,Y)) %>%
      sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs(alpha_hull))

  if(verbose) cat('> [...keep biomes with at least ',min_occ_number,' occurrence records inside...]\n')

  # merge ecoregions into biomes
  #biomes <- RGeodata::ecoregions %>% dplyr::group_by(BIOME_NAME) %>% dplyr::summarise()

  num.point.by.biome <- suppressMessages(lengths(st_par_intersects(biomes, sf_loc_data, nchunks=100)))
  # keep biomes with more than one occurrence records inside
  biomes.with.more.than.n.point <- biomes[num.point.by.biome >= min_occ_number,]
  if(verbose) cat('> [...extract polygons from biomes containing at least ',min_occ_number,' record(s) and intersecting with the alpha hull...]\n')
  # extract polygons whose biomes contains more than one records and that intersect with the alpha hull (https://rpubs.com/sogletr/sf-ops)
  biomes_more_than_n_point_geom <- suppressMessages(biomes.with.more.than.n.point %>%
                                                        sf::st_cast("POLYGON", warn=FALSE))
  biome_alpha_hull <- alpha_hull %>%
    sf::st_cast("POLYGON", warn=FALSE) %>%
    sf::st_geometry()

  cond <- suppressMessages(st_par_intersects(sf::st_geometry(biomes_more_than_n_point_geom), biome_alpha_hull, nchunks=100) %>%
                      purrr::map_lgl(., function(x) length(x) > 0L))

  occupied_polygons <- biomes_more_than_n_point_geom %>%
    dplyr::filter(cond) %>%
    dplyr::select(BIOME_NAME) %>%
    dplyr::rename(REGION_NAME=BIOME_NAME)

  geographic_domain <- occupied_polygons

  if(verbose) cat('> [...Step 1 [Done]...]\n')

  total_area_cover<- try(as.vector(sum(sf::st_area(alpha_hull))/sum(sf::st_area(geographic_domain))), silent=TRUE)
  if(inherits(total_area_cover,"try-error")){
    alpha_hull.sp = as(alpha_hull,"Spatial")
    geographic_domain.sp = as(geographic_domain,"Spatial")
    total_area_cover <- geosphere::areaPolygon(alpha_hull.sp)/geosphere::areaPolygon(geographic_domain.sp)
  }

  if(total_area_cover < min_area_cover){

    if(verbose) {
      cat('> ... % of area cover of the point distribution < % minimum of area cover specified...\n')
      cat('> ...[making background from ecoregions]...\n')
    }
    bg_args <- list(loc_dat=loc_dat, coordHeaders=coordHeaders, output_dir=output_dir, output_name=output_name, do.alpha_hull=do.alpha_hull, dissolve=dissolve, save.outputs=FALSE, verbose=FALSE)
    dot.args <- list(...)
    if(length(dot.args)>0L){
      arg.names <- names(dot.args)
      if("land_file" %in% arg.names)
        bg_args[["land_file"]] <- dot.args[["land_file"]]
      else
        bg_args[["land_file"]] <- RGeodata::terrestrial_lands
    }
    return(try(do.call(make_ecoregion_domain, bg_args), silent=TRUE))
  }

  is.outside.geographic_domain <- suppressMessages(st_par_intersects(sf_loc_data, geographic_domain, nchunks=10) %>%
                                                 purrr::map_lgl(., function(x) length(x) == 0L ))
  if(sum(is.outside.geographic_domain)==0L){
    if(dissolve){
      if(verbose){
        cat('\n')
        cat('#--------------------------------------------\n')
        cat('#= 2.b Step 2: Merge estimated domain \n')
        cat('#--------------------------------------------\n')
      }
      suppressMessages(
        geographic_domain <- geographic_domain %>%
          sf::st_buffer(0) %>%
          sf::st_union()
      )
      gc()
    }
    geographic_domain #%<>%
      #st_remove_holes()

    if(save_bg){
      saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
    }
    return(geographic_domain)
  }

  if(verbose){
    cat('\n\n')
    cat('#--------------------------------------------\n')
    cat('#= 2.b. Step 2: Fine scale domain estimation \n')
    cat('#--------------------------------------------\n')
  }
  # create working location data
  working_loc_data <- sf_loc_data %>%
    dplyr::filter(is.outside.geographic_domain)

  if(verbose) cat("> [...identify unique records in biomes...]\n")
  # identify unique records in biomes
  biomes.with.min.n.point <- biomes[num.point.by.biome < min_occ_number,]
  is.min.n.point.biome = suppressMessages(st_par_intersects(working_loc_data, biomes.with.min.n.point, nchunks=10) %>%
                                               purrr::map_lgl(., function(x) length(x) == 1L))
  if(any(is.min.n.point.biome)){
    # retrieve ecoregions of records flagged as unique within the biomes
    if(verbose) cat("> [...retrieve ecoregions of records flagged as 'unique' within the biomes...]\n")

    # get the bounding box of points, make sfc
    bounding_box <- sf::st_bbox(working_loc_data[is.min.n.point.biome,]) %>%
      sf::st_as_sfc()

    # identify best ecoregions polygons candidates (intersects)
    best_ecoreg_guess <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split, bounding_box, nchunks=100) %>%
                                            lengths() %>% `>`(.,0L))

    # get ecoregions of records flagged as 'unique' within the biomes
    ecoreg_from_single_record <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split[best_ecoreg_guess, ], working_loc_data[is.min.n.point.biome,], nchunks=100)%>%
                                                    lengths() %>% `>`(.,0L) %>%
                                                    RGeodata::ecoregions_split[best_ecoreg_guess, ][.,] %>%
                                                    dplyr::select(ECO_NAME) %>%
                                                    dplyr::rename(REGION_NAME=ECO_NAME))

    # ecoreg_from_single_record <- suppressMessages(st_par_intersection(RGeodata::ecoregions, working_loc_data[is.min.n.point.biome,], nchunks=100) %>%
    #                                                 dplyr::pull(ECO_NAME) %>%
    #                                                 match(., RGeodata::ecoregions$ECO_NAME) %>%
    #                                                 RGeodata::ecoregions[.,] %>%
    #                                                 dplyr::select(ECO_NAME) %>%
    #                                                 dplyr::rename(REGION_NAME=ECO_NAME))

    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_single_record)
    )

    # remove those records
    working_loc_data %<>%
      dplyr::filter(!is.min.n.point.biome)

    if(length(sf::st_geometry(working_loc_data))==0L){
      if(dissolve){
        if(verbose){
          cat('\n\n')
          cat('> [...Step 2 [Done]...]\n')
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        gc()
        if(verbose) cat('> [...Step 3 [Done]...]\n')
      }
      geographic_domain #%<>%
        #st_remove_holes()

      if(save_bg){
        saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
      }
      gc()
      if(verbose){
        cat('\n\n')
        cat('> [...Step 2 [Done]...]\n')
      }

      return(geographic_domain)
    }

  }else{
    if(verbose) cat("> [...No records found...]\n")
  }

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if(verbose) cat("> [...identify records excluded from the alpha hull...]\n")
  # identify records excluded from the alpha hull
  is.outside.alpha.hull=suppressMessages(sf::st_intersects(working_loc_data, biome_alpha_hull) %>%
                                              purrr::map_lgl(., function(x) length(x) == 0L))

  if(any(is.outside.alpha.hull)){
    # do the same with records excluded from the alpha hull
    if(verbose) cat("> [...retrieve ecoregions of records flagged outside of the alpha hull...]\n")

    # get the bounding box of points, make sfc
    bounding_box <- sf::st_bbox(working_loc_data[is.outside.alpha.hull,]) %>%
      sf::st_as_sfc()

    # identify best ecoregions polygons candidates (intersects)
    best_ecoreg_guess <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split, bounding_box, nchunks=100) %>%
                                            lengths() %>% `>`(.,0L))

    # retrieve ecoregions of records flagged as outside of the alpha hull
    ecoreg_from_alpha_hull_exclusion <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split[best_ecoreg_guess, ], working_loc_data[is.outside.alpha.hull,], nchunks=100)%>%
                                                           lengths() %>% `>`(.,0L) %>%
                                                           RGeodata::ecoregions_split[best_ecoreg_guess, ][.,] %>%
                                                           dplyr::select(ECO_NAME) %>%
                                                           dplyr::rename(REGION_NAME=ECO_NAME))

    # retrieve ecoregions of records flagged as outside of the alpha hull
    # ecoreg_from_alpha_hull_exclusion <- suppressMessages(st_par_intersection(RGeodata::ecoregions, working_loc_data[is.outside.alpha.hull,], nchunks=100) %>%
    #                                                        dplyr::pull(ECO_NAME) %>%
    #                                                        match(., RGeodata::ecoregions$ECO_NAME) %>%
    #                                                        RGeodata::ecoregions[.,] %>%
    #                                                        dplyr::select(ECO_NAME) %>%
    #                                                        dplyr::rename(REGION_NAME=ECO_NAME))

    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_alpha_hull_exclusion)
    )

    # remove those records
    working_loc_data %<>%
      dplyr::filter(!is.outside.alpha.hull)

    if(length(sf::st_geometry(working_loc_data))==0L){

      if(dissolve){
        if(verbose){
          cat('\n\n')
          cat('> [...Step 2 [Done]...]\n')
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        gc()
        if(verbose) cat('> [...Step 3 [Done]...]')
      }
      geographic_domain #%<>%
        #st_remove_holes()

      if(save_bg){
        saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
      }
      gc()
      if(verbose){
        cat('\n\n')
        cat('> [...Step 2 [Done]...]\n')
      }
      return(geographic_domain)
    }

  }else{
    if(verbose) cat("> [...No records found...]\n")
  }

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if(verbose) cat("> [...identify geographic outliers...]\n")
  # identify geographic outliers
  is.geographic.outlier = sf::st_coordinates(working_loc_data) %>%
      data.frame(row.names=NULL) %>%
      dplyr::mutate(kdist = geosphere::distGeo(cbind(X,Y), colMeans(cbind(X,Y)))) %>%
      dplyr::mutate_at(.vars= dplyr::vars(kdist), .funs= list(is.outlier = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %$% is.outlier

  if(any(is.geographic.outlier)){
    if(verbose) cat("> [...retrieve ecoregions of records flagged as geographic outliers...]\n")


    # get the bounding box of points, make sfc
    bounding_box <- sf::st_bbox(working_loc_data[is.geographic.outlier,]) %>%
      sf::st_as_sfc()

    # identify best ecoregions polygons candidates (intersects)
    best_ecoreg_guess <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split, bounding_box, nchunks=100) %>%
                                            lengths() %>% `>`(.,0L))

    # retrieve ecoregions of records flagged as geographic outliers
    ecoreg_from_outliers <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split[best_ecoreg_guess, ], working_loc_data[is.geographic.outlier,], nchunks=100)%>%
                                               lengths() %>% `>`(.,0L) %>%
                                               RGeodata::ecoregions_split[best_ecoreg_guess, ][.,] %>%
                                               dplyr::select(ECO_NAME) %>%
                                               dplyr::rename(REGION_NAME=ECO_NAME))

    # retrieve ecoregions of records flagged as unique within the biomes
    # ecoreg_from_outliers <- suppressMessages(st_par_intersection(RGeodata::ecoregions, working_loc_data[is.geographic.outlier,], nchunks=100) %>%
    #                                            dplyr::pull(ECO_NAME) %>%
    #                                            match(., RGeodata::ecoregions$ECO_NAME) %>%
    #                                            RGeodata::ecoregions[.,] %>%
    #                                            dplyr::select(ECO_NAME) %>%
    #                                            dplyr::rename(REGION_NAME=ECO_NAME))

    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_outliers)
    )

    # remove those records
    working_loc_data %<>%
      dplyr::filter(!is.geographic.outlier)

    if(length(sf::st_geometry(working_loc_data))==0L){
      if(dissolve){
        if(verbose){
          cat('\n\n')
          cat('> [...Step 2 [Done]...]\n')
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        gc()
        if(verbose) cat('> [...Step 3 [Done]...]\n')
      }
      geographic_domain #%<>%
        #st_remove_holes()

      if(save_bg){
        saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
      }

      gc()
      if(verbose){
        cat('\n\n')
        cat('> [...Step 2 [Done]...]\n')
      }

      return(geographic_domain)
    }
  }else{
    if(verbose) cat("> [...No records found...]\n")
    warning('[WARNING] ',length(sf::st_geometry(working_loc_data)),' record(s) remaining outside of the geographic range !')
  }
  geographic_domain #%<>%
    #st_remove_holes()

  if(save_bg){
    saveRDS(geographic_domain, file.path(output_dir, paste0(species_name,".rds")))
  }

  return(geographic_domain)
}

#'                                 Estimating the potential geographic distribution area of a species
#'
#' Works exactly as \code{make_geographic_domain} but is only applied at the ecoregion level (instead of biome + ecoregion level)
#'
#' @param loc_dat A two-column matrix or data.frame or data.table or a path to a csv file with longitude and latitude coordinates.
#' @param coordHeaders A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param output_dir An optional directory where to save the domain defined.
#' @param output_name An optional character string specifying the name to be given to the domain.
#' @param do.alpha_hull A logical. Should an alpha hull be built from the set of points ? Default is \code{TRUE}.
#' @param dissolve A logical. Should the borders between adjacent polygons be dissolved ? Default is \code{FALSE}.
#' @param path_to_alpha_hull A character string specifying the path to the alpha-hull shapefile. Ignored if \code{do.alpha_hull=TRUE}.
#' @param crs_proj A character string specifying the Coordinate Reference System of the locations. Default is set to 4326 and should not be changed.
#' @param ... Parameters to be passed to \code{make_alpha_hulls} function
#' @return An sf object
#' @export
make_ecoregion_domain <- function(loc_dat, coordHeaders=NULL, output_dir=NULL, output_name=NULL, do.alpha_hull=TRUE, dissolve=FALSE, path_to_alpha_hull=NULL, crs_proj="+init=epsg:4326", verbose=TRUE, ...) {

  if(verbose){
    cat('#==================\n')
    cat('#= 1. Check inputs\n')
    cat('#==================\n\n')
  }
  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  data_flag = !file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!file_flag & !data_flag)
    stop('Unable to read input data. Please provide valid input data')

  if(file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(output_name)){
    if(file_flag){
      output_name=gsub(".csv","",basename(loc_dat))
      if(nchar(output_name)==0L)
        output_name="Unknown"
    }
    else
      output_name="Unknown"
  }

  save_bg <- !is.null(output_dir)
  # ensure that the output directory provided exists
  if(!is.null(output_dir) && !dir.exists(output_dir)){
    warning(paste("Output directory:", output_dir,"does not exist or could not be found."));flush.console()
    save_bg = FALSE
  }

  do.convex_hull = FALSE
  if(!is.null(path_to_alpha_hull)){
    if(!file.exists(path_to_alpha_hull)){
      warning(paste("Path:", path_to_alpha_hull,"does not exist or could not be found."));flush.console()
      do.alpha_hull =TRUE
    }
  }else if (!do.alpha_hull){
    do.convex_hull = TRUE
  }

  if(verbose){
    cat('\n')
    cat('#==================\n')
    cat('#= 2. Make domain \n')
    cat('#==================\n\n')
  }

  if(do.alpha_hull){
    if(verbose) cat('> [...building alpha-hull...]\n')
    dot.args <- list(...)
    if(length(dot.args)>0L){
      dot.args <-dot.args[!names(dot.args)=="save.outputs"]
    }
    # try to guess species column
    alpha_hull <- do.call(AlphaHullRangeModeller::make_alpha_hulls,
                          c(list(loc_data=loc_dat, save.outputs=FALSE, verbose=FALSE), dot.args))
    if(is.null(alpha_hull) | inherits(alpha_hull,'try-error'))
      stop("Unable to build an alpha hull from this set of points.")
  }else if(!do.convex_hull){
    if(verbose) cat('> [...read alpha-hull from file...]\n')
    alpha_hull <- try(suppressMessages(sf::st_read(path_to_alpha_hull)), silent=TRUE)
    if(inherits(alpha_hull,'try-error'))
      stop("Unable to read alpha hull from this file specified.")
  }
  if(verbose) cat('> [...Done...]')

  if(verbose){
    cat('\n\n')
    cat('#----------------------------------------------\n')
    cat('#= 2.a. Step 1: Small scale domain estimation \n')
    cat('#----------------------------------------------\n')
  }

  if(file_flag){
    loc_dat <- try(read.csv(loc_dat,h=TRUE),silent=TRUE)
    if(inherits(loc_dat,"try-error"))
      stop("Unable to read the csv file of the input dataset.")
    sf_loc_data <- loc_dat %>%
      as.data.frame(row.names=NULL)
  }else{
    sf_loc_data <- loc_dat %>%
      as.data.frame(row.names=NULL)
  }

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(sf_loc_data), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(sf_loc_data), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }
  sf_loc_data <- sf_loc_data %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    dplyr::select(c(X,Y)) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs(RGeodata::terrestrial_lands))

  if(do.convex_hull){
    if(verbose) cat('> [...building convex-hull...]\n')
    # TODO: rename alpha_hull to hull to be more general
    alpha_hull <- try(sf::st_convex_hull(sf::st_union(sf_loc_data)),silent=TRUE)
    if(inherits(alpha_hull,'try-error')){
      cat(alpha_hull)
      stop("Unable to build convex hull from this set of points.")
    }else{

      # TODO: project alpha hull in easting/northing and reproject back to lon/lat
      # calculate a buffer distance in meters
      buffer_size <- sf_loc_data %>%
        sf::st_coordinates() %>%
        as.data.frame(row.names = NULL) %>%
        AlphaHullRangeModeller::get_OneTenth_distmax()

      # convert the buffer distance in degrees
      dist_in_degrees <- meters2degrees(x=sf_loc_data, dist=buffer_size)

      # buffer the convex_hull
      alpha_hull <- suppressWarnings(suppressMessages(sf::st_buffer(alpha_hull, dist=dist_in_degrees[[which.max(dist_in_degrees)]])))

      # suppressWarnings(
      #     # Reproject to same projection as coastline
      #     alpha_hull = sf::st_transform(alpha_hull, sf::st_crs(RGeodata::terrestrial_lands))
      #     # Intersect
      #     alpha_hull = suppressMessages(sf::st_intersection(alpha_hull, sf::st_geometry(RGeodata::terrestrial_lands)))
      # )
    }
  }

  if(verbose) cat('> [...keep ecoregions with more than one occurrence records inside...]\n')
  num.point.by.ecoregion <- suppressMessages(lengths(st_par_intersects(RGeodata::ecoregions, sf_loc_data, nchunks=100)))
  # keep ecoregions with more than one occurrence records inside
  ecoregions.with.more.than.one.point <- RGeodata::ecoregions[num.point.by.ecoregion > 1,]
  if(verbose) cat('> [...extract polygons whose ecoregions contain more than one record and that intersect with the alpha hull...]\n')
  # extract polygons whose biomes contains more than one records and that intersect with the alpha hull (https://rpubs.com/sogletr/sf-ops)
  ecoregions_more_than_one_point_geom <- RGeodata::ecoregions_split %>%
    dplyr::filter(ECO_NAME %in% unique(ecoregions.with.more.than.one.point$ECO_NAME))

  ecoregion_alpha_hull <- alpha_hull %>%
    sf::st_cast(., "POLYGON", warn=FALSE) %>%
    sf::st_geometry()

  cond <- suppressMessages(st_par_intersects(sf::st_geometry(ecoregions_more_than_one_point_geom), ecoregion_alpha_hull, nchunks=100) %>%
                             purrr::map_lgl(., function(x) length(x) > 0L))

  occupied_polygons <- ecoregions_more_than_one_point_geom %>%
    dplyr::filter(cond) %>%
    dplyr::select(ECO_NAME) %>%
    dplyr::rename(REGION_NAME=ECO_NAME)

  geographic_domain <- occupied_polygons

  if(verbose) cat('> [...Step 1 [Done]...]\n')

  is.outside.geographic_domain <- suppressMessages(sf::st_intersects(sf_loc_data, geographic_domain) %>%
                                                     purrr::map_lgl(., function(x) length(x) == 0L ))
  if(sum(is.outside.geographic_domain)==0L){

    if(dissolve){
      if(verbose){
        cat('\n')
        cat('#--------------------------------------------\n')
        cat('#= 2.b Step 2: Merge estimated domain \n')
        cat('#--------------------------------------------\n')
      }
      suppressMessages(
        geographic_domain <- geographic_domain %>%
          sf::st_buffer(0) %>%
          sf::st_union()
      )
      gc()
      if(verbose) cat('> [...Step 2 [Done]...]\n')
    }
    geographic_domain  %<>%
      sf::st_make_valid()#%<>%
      #st_remove_holes()

    if(save_bg){
      saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
    }
    return(geographic_domain)
  }

  if(verbose){
    cat('\n\n')
    cat('#--------------------------------------------\n')
    cat('#= 2.b. Step 2: Refine scale domain estimation \n')
    cat('#--------------------------------------------\n')
  }
  # create working location data
  working_loc_data <- sf_loc_data %>%
    dplyr::filter(is.outside.geographic_domain)

  if(verbose) cat("> [...identify unique records in ecoregions...]\n")
  # identify unique records in ecoregions
  ecoregions.with.one.point <- RGeodata::ecoregions[num.point.by.ecoregion == 1,]
  is.unique.point.in.ecoregion = suppressMessages(sf::st_intersects(working_loc_data, ecoregions.with.one.point) %>%
                                                purrr::map_lgl(., function(x) length(x) > 0L))
  if(any(is.unique.point.in.ecoregion)){
    # retrieve ecoregions of records flagged as unique within the ecoregions
    if(verbose) cat("> [...retrieve ecoregions of records flagged as unique within the ecoregions...]\n")

    # get the bounding box of points, make sfc
    bounding_box <- sf::st_bbox(working_loc_data[is.unique.point.in.ecoregion,]) %>%
      sf::st_as_sfc()

    # identify best ecoregions candidates (intersects)
    best_ecoreg_guess <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split, bounding_box, nchunks=100) %>%
                                            lengths() %>% `>`(.,0L))

    # get ecoregions of records flagged as 'unique' within the biomes
    ecoreg_from_single_record <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split[best_ecoreg_guess, ], working_loc_data[is.unique.point.in.ecoregion,], nchunks=100)%>%
                                                    lengths() %>% `>`(.,0L) %>%
                                                    RGeodata::ecoregions_split[best_ecoreg_guess, ][.,] %>%
                                                    dplyr::select(ECO_NAME) %>%
                                                    dplyr::rename(REGION_NAME=ECO_NAME))

    # ecoreg_from_single_record <- suppressMessages(st_par_intersection(RGeodata::ecoregions, working_loc_data[is.unique.point.in.ecoregion,], nchunks=100) %>%
    #                                                 dplyr::pull(ECO_NAME) %>%
    #                                                 match(., RGeodata::ecoregions$ECO_NAME) %>%
    #                                                 RGeodata::ecoregions[.,] %>%
    #                                                 dplyr::select(ECO_NAME) %>%
    #                                                 dplyr::rename(REGION_NAME=ECO_NAME))

    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_single_record)
    )

    # remove those records
    working_loc_data %<>%
      dplyr::filter(!is.unique.point.in.ecoregion)

    if(length(sf::st_geometry(working_loc_data))==0L){

      if(dissolve){
        if(verbose){
          cat('\n\n')
          cat('> [...Step 2 [Done]...]\n')
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        gc()
        if(verbose) cat('> [...Step 3 [Done]...]\n')
      }
      geographic_domain %<>%
        sf::st_make_valid() #%<>%
        #st_remove_holes()

      if(save_bg){
        saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
      }
      return(geographic_domain)
    }

  }else{
    if(verbose) cat("> [...No records found...]\n")
  }

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if(verbose) cat("> [...identify records excluded from the alpha hull...]\n")
  # identify records excluded from the alpha hull
  is.outside.alpha.hull=suppressMessages(sf::st_intersects(working_loc_data,ecoregion_alpha_hull) %>%
                                           purrr::map_lgl(., function(x) length(x) == 0))

  if(any(is.outside.alpha.hull)){
    # do the same with records excluded from the alpha hull
    if(verbose) cat("> [...retrieve ecoregions of records flagged outside of the alpha hull...]\n")

    # get the bounding box of points, make sfc
    bounding_box <- sf::st_bbox(working_loc_data[is.outside.alpha.hull,]) %>%
      sf::st_as_sfc()

    # identify best ecoregions polygons candidates (intersects)
    best_ecoreg_guess <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split, bounding_box, nchunks=100) %>%
                                            lengths() %>% `>`(.,0L))

    # retrieve ecoregions of records flagged as outside of the alpha hull
    ecoreg_from_alpha_hull_exclusion <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split[best_ecoreg_guess, ], working_loc_data[is.outside.alpha.hull,], nchunks=100)%>%
                                                           lengths() %>% `>`(.,0L) %>%
                                                           RGeodata::ecoregions_split[best_ecoreg_guess, ][.,] %>%
                                                           dplyr::select(ECO_NAME) %>%
                                                           dplyr::rename(REGION_NAME=ECO_NAME))

    # retrieve ecoregions of records flagged as outside of the alpha hull
    # ecoreg_from_alpha_hull_exclusion <- suppressMessages(st_par_intersection(RGeodata::ecoregions, working_loc_data[is.outside.alpha.hull,], nchunks=100) %>%
    #                                                        dplyr::pull(ECO_NAME) %>%
    #                                                        match(., RGeodata::ecoregions$ECO_NAME) %>%
    #                                                        RGeodata::ecoregions[.,] %>%
    #                                                         dplyr::select(ECO_NAME) %>%
    #                                                         dplyr::rename(REGION_NAME=ECO_NAME))
    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_alpha_hull_exclusion)
    )

    # remode those records
    working_loc_data %<>%
      dplyr::filter(!is.outside.alpha.hull)

    if(length(sf::st_geometry(working_loc_data))==0L){
      if(dissolve){
        if(verbose){
          cat('\n\n')
          cat('> [...Step 2 [Done]...]\n')
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        if(verbose) cat('> [...Step 3 [Done]...]')
      }
      geographic_domain %<>%
        sf::st_make_valid() #%>% #%<>%
        #st_remove_holes()

      if(save_bg){
        saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
      }
      return(geographic_domain)
    }

  }else{
    if(verbose) cat("> [...No records found...]\n")
  }

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if(verbose) cat("> [...identify geographic outliers...]\n")
  # identify geographic outliers
  is.geographic.outlier = sf::st_coordinates(working_loc_data) %>%
    data.frame(row.names=NULL) %>%
    dplyr::mutate(kdist = tryCatch(spatialrisk::haversine(Y, X, mean(Y), mean(X)),
                                   error=function(err) geosphere::distGeo(p1=cbind(X,Y), p2=cbind(mean(X),mean(Y)))) ) %>%
    dplyr::mutate_at(.vars= dplyr::vars(kdist), .funs= list(is.outlier = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %$% is.outlier

  if(any(is.geographic.outlier)){
    if(verbose) cat("> [...retrieve ecoregions of records flagged as geographic outliers...]\n")

    # get the bounding box of points, make sfc
    bounding_box <- sf::st_bbox(working_loc_data[is.geographic.outlier,]) %>%
      sf::st_as_sfc()

    # identify best ecoregions polygons candidates (intersects)
    best_ecoreg_guess <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split, bounding_box, nchunks=100) %>%
                                            lengths() %>% `>`(.,0L))

    # retrieve ecoregions of records flagged as geographic outliers
    ecoreg_from_outliers <- suppressMessages(st_par_intersects(RGeodata::ecoregions_split[best_ecoreg_guess, ], working_loc_data[is.geographic.outlier,], nchunks=100)%>%
                                               lengths() %>% `>`(.,0L) %>%
                                               RGeodata::ecoregions_split[best_ecoreg_guess, ][.,] %>%
                                               dplyr::select(ECO_NAME) %>%
                                               dplyr::rename(REGION_NAME=ECO_NAME))

    # retrieve ecoregions of records flagged as unique within the biomes
    # ecoreg_from_outliers <- suppressMessages(st_par_intersection(RGeodata::ecoregions, working_loc_data[is.geographic.outlier,], nchunks=100) %>%
    #                                            dplyr::pull(ECO_NAME) %>%
    #                                            match(., RGeodata::ecoregions$ECO_NAME) %>%
    #                                            RGeodata::ecoregions[.,] %>%
    #                                            dplyr::select(ECO_NAME) %>%
    #                                            dplyr::rename(REGION_NAME=ECO_NAME))

    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_outliers)
    )

    # remode those records
    working_loc_data %<>%
      dplyr::filter(!is.geographic.outlier)

    if(length(sf::st_geometry(working_loc_data))==0L){
      if(dissolve){
        if(verbose){
          cat('\n\n')
          cat('> [...Step 2 [Done]...]\n')
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        gc()
        if(verbose) cat('> [...Step 3 [Done]...]\n')
      }
      geographic_domain %<>%
        sf::st_make_valid() #%<>%
        #st_remove_holes()

      if(save_bg){
        saveRDS(geographic_domain, file.path(output_dir, paste0(output_name,".rds")))
      }
      return(geographic_domain)
    }
  }else{
    if(verbose) cat("> [...No records found...]\n")
    warning('[WARNING] ',length(sf::st_geometry(working_loc_data)),' records remaining outside of the geographic range !')
  }

  return(geographic_domain)
}

#'                                 Estimating the projection area of a species
#'
#' Extract all ecoregions falling inside the bounding box of a raster or polygon
#'
#' @param bg_dat A spatial object. It can be a raster or a sf surface object (i.e. polygon).
#' @param output_name An optional character string specifying the name to be given to the domain.
#' @param output_dir An optional directory where to save the domain defined.
#' @param dissolve A logical. Should the borders between adjacent polygons be dissolved ? Default is \code{FALSE}.
#' @param use_ecoregions A logical. Should ecoregions be used instead of biomes ? Default is TRUE (i.e. ecoregions are used by default).
#' @return An sf object
#' @export
make_projection_domain <- function(bg_dat, output_name=NULL, output_dir=NULL, dissolve=FALSE, extent_only=FALSE, use_ecoregions=TRUE, verbose=TRUE) {

  if(verbose){
    cat('#==================\n')
    cat('#= 1. Check inputs\n')
    cat('#==================\n\n')
  }
  if(missing(bg_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  file_flag = tryCatch(file.exists(bg_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  data_flag = !file_flag & any(inherits(bg_dat, c("RasterLayer","stars", "sf")))

  if(!file_flag & !data_flag)
    stop('Unable to read input data. Please provide valid input data')

  if(data_flag && inherits(bg_dat,"sf") && !sf::st_dimension(bg_dat)==2)
    stop("bg_dat must be a sf (multi)polygon object.")

  if(is.null(output_name)){
    if(file_flag){
      output_name=strip_extension(basename(bg_dat))
      if(nchar(output_name)==0L)
        output_name="Unknown"
    }
    else
      output_name="Unknown"
  }

  if(verbose){
    cat('\n')
    cat('#==================\n')
    cat('#= 2. Make domain \n')
    cat('#==================\n\n')
  }

  if(file_flag){
    bg_dat <- try(sf::st_read(bg_dat),silent=TRUE)
    if(inherits(bg_dat,"try-error"))
      stop("Unable to read the input dataset bg_dat.")

    if(!sf::st_dimension(bg_dat)==2)
      stop("bg_dat must be a sf (multi)polygon object.")
  }


  if(verbose){
    cat('\n')
    cat('#--------------------------------------------------------------\n')
    cat('#= 2.a Extract polygons falling within the given bounding box\n')
    cat('#--------------------------------------------------------------\n')
  }
  bbx <- suppressWarnings(suppressMessages(if(inherits(bg_dat,"RasterLayer")) st_rectangle(sf::st_bbox(sf::st_as_sf(bg_dat))) else st_rectangle(sf::st_bbox(bg_dat))))
  if(!extent_only){
    if(inherits(bg_dat,"RasterLayer")){
      bg_dat <- bg_dat %>%
      sf::st_as_stars() %>%
      sf::st_as_sf(as_points=FALSE, merge=TRUE)
    }
  }
  suppressMessages(
    if(!use_ecoregions){
      polygons_ <- biomes %>%
        sf::st_cast(., "POLYGON", warn=FALSE)
    }else{
      polygons_ <- RGeodata::ecoregions %>%
        sf::st_cast(., "POLYGON", warn=FALSE)
    }
  )
  # condition 1: polygons must be inside the bounding box of the input spatial object
  cond <- suppressMessages(st_par_within(polygons_, bbx, nchunks=100) %>% lengths() > 0)
  # condition 2: polygons must overlap with the input spatial object
  if(!extent_only){
    cond <- suppressMessages(cond & (st_par_overlaps(polygons_, bg_dat, nchunks=100) %>% lengths() > 0 |
                                     st_par_covered_by(polygons_, bg_dat, nchunks=100) %>% lengths() > 0 ))
  }

  if(sum(cond)>0L){

    inner_polygons <- polygons_ %>%
      dplyr::filter(cond)

    if(!use_ecoregions){
      inner_polygons %<>%
      dplyr::select(BIOME_NAME) %>%
      dplyr::rename(REGION_NAME=BIOME_NAME)
    }else{
      inner_polygons %<>%
        dplyr::select(ECO_NAME) %>%
        dplyr::rename(REGION_NAME=ECO_NAME)
    }
    projection_domain <- inner_polygons

    if(dissolve){
      if(verbose){
        cat('\n')
        cat('#--------------------------------------------\n')
        cat('#= 2.b Step 2: Merge all ecoregions\n')
        cat('#--------------------------------------------\n')
      }
      suppressWarnings(suppressMessages(
        projection_domain <- projection_domain %>%
          sf::st_buffer(0) %>%
          sf::st_union()
      ))
    }
    gc()
    projection_domain %<>%
      sf::st_make_valid() %>%
      st_remove_holes()

    # ensure that the output directory provided exists
    if(!is.null(output_dir)){
      if(!dir.exists(output_dir)){
        warning(paste("Output directory:", output_dir,"does not exist or could not be found."));flush.console()
      }else
        saveRDS(projection_domain, file.path(output_dir,paste0(output_name,".rds")))
    }

    return(projection_domain)

  }else{
    warning("No ecoregions found inside the given bounding box !")
  }
  return(NULL)
}

#' utility function that creates a rectangle (polygon) from a bounding box
#' @param bbx A bbox object
#' @param tolerance A numeric number specifying the buffer distance in decimal degree allowed around the rectangle
#' @export
st_rectangle <- function(bbx, tolerance=0.1){

  if(!inherits(bbx,"bbox")){
    stop("bbx must be a bbox object")
  }

  if(tolerance<0 || tolerance>1)
    stop("tolerance must be between 0 and 1.")
  # outer <- matrix(c(bbx["xmin"], bbx["ymin"], bbx["xmin"], bbx["ymax"],  bbx["xmax"], bbx["ymax"], bbx["xmax"], bbx["ymin"]),ncol=2,byrow=TRUE)
  #
  # rect <- outer %>%
  #   as.data.frame(row.names=NULL) %>%
  #   sf::st_as_sf(coords=1:2, crs=sf::st_crs("+init=epsg:4326")) %>%
  #   sf::st_combine() %>%
  #   sf::st_cast("POLYGON", warn=FALSE) %>%
  #   st_buffer(tolerance)
  xrange <- bbx["xmax"] - bbx["xmin"]
  yrange <- bbx["ymax"] - bbx["ymin"]

  bbx_new <- bbx

  bbx_new[1] <- bbx["xmin"] - (tolerance * xrange) # xmin - left
  bbx_new[2] <- bbx["ymin"] - (tolerance * yrange) # ymin - bottom
  bbx_new[3] <- bbx["xmax"] + (tolerance * xrange) # xmax - right
  bbx_new[4] <- bbx["ymax"] + (tolerance * yrange) # ymax - top

  sf::st_make_grid(bbx_new, n=1)

}

#'                                 Extent the bounding box of a sf point object
#' https://www.jla-data.net/eng/adjusting-bounding-box-of-a-tmap-map/
#' @param loc_dat A sf object or a two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param coordHeaders A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param bbox_ext A numeric value specifying how much the bounding box should be extended in percentage.
#' @export
make_bbox_domain <- function(loc_dat, coordHeaders=NULL, bbox_ext=0.5, min_xrange=0.1666667, min_yrange=0.1666667, left=TRUE, right=TRUE, bottom=TRUE, top=TRUE){

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(is.matrix(loc_dat))
    loc_dat %<>%
    as.data.frame(row.names=NULL)

  if(inherits(loc_dat, "data.frame") & !inherits(loc_dat, "sf")){
    if(is.null(coordHeaders)){
      id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
      id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
      coordHeaders <- c(id_x_lon, id_y_lat)
    }else if(length(coordHeaders)!=2 || all(!coordHeaders %in% names(loc_dat))){
      stop("coordHeaders must be a vector of length 2 of names from loc_dat.")
    }

    loc_dat %<>%
      dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
      sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84'))
  }else if(!inherits(loc_dat,"sf"))
    stop("loc_dat must be a sf, a data.frame or a matrix object.")

  bbox_new <- sf::st_bbox(loc_dat) # current bounding box
  bbox_max<- sf::st_bbox(RGeodata::terrestrial_lands) # maximum bounding box

  xrange <- max(abs(bbox_new$xmax - bbox_new$xmin), min_xrange)# range of x values
  yrange <- max(abs(bbox_new$ymax - bbox_new$ymin), min_yrange)# range of y values

  if(left)
    bbox_new[1] <- min(bbox_new[1] - (bbox_ext * xrange), bbox_max[1]) # xmin - left
  if(right)
    bbox_new[3] <- min(bbox_new[3] + (bbox_ext * xrange), bbox_max[3]) # xmax - right
  if(bottom)
    bbox_new[2] <- min(bbox_new[2] - (bbox_ext * yrange), bbox_max[2]) # ymin - bottom
  if(top)
    bbox_new[4] <- min(bbox_new[4] + (bbox_ext * yrange), bbox_max[4]) # ymax - top

  bbox_new <- bbox_new %>%  # take the bounding box ...
    sf::st_as_sfc() %>%
    sf::st_as_sf() # ... and make it a sf polygon

  return(bbox_new)
}
