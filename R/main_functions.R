#'                                 Run maxent model workflow
#'
#' Uses masked geographically structured cross-validation method for evaluating and selecting a maxent model
#'
#' @param loc_dat A character string specifying the path to a csv file with longitude and latitude coordinates of occurrence records.
#' @param outputdir A character string specifying the directory where to save the raster of interpolated values.
#' @param newdata A character string specifying the directory where the new environmental raster layers for projection are stored.
#' Default is `NULL` i.e. projection will use environmental layers provided for calibration.
#' @param ... Additional parameters to be passed to \code{block_cv_maxent} or \code{make_geographic_block}
#' @return An integer number. `0` in case of success otherwise return `-1`
#' @export
run_maxentModel <- function(loc_dat, species_name=NULL, outputdir, newdata=NULL, verbose=TRUE, maxent_settings=list(), do.map=FALSE, ...){

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  loc_file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  loc_data_flag = !loc_file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!loc_file_flag & !loc_data_flag)
    stop('Unable to read location data. Please provide valid location data')

  if(loc_file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(species_name)){
    if(loc_file_flag){
      species_name=gsub(".csv","",basename(loc_dat),fixed=TRUE)
      if(nchar(species_name)==0L){
        warning("species_name is null. Will assign absolute data-time as species name.")
        species_name=paste0("Unknown_", gsub(" ","_",Sys.time()))
      }
    }
    else
      species_name=paste0("Unknown_", gsub(" ","_",Sys.time()))
  }

  list_train_args <- list(loc_dat=loc_dat, species_name=species_name, outputdir=outputdir)

  extra_args <- list(...)
  if(length(extra_args)>0L){
    arg.names <- names(extra_args)

    # check for maxent arguments
    extra_maxent_args <- extra_args[arg.names %in% c("k","coordHeaders","bg_masks","eval_metrics","varying_parameter_name","do.parallel","ncores")]
    if(length(extra_maxent_args)>0L)
      list_train_args <- append(list_train_args, extra_maxent_args)

    # check for background arguments
    bg_args <- extra_args[arg.names %in% c("land_file","do.alpha_hull","dissolve", "min_area_cover","grid_res", "coastline", "method")]
    if(length(bg_args)>0L)
      list_train_args <- append(list_train_args, bg_args)
  }
  env_dn <- file.path(outputdir,"modelsVariables")
  if(!dir.exists(env_dn))
    stop("Unable to locate directory: ", env_dn)

  list_train_args[['env_dat']] <- env_dn
  list_train_args <- append(list_train_args, list(maxent_settings=maxent_settings))

  if(verbose){
    cat(">> [[ environmentalModel: Maxent ]] <<\n\n")
    cat(sprintf("> ...[Species name: %s]\n", species_name));flush.console()
    cat("> ...[Training]...");flush.console()
  }

  trained_model <- try(do.call(train_Maxent, list_train_args), silent=TRUE)

  if(!inherits(trained_model,"try-error")){
    if(verbose) cat("OK\n");flush.console()

    if(verbose) cat("> ...[Projecting]...");flush.console()
    domain<- try(readRDS(file.path(outputdir,"modelsBackground",paste0(species_name,".rds"))),silent=TRUE)
    if(inherits(domain,"try-error")){
      cat(domain)
      stop("Unable to read background data. See error.")
    }
    out_proj_dn <- file.path(outputdir,"modelsProjectedDomain")
    if(!dir.exists(out_proj_dn)) dir.create(out_proj_dn, recursive=TRUE)
    projection_domain <- make_projection_domain(domain, output_name=species_name, output_dir=out_proj_dn, verbose=FALSE)

    if(is.null(projection_domain)){
      # try using biomes
      projection_domain <- make_projection_domain(domain, output_name=species_name, output_dir=out_proj_dn, use_ecoregions=FALSE, verbose=FALSE)
    }
    if(is.null(projection_domain) || inherits(projection_domain,"try-error")){
      cat(projection_domain)
      stop("Unable to create projection background. See error")
    }
    # prepare new environmental variables
    if(!is.null(newdata) && dir.exists(newdata))
      env_dn <- newdata

    env_bg0 <- try(env_dn %>%
                    read_layers(),silent=TRUE)

    if(inherits(env_bg0,"try-error")){
      cat(env_bg0)
      stop("Unable to read environmental raster layers.")
    }

    newdata <-try({

      names(env_bg0) <- if(any(is.na(strip_extension(names(env_bg0)))))  names(env_bg0) else strip_extension(names(env_bg0))

      if(inherits(env_bg0,"stars")){

        # collapse attributes
        if(length(stars::st_dimensions(env_bg0))<3)
          env_bg0 %<>%
          merge() %>%
          stars::st_set_dimensions(3, values = names(env_bg0)) %>%
          stars::st_set_dimensions(names = c("x", "y", "band"))

        newdata_stars <- suppressMessages( env_bg0 %>%
                                             as("Raster") %>%
                                             raster::stack() %>%
                                             raster::crop(y=as(projection_domain,"Spatial")) %>%
                                             stars::st_as_stars()) # `[`(projection_domain))

        st_dim <- attr(newdata_stars,"dimension")
        # ugly work-around when crop returns crap
        if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){

          env.prj <- as(newdata_stars,"Raster") %>%
            raster::crop(y=as(projection_domain,"Spatial")) %>%
            setNames(stars::st_dimensions(env_bg0)$band$values)

          prj.dom0 <- projection_domain %>%
            as("Spatial") %>%
            raster::rasterize(env.prj, getCover=TRUE)
          prj.dom0[prj.dom0==0] <- NA

          newdata_stars <- env.prj %>%
            raster::mask(mask=prj.dom0) %>%
            raster::stack() # %>%
            #stars::st_as_stars()
          # newdata_stars <- 1:length(newdata_stars) %>%
          #   lapply(function(l){
          #            env.prj <- as(newdata_stars[l],"Raster") %>%
          #              raster::crop(y=as(projection_domain,"Spatial"))
          #
          #            prj.dom0 <- projection_domain %>%
          #              as("Spatial") %>%
          #              raster::rasterize(env.prj, getCover=TRUE)
          #            prj.dom0[prj.dom0==0] <- NA
          #
          #            env.prj %>%
          #            raster::mask(mask=prj.dom0) %>%
          #            stars::st_as_stars()
          #          }
          #     ) %>%
          #   do.call(c,.) %>%
          #   setNames(names(env_bg))
        }else{
          # create a mask with cells that covers the training block
          prj.dom0 <- projection_domain %>%
            as("Spatial") %>%
            raster::rasterize(as(newdata_stars,"Raster"), getCover=TRUE)
          prj.dom0[prj.dom0==0] <- NA

          newdata_stars %>%
            as("Raster") %>%
            setNames(stars::st_dimensions(env_bg0)$band$values) %>%
            raster::mask(mask=prj.dom0)
        }

        # n_layers <- length(newdata_stars)
        #
        # 1:n_layers %>%
        #   as.list() %>%
        #   purrr::map(function(x) as(`[`(newdata_stars,x), "Raster")) %>%
        #   raster::stack()

      }else{

          # set environmental layers extent to projection domain extent
          env_bg <- env_bg0 %>%
            raster::crop(y=as(projection_domain,"Spatial"))

          # create a mask with cells that covers the projection domain
          prj.dom0 <- projection_domain %>%
            as("Spatial") %>%
            raster::rasterize(env_bg, getCover=TRUE)
          prj.dom0[prj.dom0==0] <- NA

          env_bg %>%
          raster::mask(mask=prj.dom0) %>%
          raster::stack()
      }

    }, silent=TRUE)

    if(inherits(newdata,"try-error")){
      cat(newdata)
      stop("Unable to create newdata layers for projection. See error.")
    }

    names(newdata) <- stars::st_dimensions(env_bg0)$band$values

    out_pred_dn <- file.path(outputdir,"modelsPredictions")
    if(!dir.exists(out_pred_dn)) dir.create(out_pred_dn, recursive=TRUE)
    projected_model <- try(project_maxent(lambdas=trained_model, newdata=newdata, quiet=TRUE), silent=TRUE)

    if(!inherits(projected_model, "try-error")){
      if(verbose){
        cat("OK\n");flush.console()
        cat(">...[writing projection]...\n")
      }
      raster::writeRaster(projected_model$prediction_cloglog, file.path(out_pred_dn,paste0(species_name,".tif")), format="GTiff", overwrite=TRUE)

      if(do.map){
        if(verbose) cat(">...[mapping projection]...\n")
        out_map_dn <- file.path(outputdir,"modelsMaps")
        if(!dir.exists(out_map_dn)) dir.create(out_map_dn, recursive=TRUE)

        if(loc_file_flag){
          loc_dat = read.csv(loc_dat, header=TRUE)
          id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
          id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
          coordHeaders <- c(id_x_lon, id_y_lat)
          loc_dat %<>%
            dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
            sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84'))
        }
        geo_code <- TDWG::get_geocode(gsub("_", " ", species_name), level=2)
        if(is.null(geo_code))
          inset=NULL
        else
          inset <- TDWG:::tdwg_level2 %>%
          sf::st_as_sf() %>%
          dplyr::filter(LEVEL2_COD %in% geo_code)

        map_raster(x=setNames(projected_model$prediction_cloglog, species_name), y=loc_dat, inset=inset, inset_shape="circular", outputdir=out_map_dn, verbose=FALSE)
      }
      if(verbose) cat(">>[[COMPLETED]]<<\n\n");flush.console()
      return(0)
    }else{
      if(verbose) cat(projected_model);flush.console()
      if(verbose) cat("FAILED\n");flush.console()
      if(verbose) cat(">>[[FAILURE]]<<\n\n");flush.console()
      return(-1)
    }
  }else{
    cat(trained_model)
    stop("Model training failed. See errors.")
  }

  return(-1)
}

#'                                 Run geographic distance-based model workflow
#'
#' Uses only known occurrence locations to predict the likelihood of finding species in an area surrounding its known presence.
#' The likelihood of finding the species is thus inversely proportional to the distance to the nearest known presence record.
#'
#' @param loc_dat A character string specifying the path to a csv file with longitude and latitude coordinates of occurrence records.
#' @param algorithm A character string specifying the name of the algorithm to use.
#' @param outputdir A character string specifying the directory where to save the raster of interpolated values
#' @param train_res A numeric value specifying the resolution of the grid for the training in decimal degree. Default is 0.01666667 i.e ~ 2km
#' @param project_res A numeric value specifying the resolution of the projected raster in decimal degree. Default is 0.1666667 i.e ~ 20km
#' @return An integer number. `0` in case of success otherwise return `-1`
#' @export
run_geoModel <- function(loc_dat, species_name=NULL, algorithm="idw", outputdir, train_res = 0.01666667, project_res= 0.1666667, verbose=TRUE, ...){

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  if(!is.numeric(train_res)|| is.na(train_res))
    stop("The output resolution must be numeric")

  if(!is.numeric(project_res)|| is.na(project_res))
    stop("The output resolution must be numeric")

  if(is.null(species_name)){
    warning("species_name is null. Will assign absolute data-time as species name.")
    species_name=paste0("Unknown_", gsub(" ","_",Sys.time()))
  }

  if(verbose){
    switch(algorithm,
      idw = cat(">> [[ GeographicModel: Inverse-Distance Weighted ]] <<\n\n"),
      dist = cat(">> [[ GeographicModel: Geographic Distance ]] <<\n\n")
    )
    cat(sprintf("> ...[Species name: %s]\n", species_name));flush.console()
    cat("> ...[Training]...");flush.console()
  }

  train_args <- list(loc_dat=loc_dat, algorithm=algorithm, outputdir=outputdir, species_name=species_name, grid_res = train_res)
  dot.args <- list(...)
  if(length(dot.args)>0L){
    arg.names <- names(dot.args)
    if(any(c("land_file","do.alpha_hull","dissolve") %in% arg.names))
      switch(arg.names,
             land_file = {train_args[["land_file"]] <- dot.args[["land_file"]]},
             do.alpha_hull = {train_args[["do.alpha_hull"]] <- dot.args[["do.alpha_hull"]]},
             dissolve = {train_args[["dissolve"]] <- dot.args[["dissolve"]]}
      )
  }
  eval_trained_model <- try(do.call(train_geoModel, train_args),silent=TRUE)

  if(!inherits(eval_trained_model,"try-error")){
    out_stat_dn <- file.path(outputdir,"modelsStats")
    out_model_fn <- file.path(out_stat_dn, paste0(species_name,".rds"))
    if(!file.exists(out_model_fn))
      stop("unable to locate file: ", out_model_fn)
    if(verbose) cat("OK\n");flush.console()
    trained_model <- readRDS(out_model_fn)$model.object

    if(verbose) cat("> ...[Preparing data for projection]...\n");flush.console()
    out_bg_fn <- file.path(outputdir,"modelsBackground",paste0(species_name,".rds"))
    if(!file.exists(out_bg_fn)){
      cat("No region available for projecting the model")
      return(-1)
    }
    bg_domain <- readRDS(out_bg_fn)
    bg_domain %<>%
      sf::st_buffer(project_res)
    r  <- tryCatch(raster::raster(bg_domain , res = train_res), error=function(err) raster::raster(as(bg_domain,"Spatial"), res = train_res))
    projected_domain <- suppressMessages(fasterize::fasterize(sf::st_cast(bg_domain,"MULTIPOLYGON"), r, fun="count"))
    #r <- raster::raster(bg_domain, res = train_res)
    #projected_domain <- fasterize::fasterize(sf::st_cast(bg_domain,"MULTIPOLYGON"), r, fun="count")

    if(verbose) cat("> ...[Projecting]...");flush.console()
    out_pred_dn <- file.path(outputdir,"modelsPredictions")
    if(!dir.exists(out_pred_dn))
      dir.create(out_pred_dn, recursive=TRUE)

    projected_model <- try(project_geoModel(domain=projected_domain,
                                            model=trained_model,
                                            output_directory=out_pred_dn,
                                            project_res=project_res,
                                            output_name=species_name,
                                            save_output=TRUE), silent=TRUE)

    if(!inherits(projected_model, "try-error")){
      if(verbose){
        cat("OK\n\n");flush.console()
        cat(sprintf(">>[[COMPLETED]]<<\n\n"));flush.console()
      }
      return(0)
    }else{
      if(verbose){
        cat("FAILED\n\n");flush.console()
        cat(sprintf(">>[[FAILURE]]<<\n\n"));flush.console()
        cat(projected_model);flush.console()
      }
      return(-1)
    }
  }
  if(verbose){
    cat("FAILED\n\n");flush.console()
    cat(sprintf(">>[[FAILURE]]<<\n\n"));flush.console()
    cat(eval_trained_model);flush.console()
  }

  return(-1)
}

#'                                 Run point-based model
#'
#' Uses known occurrence locations as its own.
#'
#'
#' @param loc_dat A character string specifying the path to a csv file with longitude and latitude coordinates of occurrence records.
#' @param outputdir A character string specifying the directory where to save the output raster.
#' @param project_res A numeric value specifying the resolution of the output raster in decimal degree. default is 0.1666667 i.e ~ 20km
#' @return An integer number. `0` in case of success otherwise `-1` in case of failure.
#' @export
run_pointModel <- function(loc_dat, species_name=NULL, outputdir, project_res= 0.1666667, verbose=TRUE, ...){

  Require("magrittr")

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  if(!dir.exists(outputdir))
    stop("Unable to find directory :", outputdir)

  if(!is.numeric(project_res)|| is.na(project_res))
    stop("The projecting resolution must be numeric")

  if(is.null(species_name)){
    warning("species_name is null. Will assign absolute data-time as species name.")
    species_name=paste0("Unknown_", gsub(" ","_",Sys.time()))
  }

  dat <- try(read.csv(loc_dat), silent=TRUE)

  if(inherits(dat,"try-error")){
    cat("Unable to read the data")
    return(-1)
  }

  if(verbose){
    cat(">> [[ PointModel ]] <<\n\n")
    cat(sprintf("> ...[Species name: %s]\n", species_name));flush.console()
  }
  #---------------------------------------------------------------------------
  #= 1. get coords (with a regular expression to look for longitude, latitude)
  #---------------------------------------------------------------------------
  x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
  y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
  coordHeaders <- c(x_lon,y_lat)

  if(all(!coordHeaders %in% names(dat))){
    cat("Unable to extract coordinates headers from the data")
    return(-1)
  }

  sf_loc_data <- dat %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84'))

  #---------------------------------------------------------------------------
  #= 2. fasterize
  #---------------------------------------------------------------------------
  list_args <- list(loc_dat=sf_loc_data)
  dot.args <- list(...)
  if(length(dot.args)>0L){
    dot.args <- dot.args[names(dot.args) %in% c("bbox_ext","min_xrange","min_yrange","left","right","bottom","top")]
    list_args <- append(list_args, dot.args)
  }
  if(verbose) cat("> ...[Extracting grid cell location(s) from point(s)]...");flush.console()

  bbox_domain <- do.call(make_bbox_domain, list_args)
  r <- tryCatch(raster::raster(bbox_domain, res = project_res), error=function(err) raster::raster(as(bbox_domain,"Spatial"), res = project_res))
  projected_domain <- raster::rasterize(sf::st_coordinates(sf_loc_data), r)

  if(verbose) cat("OK\n")

  if(verbose) cat(sprintf(">>[[COMPLETED]]<<\n\n"));flush.console()
  #---------------------------------------------------------------------------
  #= 3. save output
  #---------------------------------------------------------------------------
  out_dn_path <- file.path(outputdir,"modelsPredictions")
  if(!dir.exists(out_dn_path))
    dir.create(out_dn_path,recursive=TRUE)

  raster::writeRaster(projected_domain, file.path(out_dn_path, paste0(species_name,".tif")), overwrite=TRUE)

  return(0)
}

