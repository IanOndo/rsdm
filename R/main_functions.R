#'                                 Run maxent model workflow
#'
#' Uses masked geographically structured cross-validation method for evaluating and selecting a maxent model
#'
#' @param loc_dat A two-column data.frame or matrix, or a path to a csv file with longitude and latitude coordinates of occurrence records.
#' @param envir_dir A path to a directory with the environmental raster layers.
#' @param outputdir A path to a directory where to save the outputs of the species distribution model.
#' @param newdata A character string specifying the directory with the new environmental raster layers for projection. Default is `NULL` i.e. projection will use environmental layers provided for calibration.
#' @param maxent_settings A list of Maxent flags.
#' @param do.map A logical. Should the species distribution be mapped ?  default is `FALSE`.
#' @param do.extrapolation.map A logical. Should extrapolation uncertainty be mapped ? default is `TRUE`.
#' @param do.mic.map A logical. Should the spatially distributed most influential covariate (mic) be mapped ? default is `FALSE`.
#' @param ... Additional parameters to be passed to \code{block_cv_maxent} or \code{make_geographic_block}
#' @return An integer number. `0` in case of success otherwise return `-1`
#' @export
run_maxentModel <- function(loc_dat,
                            envir_dir=NULL,
                            species_name=NULL,
                            outputdir,
                            newdata=NULL,
                            verbose=TRUE,
                            maxent_settings=list(),
                            do.map=FALSE, ...){

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  loc_file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  loc_data_flag = !loc_file_flag & inherits(loc_dat, c("data.frame","matrix"))

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

  list_train_args <- list(loc_dat=loc_dat,
                          species_name=species_name,
                          outputdir=outputdir)

  extra_args <- list(...)
  if(length(extra_args)>0L){
    arg.names <- names(extra_args)

    # check for maxent arguments
    extra_maxent_args <- extra_args[arg.names %in% c("k","coordHeaders","bg_masks","eval_metrics","varying_parameter_name","do.parallel","ncores")]
    if(length(extra_maxent_args)>0L)
      list_train_args <- append(list_train_args, extra_maxent_args)

    # check for training area arguments
    bg_args <- extra_args[arg.names %in% c("land_file","do.alpha_hull","dissolve", "min_area_cover","grid_res", "coastline", "method")]
    if(length(bg_args)>0L)
      list_train_args <- append(list_train_args, bg_args)
  }

  env_dn <- file.path(outputdir,"modelsVariables")
  if(!dir.exists(env_dn)) dir.create(env_dn,recursive=TRUE)

  if(inherits(try(read_layers(env_dn),silent=TRUE),"try-error")){
    if(is.null(envir_dir)){
      #unlink(env_dn,recursive=TRUE)
      stop("Please provide a directory with predictor raster layers")
    }
    rlyr <- read_layers(envir_dir)
    terra::writeRaster(rlyr,
                       file.path(env_dn,paste0(names(rlyr),".tif")),
                       overwrite=TRUE)
  }

  list_train_args[['env_dat']] <- env_dn
  list_train_args <- append(list_train_args,
                            list(maxent_settings=maxent_settings, verbose=verbose))

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
      message(domain)
      stop("Unable to load the training area dataset. See error.")
    }
    out_proj_dn <- file.path(outputdir,"modelsProjectedDomain")
    if(!dir.exists(out_proj_dn)) dir.create(out_proj_dn, recursive=TRUE)
    projection_domain <- make_projection_domain(domain,
                                                output_name=species_name,
                                                output_dir=out_proj_dn,
                                                verbose=verbose)

    if(is.null(projection_domain)){
      # try using biomes
      projection_domain <- make_projection_domain(domain,
                                                  output_name=species_name,
                                                  output_dir=out_proj_dn,
                                                  use_ecoregions=FALSE,
                                                  verbose=verbose)
    }
    if(is.null(projection_domain) || inherits(projection_domain,"try-error")){
      message(projection_domain)
      stop("Unable to create a projection area. See error")
    }

    projection_domain <- terra::vect(projection_domain)

    # prepare new environmental variables
    if(!is.null(newdata) && dir.exists(newdata))
      env_dn <- newdata

    env_bg0 <- try(read_layers(env_dn),silent=TRUE)

    if(inherits(env_bg0,"try-error")){
      message(env_bg0)
      stop("Unable to read environmental raster layers.")
    }

    newdata <-try({
      names(env_bg0) <- if(anyNA(strip_extension(names(env_bg0))))  names(env_bg0) else strip_extension(names(env_bg0))
      env_bg <- terra::crop(env_bg0, projection_domain)
      # create a mask with cells that covers the training block
      prj.dom0 <- terra::rasterize(projection_domain, env_bg, cover=TRUE)
      prj.dom0[prj.dom0==0] <- NA
      terra::mask(env_bg, mask=prj.dom0)
    }, silent=TRUE)

    if(inherits(newdata,"try-error")){
      message(newdata)
      stop("Unable to create newdata layers for projection. See error.")
    }

    out_pred_dn <- file.path(outputdir,"modelsPredictions")
    if(!dir.exists(out_pred_dn)) dir.create(out_pred_dn, recursive=TRUE)
    projected_model <- try(project_maxent(lambdas=trained_model,
                                          newdata=raster::stack(newdata),
                                          quiet=TRUE),
                           silent=TRUE)

    if(!inherits(projected_model, "try-error")){
      if(verbose){
        cat("OK\n");flush.console()
        cat(">...[writing projection]...\n")
      }
      terra::writeRaster(terra::rast(projected_model$prediction_cloglog),
                          file.path(out_pred_dn,paste0(species_name,".tif")),
                          overwrite=TRUE)

      if(do.map){

        if(verbose) cat(">...[mapping projection]...\n")
        out_map_dn <- file.path(outputdir,"modelsMaps")
        if(!dir.exists(out_map_dn)) dir.create(out_map_dn, recursive=TRUE)

        if(loc_file_flag){
          loc_dat = read.csv(loc_dat, header=TRUE)
        }

        id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|^[Xx]$",x = names(loc_dat), value=TRUE)[1]
        id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|^[Yy]$",x = names(loc_dat), value=TRUE)[1]
        coordHeaders <- c(id_x_lon, id_y_lat)
        loc_dat %<>%
          dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
          sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs(4326))

        geo_code <- TDWG::get_geocode(gsub("_", " ", species_name), level=2)
        if(is.null(geo_code))
          inset=NULL
        else
          inset <- TDWG:::tdwg_level2 %>%
          dplyr::filter(LEVEL2_COD %in% geo_code)

        map_raster(x=setNames(projected_model$prediction_cloglog, species_name),
                   y=loc_dat,
                   inset=inset,
                   inset_shape="circular",
                   outputdir=out_map_dn,
                   verbose=verbose)
      }


      if(do.extrapolation.map | do.mic.map){
        olddata <-try({
          names(env_bg0) <- if(anyNA(strip_extension(names(env_bg0))))  names(env_bg0) else strip_extension(names(env_bg0))
          env_bg <- terra::crop(env_bg0, domain)
          # create a mask with cells that covers the training block
          dom0 <- terra::rasterize(domain, env_bg, cover=TRUE)
          dom0[dom0==0] <- NA
          terra::mask(env_bg, mask=dom0)
        }, silent=TRUE)
        do.extra.maps <- !inherits(olddata,"try-error")
      }

      if(!do.extra.maps){
        message(olddata)
        warning("An error occured while processing the training dataset. Unable to create extrapolation uncertainty layers for projection. See error.")
      }

      if(do.extra.maps && do.extrapolation.map){
        if(verbose) cat(">...[writing extrapolation uncertainty map]...\n")
        trg <- raster::brick(newdata)
        ref <- raster::brick(olddata)
        extrap_map  <- try(exdet(trg,ref,compute_mic = do.mic.map),silent=TRUE)
        if(inherits(extrap_map,"try-error")) stop(extra_map)
        out_dirs <- c("modelsPredictionsUncertainty","modelsMIC")[c(do.extrapolation.map,do.mic.map)]
        out_predu_dn <- file.path(outputdir,out_dirs)
        if(any(!dir.exists(out_predu_dn))) sapply(out_predu_dn, dir.create, recursive=TRUE)

        for(k in 1:length(out_predu_dn)){
          terra::writeRaster(terra::rast(extrap_map[[k]]),
                             file.path(out_predu_dn[k],paste0(species_name,".tif")),
                             overwrite=TRUE)
        }

      }else if(do.extra.maps && do.mic.map){
        if(verbose) cat(">...[writing most influential covariate map]...\n")
        trg <- raster::brick(newdata)
        ref <- raster::brick(olddata)
        mic_map <- micdet(trg,ref)
        if(inherits(mic_map,"try-error")) stop(mic_map)

        out_mic_dn <- file.path(outputdir,"modelsMIC")
        if(!dir.exists(out_mic_dn)) dir.create(out_mic_dn, recursive=TRUE)

        terra::writeRaster(terra::rast(mic_map),
                           file.path(out_mic_dn,paste0(species_name,".tif")),
                           overwrite=TRUE)
      }

      if(verbose) cat(">>[[COMPLETED]]<<\n\n");flush.console()
      return(0)
    }else{
      if(verbose) message(projected_model)
      if(verbose) cat("FAILED\n");flush.console()
      if(verbose) cat(">>[[FAILURE]]<<\n\n");flush.console()
      return(-1)
    }
  }else{
    message(trained_model)
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
    r  <- try(terra::raster(bg_domain , res = train_res),silent=TRUE)
    if(inherits(r,"try-error")) {message(r); stop("Cannot create raster from the training area")}
    projected_domain <- suppressMessages(fasterize::fasterize(sf::st_cast(bg_domain,"MULTIPOLYGON"),
                                                              r,
                                                              fun="count"))
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
        message(projected_model)
      }
      return(-1)
    }
  }
  if(verbose){
    cat("FAILED\n\n");flush.console()
    cat(sprintf(">>[[FAILURE]]<<\n\n"));flush.console()
    message(eval_trained_model)
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

