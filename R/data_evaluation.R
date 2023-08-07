#'                                 Train maxent model
#'
#' Calibrate and Evaluate a maxent model using the masked geographically structured cross-validation method for evaluation
#'
#' Partition both occurrence records and background into evaluation bins based on some spatial rules.
#' @param loc_dat A two-column matrix,a data.frame, a data.table or a csv file with longitude and latitude coordinates of occurrence records.
#' @param env_dat A RasterStack object or a character string specifying the path to a directory with environmental raster layers.
#' @param outputdir A character string specifying the directory where to save the outputs of the model.
#' @param coordHeaders [Optional] A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param species_name [Optional] A character string specifying the name of the species being modelled.
#' @param save_bg A logical. Should the background be saved ? Default is `TRUE`.
#' @param ... Additional parameters to be passed to \code{block_cv_maxent} or \code{make_geographic_block}
#' @export
train_Maxent <- function(loc_dat,
                         env_dat,
                         outputdir,
                         coordHeaders=NULL,
                         species_name=NULL,
                         save_bg=TRUE, ...){

  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  if(missing(outputdir))
    stop("Output directory is missing. Please select a directory for maxent model output.")

  loc_file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  loc_data_flag = !loc_file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!loc_file_flag & !loc_data_flag)
    stop('Unable to read location data. Please provide valid location data')

  if(loc_file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(species_name)){
    if(loc_file_flag){
      species_name=strip_extension(basename(loc_dat))
      if(nchar(species_name)==0L)
        species_name="Unknown"
    }
    else
      species_name="Unknown"
  }

  if(loc_file_flag)
    loc_dat = read.csv(loc_dat, header=TRUE)

  if(!dir.exists(outputdir))
    stop("The directory ",outputdir," does not exist or is not accessible.")

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }

  env_dir_flag = tryCatch(dir.exists(env_dat), error=function(err) FALSE)
  env_data_flag = !env_dir_flag & any(inherits(loc_dat, c("RasterStack","RasterLayer")))

  sf_loc_data <- loc_dat %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84')) %>%
    sf::st_coordinates() %>%
    as.data.frame(row.names=NULL)

  #=================================
  #= 2. Train & Evaluate the model #
  #=================================
  outputsubdir<- file.path(outputdir,"modelsDirectories",species_name)
  if(!dir.exists(outputsubdir)) dir.create(outputsubdir, recursive=TRUE)

  #----------------------------------------------------------------
  #= 2.a Create a background area with geographically masked blocks
  #----------------------------------------------------------------
  list_cv_args <- list(loc_dat=sf_loc_data, env_dat=env_dat, outputdir=outputsubdir, species_name=species_name)

  dot.args <- list(...)
  if(length(dot.args)>0L){

    arg.names=names(dot.args)

    bg_args <- dot.args[arg.names %in% c("k", "do.parallel", "land_file", "output_dir", "do.alpha_hull", "dissolve", "min_area_cover", "grid_res", "coastline", "method")]
    if(length(bg_args)>0L){
      bg_masks <- try(do.call(make_geographic_block, append(list(kdat=sf_loc_data, output_name=species_name, verbose=FALSE), bg_args)), silent=TRUE)
      if(inherits(bg_masks,"try-error")){
        cat(bg_masks)
        stop("Unable to create a geographically masked background dataset.")
      }
      list_cv_args <- append(list_cv_args, list(bg_masks=bg_masks,k=dplyr::n_distinct(bg_masks$layer)))
    }else{
      bg_masks <- try(do.call(make_geographic_block, list(kdat=sf_loc_data, output_name=species_name, verbose=FALSE)), silent=TRUE)
      if(inherits(bg_masks,"try-error")){
        cat(bg_masks)
        stop("Unable to create a geographically masked background dataset.")
      }
      list_cv_args <- append(list_cv_args, list(bg_masks=bg_masks,k=dplyr::n_distinct(bg_masks$layer)))
    }
    maxent_cv_args <- dot.args[arg.names %in% c("k", "maxent_settings", "eval_metrics", "varying_parameter_name", "do.parallel", "ncores")]
    if(length(maxent_cv_args)>0L)
      list_cv_args <- c(list_cv_args, maxent_cv_args)

  }else{
    bg_masks <- try(do.call(make_geographic_block, list(kdat=sf_loc_data, output_name=species_name, verbose=FALSE)), silent=TRUE)
    if(inherits(bg_masks,"try-error"))
      stop("Unable to create a geographically masked background dataset.")
    list_cv_args <- append(list_cv_args, list(bg_masks=bg_masks,k=dplyr::n_distinct(bg_masks$layer)))
  }

  if(save_bg){
    out_bg_dn <- file.path(outputdir,"modelsBackground")
    if(!dir.exists(out_bg_dn)) dir.create(out_bg_dn)
    saveRDS(bg_masks, file.path(out_bg_dn, paste0(species_name,".rds")))
  }
  #----------------------------------------------------------------
  #= 2.b Perform block cross-validation
  #----------------------------------------------------------------
  model.output <- try(do.call(block_cv_maxent, list_cv_args),silent=TRUE)

  if(inherits(model.output,"try-error")){
    cat(model.output)
    stop("Model training produced an error.")
  }
  #============================
  #= 3. Select the best model #
  #============================
  metrics <- if(!is.null(list_cv_args[["eval_metrics"]])) list_cv_args[["eval_metrics"]] else c("auc","omission_rate","tss","ic")
  model.performance <- try(get_best_maxent_model(model_output = model.output, eval_metrics = metrics),silent=TRUE)

  if(inherits(model.performance,"try-error")){
    cat(model.performance)
    stop("An error occurred when selecting the best model.")
  }
  # save outputs
  out_stat_dn <- file.path(outputdir,"modelsStats")
  if(!file.exists(out_stat_dn))
    dir.create(out_stat_dn, recursive=TRUE)
  model.summary <- append(list(model_output=model.output), model.performance)

  saveRDS(model.summary, file.path(out_stat_dn, paste0(species_name,".rds")))

  #===========================
  #= 4. train the best model #
  #===========================
  # training occurrence records
  K <- if("k" %in% names(dot.args)) dot.args[["k"]] else 3
  out_model_dn<- file.path(outputsubdir, paste(letters[1:K],collapse=""))
  if(!file.exists(out_model_dn))
    dir.create(out_model_dn, recursive=TRUE)

  training_loc <- data.frame(species=rep(species_name, nrow(sf_loc_data))) %>%
    dplyr::bind_cols(sf_loc_data)
  trainingsampfile = file.path(out_model_dn,'sptrain.csv')
  write.csv(training_loc, trainingsampfile, row.names=FALSE)

  # environmental layers
  if(env_data_flag){
    env_bg0 <- stars::st_as_stars(env_dat)

    # collapse attributes
    if(length(stars::st_dimensions(env_bg0))<3)
      env_bg0 %<>%
      merge() %>%
      stars::st_set_dimensions(3, values = names(env_bg0)) %>%
      stars::st_set_dimensions(names = c("x", "y", "band"))

    env_bg <- env_bg0 %>%
      as("Raster") %>%
      setNames(stars::st_dimensions(env_bg0)$band$values) %>%
      raster::crop(y=as(bg_masks,"Spatial")) %>%
      stars::st_as_stars()

    bg_masks0 <- bg_masks %>%
      as("Spatial") %>%
      raster::rasterize(as(env_bg,"Raster"), getCover=TRUE)
    bg_masks0[bg_masks0==0] <- NA

    training_bg <-try({
      names(env_bg) <- strip_extension(names(env_bg))
      env_bg %>%
        as("Raster") %>%
        setNames(stars::st_dimensions(env_bg0)$band$values) %>%
        raster::mask(mask=bg_masks0) %>%
        stars::st_as_stars()
    },silent=TRUE)

    if(inherits(training_bg,"try-error"))
      stop("Unable to read environmental layers from data env_dat. Please make sure all the layers have the same dimensions.")

  }
  else if(env_dir_flag){
    env_bg0 <- try(env_dat %>%
                    read_layers(),silent=TRUE)

    if(inherits(env_bg0,"try-error")){
      cat(env_bg)
      stop("Unable to read environmental raster layers.")
    }

    training_bg <-try({
      names(env_bg0) <- if(any(is.na(strip_extension(names(env_bg0)))))  names(env_bg0) else strip_extension(names(env_bg0))

      # if(inherits(env_bg,"stars")){
      #   tmp <- env_bg %>%
      #     `[`(bg_masks)
      #
      #   st_dim <- attr(tmp,"dimension")
      #   # ugly work-around when crop returns crap
      #   if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){
      #     1:length(tmp) %>%
      #       lapply(function(l) as(env_bg[l],"Raster") %>%
      #                raster::crop(y=as(bg_masks,"Spatial")) %>%
      #                raster::mask(mask=as(bg_masks,"Spatial")) %>%
      #                stars::st_as_stars()) %>%
      #       do.call(c,.) %>%
      #       setNames(names(env_bg))
      #   }else{
      #     tmp
      #   }
      #
      # }else{
      #   env_bg %>%
      #     raster::crop(y=as(bg_masks,"Spatial")) %>%
      #     raster::mask(mask=as(bg_masks,"Spatial")) %>%
      #     raster::stack()
      # }
      if(inherits(env_bg0,"stars")){

        # collpase attributes
        if(length(stars::st_dimensions(env_bg0))<3)
          env_bg0 %<>%
          merge() %>%
          stars::st_set_dimensions(3, values = names(env_bg0)) %>%
          stars::st_set_dimensions(names = c("x", "y", "band"))

        # set environmental layers extent to training block extent
        tmp <- suppressMessages(env_bg0 %>%
                                  as("Raster") %>%
                                  raster::crop(y=as(bg_masks,"Spatial")) %>%
                                  stars::st_as_stars())

        st_dim <- attr(tmp,"dimension")
        # ugly work-around when crop returns crap
        if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){

          # 1:length(tmp) %>%
          #   lapply(function(l){
          env_bg <- as(env_bg0,"Raster") %>%
            raster::crop(y=as(bg_masks,"Spatial")) %>%
            setNames(stars::st_dimensions(env_bg0)$band$values)

          bg_masks0 <- bg_masks %>%
            as("Spatial") %>%
            raster::rasterize(env_bg, getCover=TRUE)
          bg_masks0[bg_masks0==0] <- NA

          env_bg %>%
            raster::mask(mask=bg_masks0) %>%
            stars::st_as_stars()
          #   }) %>%
          # do.call(c,.) %>%
          # setNames(names(env_bg))
        }else{
          # create a mask with cells that covers the training block
          bg_masks0 <- bg_masks %>%
            as("Spatial") %>%
            raster::rasterize(as(tmp,"Raster"), getCover=TRUE)
          bg_masks0[bg_masks0==0] <- NA

          tmp %>%
            as("Raster") %>%
            setNames(stars::st_dimensions(env_bg0)$band$values) %>%
            raster::mask(mask=bg_masks0) %>%
            stars::st_as_stars()

          # `[`(training_block0 %>%
          #       stars::st_as_stars() %>%
          #       sf::st_as_sf(as_points=FALSE))

        }

      }else{
        # set environmental layers extent to training block extent
        env_bg <- env_bg0 %>%
          raster::crop(y=as(bg_masks,"Spatial"))

        # create a mask with cells that covers the training block
        bg_masks0 <- bg_masks %>%
          as("Spatial") %>%
          raster::rasterize(env_bg, getCover=TRUE)
        bg_masks0[bg_masks0==0] <- NA

        env_bg %>%
          raster::mask(mask=bg_masks0) %>%
          raster::stack()
      }
    })
  }
  if(inherits(training_bg,"try-error"))
    stop("Unable to read environmental layers from ",env_dat,". Please make sure all the layers have the same dimensions.")

  dn <- file.path(out_model_dn,'env')
  if(!dir.exists(dn)) dir.create(dn)

  if(inherits(training_bg,"stars"))
    training_bg %<>%
    split("band") %>%
    setNames(stars::st_dimensions(env_bg0)$band$values)

  # save layers
  n_layers <- if(inherits(training_bg,"RasterStack")) raster::nlayers(training_bg) else length(training_bg)
  for(l in 1:n_layers){
    training_bg_cvrt <- if(inherits(training_bg,"RasterStack"))  training_bg[[l]] else as(training_bg[l],"Raster")
    raster::writeRaster(training_bg_cvrt, file.path(dn, paste0(names(training_bg)[l],'.asc')), format='ascii', overwrite=TRUE)
  }

  # get maximum number of background samples
  maxbackground <- min(sapply(training_bg, function(x) sum(!is.na(x))))

  if(length(dot.args)>0L && "maxent_settings" %in% names(dot.args)){

    maxent_settings<- dot.args[["maxent_settings"]]

    if("varying_param" %in% names(model.summary$best_avg_performance)){

      # get best parameter
      best_param <- model.summary %$%
        best_avg_performance %$%
        varying_param %>%
        gsub(".+?(?<=_)","", ., perl=TRUE)

      if(grepl("[0-9]+",best_param))
        best_param <- as.numeric(best_param)

      varying_parameter_name <- dot.args[["varying_parameter_name"]]
      # set best parameter
      maxent_settings <- replace(maxent_settings, names(maxent_settings)==varying_parameter_name, best_param)
    }

    # set maximum number of background samples
    if("maximumbackground" %in% names(maxent_settings))
      maxent_settings[["maximumbackground"]] <- min(maxent_settings[["maximumbackground"]], maxbackground)
    else
      maxent_settings[["maximumbackground"]] <- min(10000, maxbackground)
  }
  else{
    maxent_settings<- list(maximumbackground=min(10000, maxbackground))
  }

  # check biasfile option
  if(!is.null(maxent_settings[["biasfile"]])){

    biasfile <- maxent_settings[["biasfile"]]
    if(file.exists(biasfile) && !dir.exists(biasfile)){
      bias_file <- try(biasfile %>%
                         raster::raster() %>%
                         raster::crop(y=as(bg_masks,"Spatial")) %>%
                         stars::st_as_stars(),silent=TRUE)

      if(!inherits(bias_file,"try-error")){
        st_dim <- attr(bias_file,"dimension")
        if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){

          biasfile0 <- biasfile %>%
            raster::raster() %>%
            raster::crop(y=as(bg_masks,"Spatial"))

          # create a mask with cells that covers the training block
          bg_masks0 <- bg_masks %>%
            as("Spatial") %>%
            raster::rasterize(biasfile0, getCover=TRUE)
          bg_masks0[bg_masks0==0] <- NA

          bias_file <- biasfile0 %>%
            raster::mask(mask=bg_masks0)
        }else{
          # create a mask with cells that covers the training block
          bg_masks0 <- bg_masks %>%
            as("Spatial") %>%
            raster::rasterize(as(bias_file,"Raster"), getCover=TRUE)
          bg_masks0[bg_masks0==0] <- NA

          bias_file %<>%
            as("Raster") %>%
            raster::mask(mask=bg_masks0)
        }
        # save bias file raster in temporary folder
        outfile <-file.path(raster::tmpDir(),paste0(species_name,"_biasfile",".asc"))
        raster::writeRaster(bias_file, outfile, format="ascii", overwrite=TRUE)
        # set the new biasfile file path
        maxent_settings[["biasfile"]] <- outfile
        maxent_settings[["biastype"]] <- 3
      }else{
        # remove biasfile option(s)
        maxent_settings<- maxent_settings[!names(maxent_settings)=="biasfile"]
      }
    }else{
      # remove biasfile option(s)
      maxent_settings<- maxent_settings[!names(maxent_settings)=="biasfile"]
    }
  }

  model.args <-Reduce(append,list(

    list(sp_points=trainingsampfile,

         env_layers=dn,

         outputdir=out_model_dn),

    maxent_settings,

    list(threads=4, pictures=FALSE, warnings = FALSE)))

  # train the best model
  model_trained <- suppressWarnings(do.call(maxent, model.args))

  if(inherits(model_trained,"try-error")){
    cat(model_trained)
    stop("Maxent model failed. See error.")
  }

  # remove temporary file if needed
  on.exit({
    tmp <- file.path(raster::tmpDir(),paste0(species_name,"_biasfile",".asc"))
    if(file.exists(tmp))
      file.remove(tmp)
  })

  #===========================
  #= 5. return model lambdas #
  #===========================
  return(rmaxent::parse_lambdas(file.path(out_model_dn, paste0(species_name,".lambdas"))))

}

#'                                 Train geographic model
#'
#' Calibrate and evaluate a geographic model using interpoint distance-based algorithm
#'
#' @param loc_dat A two-column matrix,a data.frame, a data.table or a csv file with longitude and latitude coordinates of occurrence records.
#' @param algorithm A character string specifying the name of the algorithm to use. Default is `idw` i.e. Inverse-Distance Weighted.
#' @param coordHeaders [Optional] A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param outputdir A character string specifying the directory where to save the outputs of the model.
#' @param species_name [Optional] A character string specifying the name of the species being modelled.
#' @param grid_res A numeric value specifying the resolution of the output raster layer in decimal degree. Default is 0.01666667 (~ 2 km).
#' @param save_bg A logical. Should the background be saved ? Default is `TRUE`.
#' @export
train_geoModel <- function(loc_dat,
                           bg_dat=NULL,
                           algorithm="idw",
                           coordHeaders=NULL,
                           outputdir,
                           species_name=NULL,
                           grid_res = 0.01666667, save_bg=TRUE, ...){

  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  if(missing(outputdir))
    stop("Output directory is missing. Please select a directory for inverse-distance weighted model output.")

  loc_file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  loc_data_flag = !loc_file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!loc_file_flag & !loc_data_flag)
    stop('Unable to read location data. Please provide valid location data')

  if(loc_file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(species_name)){
    if(loc_file_flag){
      species_name=gsub(".csv","",basename(loc_dat))
      if(nchar(species_name)==0L)
        species_name="Unknown"
    }
    else
      species_name="Unknown"
  }

  if(loc_file_flag)
    loc_dat = read.csv(loc_dat, header=TRUE)

  if(!dir.exists(outputdir))
    stop("The directory ",outputdir," does not exist or is not accessible.")

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }

  sf_loc_data <- loc_dat %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84')) %>%
    sf::st_coordinates()

  #==========================
  #= 1. Make the background #
  #==========================
  if(is.null(bg_dat)){

    bg_args <- list(loc_dat=sf_loc_data, do.alpha_hull=TRUE, dissolve=FALSE, land_file=RGeodata::terrestrial_lands, save.outputs=FALSE, verbose=FALSE)
    dot.args <- list(...)
    if(length(dot.args)>0L){
      arg.names <- names(dot.args)
      if(any(c("land_file","do.alpha_hull","dissolve") %in% arg.names))
        switch(arg.names,
               land_file = {bg_args[["land_file"]] <- dot.args[["land_file"]]},
               do.alpha_hull = {bg_args[["do.alpha_hull"]] <- dot.args[["do.alpha_hull"]]},
               dissolve = {bg_args[["dissolve"]] <- dot.args[["dissolve"]]}
        )
    }
    invisible(capture.output(domain <- try(do.call(make_ecoregion_domain, bg_args), silent=TRUE)))

    if(inherits(domain,'try-error')){
      cat(domain)
      stop("Unable to build a background from the set of points")
    }

    if(save_bg){
      out_bg_dn <- file.path(outputdir,"modelsBackground")
      if(!dir.exists(out_bg_dn)) dir.create(out_bg_dn)
      saveRDS(domain, file.path(out_bg_dn, paste0(species_name,".rds")))
    }

    # get background coordinates
    raster_template  <- tryCatch(raster::raster(domain, res = grid_res), error=function(err) raster::raster(as(domain,"Spatial"), res = grid_res))
    raster_template[] <- 1
    # get presence cells
    pres_cells <- raster::cellFromXY(raster_template, xy=sf_loc_data)
    all_cells <- raster::extract(x=raster_template, y=as(domain,"Spatial"),cellnumber=TRUE,df=TRUE)[['cell']]
    # get background cells
    bg_cells <- base::setdiff(all_cells,pres_cells)
    # sample background cells to achieve at least 0.1 of prevalence
    samp_cells <- sample(bg_cells, size=min(length(unique(pres_cells))*10, length(all_cells)))
    # get coordinates
    bg_dat <- raster::xyFromCell(raster_template, cell=samp_cells)
    bg_dat<- na.omit(bg_dat)
  }

  #======================
  #= 2. Train the model #
  #======================
  geo_model <- switch(algorithm,
                      idw  = try(dismo::geoIDW(p=sf_loc_data, a=bg_dat), silent=TRUE),
                      dist = try(dismo::geoDist(p=sf_loc_data, lonlat=TRUE), silent=TRUE)
  )

  if(inherits(geo_model,"try-error")){
    cat(geo_model)
    stop("Unable to train the geographic model: ", algorithm)
  }

  #=========================
  #= 3. Evaluate the model #
  #=========================
  colnames(bg_dat) <- colnames(geo_model@presence)
  list_args <- list(pr=sf_loc_data, bg=bg_dat, algorithm=algorithm)
  dot.args <- list(...)
  if(length(dot.args)>0L){
    if("tr" %in% names(dot.args))
      list_args <- append(list_args, dot.args[names(dot.args)=="tr"])
  }
  # leave-one-out cross-validation
  eval_dist_model <- do.call(loocv_geo, list_args)

  # save outputs
  out_stat_dn <- file.path(outputdir,"modelsStats")
  if(!file.exists(out_stat_dn))
    dir.create(out_stat_dn, recursive=TRUE)
  saveRDS(list(model.object=geo_model, model.eval=eval_dist_model), file.path(out_stat_dn, paste0(species_name,".rds")))

  return(eval_dist_model)
}

#'                                 Cross-validation with geographically masked blocks of occurrence records for evaluation
#'
#' Partition both occurrence records and background into evaluation bins based on some spatial rules.
#'
#'
#' @param loc_dat A two-column matrix,a data.frame, a data.table or a csv file with longitude and latitude coordinates of occurrence records.
#' @param env_dat A RasterStack object or a character string specifying the path to a directory with environmental raster layers.
#' @param k A numeric integer specifying the number of bins required.
#' @param coordHeaders A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param bg_masks A RasterStack object of the background to be partitionned.
#' @param species_name [Optional] A character string specifying the name of the species. If \code{NULL}, will be set to `Unknown`.
#' @param maxent_settings A list of maxent flags.
#' @param eval_metrics A vector of evaluation metric names. Evaluation metrics available are: `auc`, `omission_rate`,`tss`, `ic`.
#' @param do.mask A logical
#' @param ... Additional parameters to be passed to \code{make_geographic_block}
#' @param return A data.frame reporting the evaluation metrics for each k bins.
#' @export
block_cv_maxent <- function(loc_dat, env_dat,
                            k=3,
                            coordHeaders=NULL,
                            bg_masks=NULL,
                            outputdir,
                            species_name=NULL,
                            maxent_settings=list(),
                            eval_metrics=c("auc","omission_rate","tss","ic"),
                            varying_parameter_name=NULL,
                            do.mask=TRUE, do.parallel=FALSE, ncores=NULL, ...){

  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  if(missing(env_dat))
    stop("Please provide a path to your environmental layers.")

  if(missing(outputdir))
    stop("Output directory is missing. Please select a directory for maxent outputs.")

  loc_file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  loc_data_flag = !loc_file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!loc_file_flag & !loc_data_flag)
    stop('Unable to read location data. Please provide valid location data')

  if(loc_file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(species_name)){
    if(loc_file_flag){
      species_name=gsub(".csv","",basename(loc_dat))
      if(nchar(species_name)==0L)
        species_name="Unknown"
    }
    else
      species_name="Unknown"
  }

  if(loc_file_flag)
    loc_dat = read.csv(loc_dat, header=TRUE)

  if(!dir.exists(outputdir))
    stop("The directory ",outputdir," does not exist or is not accessible.")

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }

  env_dir_flag = tryCatch(dir.exists(env_dat), error=function(err) FALSE)
  env_data_flag = !env_dir_flag & any(inherits(loc_dat, c("RasterStack","RasterLayer")))

  if(!env_dir_flag & !env_data_flag)
    stop('Unable to read environmental data. Please provide valid environmental data.')

  if(is.matrix(loc_dat))
    loc_dat %<>%
    as.data.frame(row.names=NULL)

  sf_loc_data <- loc_dat %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84'))

  maxent_args = bg_args = list()
  dot.args <- list(...)
  if(length(dot.args)>0L){
    arg.names <- names(dot.args)
    if(!any(c("path_to_maxent","path_to_java") %in% arg.names))
      bg_args <- dot.args[!arg.names%in%c("path_to_maxent","path_to_java")]
    if(any(c("path_to_maxent","path_to_java") %in% arg.names))
      maxent_args <- dot.args[arg.names%in%c("path_to_maxent","path_to_java")]
  }
  if(length(maxent_args)>0L)
    maxent_settings <- append(maxent_settings, list(maxent_args))

  varying_parameter = FALSE
  if(!is.null(varying_parameter_name)){
    varying_parameter = TRUE
    if(!varying_parameter_name %in% names(maxent_settings)){
      warning(varying_parameter_name, "is not a valid maxent parameter name or maxent_settings list is empty !")
      varying_parameter = FALSE
    }
    if(length(maxent_settings[[varying_parameter_name]])<2){
      warning("The parameter ",varying_parameter_name, " contains one single value or maxent_settings list is empty !")
      varying_parameter = FALSE
    }
  }else{
    maxent_settings <- replace(maxent_settings, lengths(maxent_settings)>1, maxent_settings[lengths(maxent_settings)>1][[1]][1])
    varying_parameter_name=""
  }

  if(!is.null(bg_masks)){
    if(!inherits(bg_masks,c('RasterStack','sf')))
      stop('Argument bg_masks must be a RasterStack or sf object.')
    if(inherits(bg_masks,'RasterStack')){
      if(raster::nlayers(bg_masks)!=k)
      stop('Number of layers in bg_masks must match the number of cluster k.')
    }else{
      if(length(unique(bg_masks %$% layer)) !=k)
        k <- length(unique(bg_masks %$% layer)) #stop('Number of cluster in bg_masks must match the number of cluster k.')
    }

  }
  else {
    #==================================
    #= 0. Create a background area with
    #     geographically masked blocks
    #==================================
    if(length(bg_args)>0L)
      bg_masks <- try(do.call(make_geographic_block, append(list(kdat=sf::st_coordinates(sf_loc_data), k=k), bg_args)), silent=TRUE)
    else
      bg_masks <- try(make_geographic_block(kdat=sf::st_coordinates(sf_loc_data), k=k), silent=TRUE)

    if(inherits(bg_masks,"try-error")){
      cat(bg_masks)
      stop("Unable to create a geographically masked background dataset.")
    }

  }

  #================================================
  #= 1. Create directories for each block partition
  #================================================
  comb <- combn(1:k,k-1)
  k_fold_names <- apply(comb,2,function(x) paste(letters[x],collapse="_"))

  dirs_to_create 		<- file.path(outputdir,k_fold_names)
  new.dirs.created 	<- sapply(dirs_to_create, function(x) if(!dir.exists(x)) dir.create(x) else TRUE)

  if(!any(new.dirs.created))
    stop("Unable to write output directories.")

  #===========================================
  #= 2. Cross-validate the model across blocks
  #===========================================
  on.exit({
    doParallel::stopImplicitCluster()
    unregister()
  })
  if(do.parallel){
    # Parallel processing
    ncores = if(is.null(ncores)) parallel::detectCores()-1 else min(ncores,parallel::detectCores()) # number of cores to use
    if(inherits(doParallel::registerDoParallel(ncores),"error"))
      doParallel::registerDoSEQ()
  }else{
    # Sequential processing
    foreach::registerDoSEQ()
  }

  cv_output <- vector(mode="list", length=k)

  cv_output <- tryCatch({

    foreach::foreach(j=1:k, .packages=c("rsdm","rmaxent","dplyr"), .combine=rbind) %dopar% {

    #--------------------------
    #= a. get the training data
    #--------------------------

    # training occurrence records
    training_block <- bg_masks %>%
      subset(layer %in% comb[,j]) %>%
      sf::st_make_valid()

    is_training_loc <- suppressMessages(sf_loc_data %>%
                                          sf::st_intersects(training_block) %>%
                                          purrr::map_lgl(., function(x) length(x) > 0L))
    training_loc <- data.frame(species=rep(species_name, sum(is_training_loc))) %>%
      dplyr::bind_cols(as.data.frame(sf::st_coordinates(sf_loc_data[is_training_loc,]))) %>%
      na.omit()

    # save samples
    trainingsampfile = file.path(dirs_to_create[j],'sptrain.csv')
    write.csv(training_loc, trainingsampfile, row.names=FALSE)

    # environmental layers
    if(env_data_flag){
      env_bg0 <- stars::st_as_stars(env_dat)

      # collapse attributes
      if(length(stars::st_dimensions(env_bg0))<3)
        env_bg0 %<>%
        merge() %>%
        stars::st_set_dimensions(3, values = names(env_bg0)) %>%
        stars::st_set_dimensions(names = c("x", "y", "band"))

      env_bg <- env_bg0 %>%
        as("Raster") %>%
        setNames(stars::st_dimensions(env_bg0)$band$values) %>%
        raster::crop(y=as(training_block,"Spatial")) %>%
        stars::st_as_stars()

      training_block0 <- training_block %>%
        as("Spatial") %>%
        raster::rasterize(as(env_bg,"Raster"), getCover=TRUE)
      training_block0[training_block0==0] <- NA

      training_bg <-try({
        names(env_bg) <- strip_extension(names(env_bg))
        env_bg %>%
          as("Raster") %>%
          setNames(stars::st_dimensions(env_bg0)$band$values) %>%
          raster::mask(mask=training_block0) %>%
          stars::st_as_stars()
      },silent=TRUE)
      if(inherits(training_bg,"try-error"))
        stop("Unable to read environmental layers from data env_dat. Please make sure all the layers have the same dimensions.")
    }
    else if(env_dir_flag){

      env_bg0 <- try(env_dat %>%
                      read_layers(),silent=TRUE)
      if(inherits(env_bg0,"try-error")){
        cat(env_bg0)
        stop("Unable to read environmental raster layers.")
      }

      training_bg <-try({

        names(env_bg0) <- if(any(is.na(strip_extension(names(env_bg0)))))  names(env_bg0) else strip_extension(names(env_bg0))

        if(inherits(env_bg0,"stars")){

          # collpase attributes
          if(length(stars::st_dimensions(env_bg0))<3)
            env_bg0 %<>%
            merge() %>%
            stars::st_set_dimensions(3, values = names(env_bg0)) %>%
            stars::st_set_dimensions(names = c("x", "y", "band"))

          # set environmental layers extent to training block extent
          tmp <- suppressMessages(env_bg0 %>%
                                    as("Raster") %>%
                                    raster::crop(y=as(training_block,"Spatial")) %>%
                                    stars::st_as_stars())

          st_dim <- attr(tmp,"dimension")
          # ugly work-around when crop returns crap
          if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){

            env_bg <- as(env_bg0,"Raster") %>%
                        raster::crop(y=as(training_block,"Spatial")) %>%
                        setNames(stars::st_dimensions(env_bg0)$band$values)

                       training_block0 <- training_block %>%
                         as("Spatial") %>%
                         raster::rasterize(env_bg, getCover=TRUE)
                       training_block0[training_block0==0] <- NA

                       env_bg %>%
                       raster::mask(mask=training_block0) %>%
                       stars::st_as_stars()

          }else{
            # create a mask with cells that covers the training block
            training_block0 <- training_block %>%
              as("Spatial") %>%
              raster::rasterize(as(tmp,"Raster"), getCover=TRUE)
            training_block0[training_block0==0] <- NA

            tmp %>%
              as("Raster") %>%
              setNames(stars::st_dimensions(env_bg0)$band$values) %>%
              raster::mask(mask=training_block0) %>%
              stars::st_as_stars()

          }

        }else{
          # set environmental layers extent to training block extent
          env_bg <- env_bg0 %>%
            raster::crop(y=as(training_block,"Spatial"))

          # create a mask with cells that covers the training block
          training_block0 <- training_block %>%
            as("Spatial") %>%
            raster::rasterize(env_bg, getCover=TRUE)
          training_block0[training_block0==0] <- NA

          env_bg %>%
            raster::mask(mask=training_block0) %>%
            raster::stack()
        }

      }, silent=TRUE)

      if(inherits(training_bg,"try-error"))
        stop("Unable to read environmental layers from ",env_dat,". Please make sure all the layers have the same dimensions.")
    }
    dn <- file.path(dirs_to_create[j],'env')
    if(!dir.exists(dn)) dir.create(dn)

    if(inherits(training_bg,"stars"))
      training_bg %<>%
      split("band") %>%
      setNames(stars::st_dimensions(env_bg0)$band$values)

    # save layers
    n_layers <- if(inherits(training_bg,"RasterStack")) raster::nlayers(training_bg) else length(training_bg)
    for(l in 1:n_layers){
      training_bg_cvrt <- if(inherits(training_bg,"RasterStack"))  training_bg[[l]] else as(training_bg[l],"Raster")
      raster::writeRaster(training_bg_cvrt, file.path(dn, paste0(names(training_bg)[l],'.asc')), format='ascii', overwrite=TRUE)
    }

    #-------------------------
    #= b. get the testing data
    #-------------------------
    testing_loc <- data.frame(species=rep(species_name, sum(!is_training_loc))) %>%
      dplyr::bind_cols(as.data.frame(sf::st_coordinates(sf_loc_data[!is_training_loc,])))

    # add environmental variables
    if(inherits(env_bg0,"stars")){
      env_bg0 %<>%
        split("band") %>%
        setNames(stars::st_dimensions(env_bg0)$band$values)
      testing_loc <- testing_loc %>%
        dplyr::bind_cols(map_lfd(env_bg0,function(x) raster::extract(as(x,"Raster"),.[,-1]))) %>%
        na.omit()
    }else{
      testing_loc <- testing_loc %>%
        dplyr::bind_cols(raster::extract(env_bg0,.[,-1]) %>% as.data.frame(row.names = NULL)) %>%
        na.omit()
    }

    # check biasfile option
    if(!is.null(maxent_settings[["biasfile"]])){
      biasfile <- maxent_settings[["biasfile"]]
      if(file.exists(biasfile) && !dir.exists(biasfile)){
        bias_file <- try(biasfile %>%
                           raster::raster() %>%
                           raster::crop(y=as(training_block,"Spatial")) %>%
                           stars::st_as_stars(),silent=TRUE)

        if(!inherits(bias_file,"try-error")){
          st_dim <- attr(bias_file,"dimension")
          if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){

            biasfile0 <- biasfile %>%
              raster::raster() %>%
              raster::crop(y=as(training_block,"Spatial"))

            # create a mask with cells that covers the training block
            training_block0 <- training_block %>%
              as("Spatial") %>%
              raster::rasterize(biasfile0, getCover=TRUE)
            training_block0[training_block0==0] <- NA

            bias_file <- biasfile0 %>%
              raster::mask(mask=training_block0)
          }else{
            # create a mask with cells that covers the training block
            training_block0 <- training_block %>%
              as("Spatial") %>%
              raster::rasterize(as(bias_file,"Raster"), getCover=TRUE)
            training_block0[training_block0==0] <- NA

            bias_file %<>%
              as("Raster") %>%
              raster::mask(mask=training_block0)
          }
          # save bias file raster in temporary folder
          outfile <-file.path(raster::tmpDir(),paste0(species_name,"_biasfile_",j,".asc"))
          raster::writeRaster(bias_file, outfile, format="ascii", overwrite=TRUE)
          # set the new biasfile file path
          maxent_settings[["biasfile"]] <- outfile
          maxent_settings[["biastype"]] <- 3
          # add sampling probs to test data
          extrct <- raster::extract(raster::raster(biasfile),testing_loc[,2:3],df=TRUE) %>%
            subset(select=2) %>%
            setNames(paste0(species_name,"_biasfile_",j))
          testing_loc <- testing_loc %>%
            dplyr::bind_cols(extrct)
        }else{
          #saveRDS(bias_file,paste0("E:/maxent_test/bias_file",j,".rds"))
          # remove biasfile option(s)
          maxent_settings<- maxent_settings[!names(maxent_settings)=="biasfile"]
        }
      }else{
        # remove biasfile option(s)
        maxent_settings<- maxent_settings[!names(maxent_settings)=="biasfile"]
      }
    }

    # save samples with data
    testingsampfile = file.path(dirs_to_create[j],'sptest.csv')
    write.csv(testing_loc, testingsampfile, row.names=FALSE)

    #-----------------------------
    #= c. run the cross-validation
    #-----------------------------
    outputdir_path = dirs_to_create[j]

    if(varying_parameter){
      # create new directories
      param_dn = paste(varying_parameter_name, maxent_settings[[varying_parameter_name]],sep="_")
      outputdir_path = file.path(outputdir_path, param_dn)
      outdirs_created = unlist(sapply(outputdir_path, function(x) if(!dir.exists(x)) dir.create(x)))
      if(!all(outdirs_created)) stop("Unable to create directories", outputdir_path[!outdirs_created])
    }

    # set the maximum number of background samples
    maxbackground <- min(sapply(training_bg, function(x) sum(!is.na(x))))

    if(length(maxent_settings)>0L){
      if("maximumbackground" %in% names(maxent_settings))
        maxent_settings[["maximumbackground"]] <- min(maxent_settings[["maximumbackground"]], maxbackground)
      else
        maxent_settings[["maximumbackground"]] <- min(10000, maxbackground)
    }else{
      maxent_settings <- list(maximumbackground = min(10000, maxbackground))
    }

    cv.model.args <-Reduce(append,list(

                                  list(sp_points=trainingsampfile,

                                  env_layers=dn,

                                  outputdir=outputdir_path[1]),

                                  list(testsamplesfile=testingsampfile),

                                  replace(maxent_settings,names(maxent_settings)==varying_parameter_name,maxent_settings[[varying_parameter_name]][1]),

                                  list(threads=ncores, pictures=FALSE, warnings = FALSE)))

    # train the model
    model_trained <- try(do.call(maxent, cv.model.args),silent=TRUE)

   if(inherits(model_trained,"try-error")){
     cat(model_trained)
     stop("Maxent model failed. Please check maxent.log file for more informations.")
   }

    # get the outputs
    model_output <- get_maxent_output(outputdir_path[1], eval_metrics=eval_metrics, index=j)

    if("ic" %in% eval_metrics){
      if(inherits(training_bg,"stars")){
        rst <- raster::stack(lapply(1:length(training_bg), function(l) as(training_bg[l], "Raster")))
        names(rst) <- names(training_bg)
      }else{
        rst <- training_bg
      }
      occ <- training_loc[,-1]
      lbds <- list.files(outputdir_path[1], pattern="\\.lambdas",full.names=TRUE)
      if(length(lbds)==0L)
        warning("Unable to find lambdas file in directory ", outputdir_path[1])

      pred_raw <- try(rmaxent::project(lambdas=lbds, newdata=rst, quiet=TRUE)$prediction_raw, silent=TRUE)

      if(inherits(pred_raw,"try-error")){
        print(geterrmessage())
        warning("An error occurred when predicting maxent model located in ", outputdir_path[1])
      }

      IC <- try(rmaxent::ic(x=pred_raw, occ=occ, lambdas = lbds),silent=TRUE)

      if(inherits(IC,"try-error"))
        IC<- data.frame(k=rmaxent::n_features(lbds),ll=NA, AIC=NA, AICc=NA, BIC=NA)

      model_output %<>%
        dplyr::bind_cols(as.data.frame(subset(IC, select=c("k","ll","AIC","AICc","BIC"))))

    }

    if(varying_parameter){

      # get the previous command line
      old_command_line <- get_maxent_command_line(file.path(outputdir_path[1],"maxent.log"))

      for(i in 2:length(maxent_settings[[varying_parameter_name]])){

        # change output directory
        new_output_directory <- outputdir_path[i]
        new_command_line <- set_command_param(old_command_line, "outputdirectory", new_output_directory)

        # change parameter value
        new_command_line <- set_command_param(new_command_line, varying_parameter_name, maxent_settings[[varying_parameter_name]][i])

        # run the model
        mem <- if("memory_allocated" %in% maxent_settings) maxent_settings$memory_allocated else 512
        new_command_line <- set_command_arg(new_command_line, "density.MaxEnt", paste0('-mx',mem,'m',' -jar ',cv.model.args$path_to_maxent))
        system(new_command_line)

        if("ic" %in% eval_metrics){
          #rst <- raster::stack(lapply(1:length(training_bg), function(l) as(training_bg[l], "Raster")))
          lbds <- list.files(new_output_directory, pattern="\\.lambdas",full.names=TRUE)
          if(length(lbds)==0L)
            stop("Unable to find lambdas file in directory ", new_output_directory)

          pred_raw <- try(rmaxent::project(lambdas=lbds, newdata=rst, quiet=TRUE)$prediction_raw, silent=TRUE)

          if(inherits(pred_raw,"try-error")){
            print(geterrmessage())
            warning("An error occurred when predicting maxent model located in ", new_output_directory)
          }

          IC <- try(rmaxent::ic(x=pred_raw, occ=occ, lambdas = lbds),silent=TRUE)

          if(inherits(IC,"try-error")){
            print(geterrmessage())
            IC<- tryCatch(data.frame(k=rmaxent::n_features(lbds),ll=NA, AIC=NA, AICc=NA, BIC=NA), error=function(err) data.frame(k=0,ll=NA, AIC=NA, AICc=NA, BIC=NA))
          }

          model_output %<>%
            dplyr::bind_rows(
            get_maxent_output(new_output_directory, eval_metrics=eval_metrics, index=j) %>%
            dplyr::bind_cols(as.data.frame(subset(IC, select=c("k","ll","AIC","AICc","BIC"))))
            )

        }else{
          # get the output
          model_output %<>%
            dplyr::bind_rows(get_maxent_output(new_output_directory,eval_metrics=eval_metrics, index=j))
        }

      }
      # add varying parameter column
      model_output %<>%
        dplyr::bind_cols(data.frame(varying_param=paste(varying_parameter_name, maxent_settings[[varying_parameter_name]],sep="_")))
    }

    # remove temporary file if needed
    on.exit({
      tmp <- file.path(raster::tmpDir(), paste0(species_name,"_biasfile_",j,".asc"))
      if(file.exists(tmp))
        file.remove(tmp)
    })

    model_output

    }

  },error=function(err){

    foreach::foreach(j=1:k, .combine=rbind) %do% {

      #--------------------------
      #= a. get the training data
      #--------------------------

      # training occurrence records
      training_block <- bg_masks %>%
        subset(layer %in% comb[,j]) %>%
        sf::st_make_valid()

      is_training_loc <- suppressMessages(sf_loc_data %>%
                                            sf::st_intersects(training_block) %>%
                                            purrr::map_lgl(., function(x) length(x) > 0L))
      training_loc <- data.frame(species=rep(species_name, sum(is_training_loc))) %>%
        dplyr::bind_cols(as.data.frame(sf::st_coordinates(sf_loc_data[is_training_loc,]))) %>%
        na.omit()

      # save samples
      trainingsampfile = file.path(dirs_to_create[j],'sptrain.csv')
      write.csv(training_loc, trainingsampfile, row.names=FALSE)

      # environmental layers
      if(env_data_flag){

        env_bg0 <- stars::st_as_stars(env_dat)

        # collpase attributes
        if(length(stars::st_dimensions(env_bg0))<3)
          env_bg0 %<>%
          merge() %>%
          stars::st_set_dimensions(3, values = names(env_bg0)) %>%
          stars::st_set_dimensions(names = c("x", "y", "band"))

        env_bg <- env_bg0 %>%
          as("Raster") %>%
          setNames(stars::st_dimensions(env_bg0)$band$values) %>%
          raster::crop(y=as(training_block,"Spatial")) %>%
          stars::st_as_stars()

        training_block0 <- training_block %>%
          as("Spatial") %>%
          raster::rasterize(as(env_bg,"Raster"), getCover=TRUE)
        training_block0[training_block0==0] <- NA

        training_bg <-try({
          names(env_bg) <- strip_extension(names(env_bg))
          env_bg %>%
            as("Raster") %>%
            setNames(stars::st_dimensions(env_bg0)$band$values) %>%
            raster::mask(mask=training_block0) %>%
            stars::st_as_stars()
        })

        if(inherits(training_bg,"try-error"))
          stop("Unable to read environmental layers from data env_dat. Please make sure all the layers have the same dimensions.")
      }
      else if(env_dir_flag){
        env_bg0 <- try(env_dat %>%
                        read_layers(),silent=TRUE)

        if(inherits(env_bg0,"try-error")){
          cat(env_bg0)
          stop("Unable to read environmental raster layers.")
        }

        training_bg <-try({
          names(env_bg0) <- if(any(is.na(strip_extension(names(env_bg0)))))  names(env_bg0) else strip_extension(names(env_bg0))

          if(inherits(env_bg0,"stars")){

            # collpase attributes
            if(length(stars::st_dimensions(env_bg0))<3)
              env_bg0 %<>%
              merge() %>%
              stars::st_set_dimensions(3, values = names(env_bg0)) %>%
              stars::st_set_dimensions(names = c("x", "y", "band"))

            # set environmental layers extent to training block extent
            tmp <- suppressMessages(env_bg0 %>%
                                      as("Raster") %>%
                                      setNames(stars::st_dimensions(env_bg0)$band$values) %>%
                                      raster::crop(y=as(training_block,"Spatial")) %>%
                                      stars::st_as_stars())

            st_dim <- attr(tmp,"dimension")
            # ugly work-around when crop returns crap
            if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){

              env_bg <- as(env_bg0,"Raster") %>%
                setNames(stars::st_dimensions(env_bg0)$band$values) %>%
                raster::crop(y=as(training_block,"Spatial"))

              training_block0 <- training_block %>%
                as("Spatial") %>%
                raster::rasterize(env_bg, getCover=TRUE)
              training_block0[training_block0==0] <- NA

              env_bg %>%
                raster::mask(mask=training_block0) %>%
                stars::st_as_stars()
            }else{
              # create a mask with cells that covers the training block
              training_block0 <- training_block %>%
                as("Spatial") %>%
                raster::rasterize(as(tmp,"Raster"), getCover=TRUE)
              training_block0[training_block0==0] <- NA

              tmp %>%
                as("Raster") %>%
                setNames(stars::st_dimensions(env_bg0)$band$values) %>%
                raster::mask(mask=training_block0) %>%
                stars::st_as_stars()
            }

          }else{

            # set environmental layers extent to training block extent
            env_bg <- env_bg0 %>%
              raster::crop(y=as(training_block,"Spatial"))

            # create a mask with cells that covers the training block
            training_block0 <- training_block %>%
              as("Spatial") %>%
              raster::rasterize(env_bg, getCover=TRUE)
            training_block0[training_block0==0] <- NA

            env_bg %>%
              raster::mask(mask=training_block0) %>%
              raster::stack()
          }

        })

        if(inherits(training_bg,"try-error"))
          stop("Unable to read environmental layers from ",env_dat,". Please make sure all the layers have the same dimensions.")
      }
      dn <- file.path(dirs_to_create[j],'env')
      if(!dir.exists(dn)) dir.create(dn)

      if(inherits(training_bg,"stars"))
        training_bg %<>%
        split("band") %>%
        setNames(stars::st_dimensions(env_bg0)$band$values)

      # save layers
      n_layers <- if(inherits(training_bg,"RasterStack")) raster::nlayers(training_bg) else length(training_bg)
      for(l in 1:n_layers){
        training_bg_cvrt <- if(inherits(training_bg,"RasterStack"))  training_bg[[l]] else as(training_bg[l],"Raster")
        raster::writeRaster(training_bg_cvrt, file.path(dn, paste0(names(training_bg)[l],'.asc')), format='ascii', overwrite=TRUE)
      }


      #-------------------------
      #= b. get the testing data
      #-------------------------
      testing_loc <- data.frame(species=rep(species_name, sum(!is_training_loc))) %>%
        dplyr::bind_cols(as.data.frame(sf::st_coordinates(sf_loc_data[!is_training_loc,])))

      # add environmental variables
      if(inherits(env_bg0,"stars")){
        env_bg0 %<>%
          split("band") %>%
          setNames(stars::st_dimensions(env_bg0)$band$values)
        testing_loc <- testing_loc %>%
          dplyr::bind_cols(map_lfd(env_bg0,function(x) raster::extract(as(x,"Raster"),.[,-1]))) %>%
          na.omit()
      }else{
        testing_loc <- testing_loc %>%
          dplyr::bind_cols(raster::extract(env_bg0,.[,-1]) %>% as.data.frame(row.names = NULL)) %>%
          na.omit()
      }

      # check biasfile option
      if(!is.null(maxent_settings[["biasfile"]])){
        biasfile <- maxent_settings[["biasfile"]]
        if(file.exists(biasfile) && !dir.exists(biasfile)){
          bias_file <- try(biasfile %>%
                             raster::raster() %>%
                             raster::crop(y=as(training_block,"Spatial")) %>%
                             stars::st_as_stars(),silent=TRUE)

          if(!inherits(bias_file,"try-error")){
            st_dim <- attr(bias_file,"dimension")
            if(st_dim$x$from==st_dim$x$to | st_dim$y$from==st_dim$y$to){

              biasfile0 <- biasfile %>%
                raster::raster() %>%
                raster::crop(y=as(training_block,"Spatial"))

              # create a mask with cells that covers the training block
              training_block0 <- training_block %>%
                as("Spatial") %>%
                raster::rasterize(biasfile0, getCover=TRUE)
              training_block0[training_block0==0] <- NA

              bias_file <- biasfile0 %>%
                #raster::crop(y=as(training_block,"Spatial")) %>%
                raster::mask(mask=training_block0)
            }else{
              # create a mask with cells that covers the training block
              training_block0 <- training_block %>%
                as("Spatial") %>%
                raster::rasterize(as(bias_file,"Raster"), getCover=TRUE)
              training_block0[training_block0==0] <- NA

              bias_file %<>%
                as("Raster") %>%
                raster::mask(mask=training_block0)
            }
            # save bias file raster in temporary folder
            outfile <-file.path(raster::tmpDir(),paste0(species_name,"_biasfile_",j,".asc"))
            raster::writeRaster(bias_file, outfile, format="ascii", overwrite=TRUE)
            # set the new biasfile file path
            maxent_settings[["biasfile"]] <- outfile
            maxent_settings[["biastype"]] <- 3
            # add sampling probs to test data
            extrct <- raster::extract(raster::raster(biasfile),testing_loc[,2:3],df=TRUE) %>%
              subset(select=2) %>%
              setNames(paste0(species_name,"_biasfile_",j))
            testing_loc <- testing_loc %>%
              dplyr::bind_cols(extrct)
          }else{
            # remove biasfile option(s)
            maxent_settings<- maxent_settings[!names(maxent_settings)=="biasfile"]
          }
        }else{
          # remove biasfile option(s)
          maxent_settings<- maxent_settings[!names(maxent_settings)=="biasfile"]
        }
      }

      # save samples with data
      testingsampfile = file.path(dirs_to_create[j],'sptest.csv')
      write.csv(testing_loc, testingsampfile, row.names=FALSE)

      #-----------------------------
      #= c. run the cross-validation
      #-----------------------------
      outputdir_path = dirs_to_create[j]

      if(varying_parameter){
        # create new directories
        param_dn = paste(varying_parameter_name, maxent_settings[[varying_parameter_name]],sep="_")
        outputdir_path = file.path(outputdir_path, param_dn)
        outdirs_created = unlist(sapply(outputdir_path, function(x) if(!dir.exists(x)) dir.create(x)))
        if(!all(outdirs_created)) stop("Unable to create directories", outputdir_path[!outdirs_created])
      }

      # set the maximum number of background samples
      maxbackground <- min(sapply(training_bg, function(x) sum(!is.na(x))))

      if(length(maxent_settings)>0L){
        if("maximumbackground" %in% names(maxent_settings))
          maxent_settings[["maximumbackground"]] <- min(maxent_settings[["maximumbackground"]], maxbackground)
        else
          maxent_settings[["maximumbackground"]] <- min(10000, maxbackground)
      }else{
        maxent_settings <- list(maximumbackground = min(10000, maxbackground))
      }

      cv.model.args <-Reduce(append,list(

        list(sp_points=trainingsampfile,

             env_layers=dn,

             outputdir=outputdir_path[1]),

        list(testsamplesfile=testingsampfile),

        replace(maxent_settings,names(maxent_settings)==varying_parameter_name,maxent_settings[[varying_parameter_name]][1]),

        list(threads=ncores, pictures=FALSE, warnings = FALSE)))

      # train the model
      model_trained <- try(do.call(maxent, cv.model.args),silent=TRUE)

      if(inherits(model_trained,"try-error")){
        cat(model_trained)
        stop("Maxent model failed. Please check maxent.log file for more informations.")
      }

      # get the outputs
      model_output <- get_maxent_output(outputdir_path[1], eval_metrics=eval_metrics, index=j)

      if("ic" %in% eval_metrics){
        if(inherits(training_bg,"stars")){
          rst <- raster::stack(lapply(1:length(training_bg), function(l) as(training_bg[l], "Raster")))
          names(rst) <- names(training_bg)
        }else{
          rst <- training_bg
        }
        occ <- training_loc[,-1]
        lbds <- list.files(outputdir_path[1], pattern="\\.lambdas",full.names=TRUE)
        if(length(lbds)==0L)
          warning("Unable to find lambdas file in directory ", outputdir_path[1])

        pred_raw <- try(rmaxent::project(lambdas=lbds, newdata=rst, quiet=TRUE)$prediction_raw, silent=TRUE)

        if(inherits(pred_raw,"try-error")){
          print(geterrmessage())
          warning("An error occurred when predicting maxent model located in ", outputdir_path[1])
        }

        IC <- try(rmaxent::ic(x=pred_raw, occ=occ, lambdas = lbds),silent=TRUE)

        if(inherits(IC,"try-error"))
          IC<- data.frame(k=rmaxent::n_features(lbds),ll=NA, AIC=NA, AICc=NA, BIC=NA)

        model_output %<>%
          dplyr::bind_cols(as.data.frame(subset(IC, select=c("k","ll","AIC","AICc","BIC"))))

      }

      if(varying_parameter){

        # get the previous command line
        old_command_line <- get_maxent_command_line(file.path(outputdir_path[1],"maxent.log"))

        for(i in 2:length(maxent_settings[[varying_parameter_name]])){

          # change output directory
          new_output_directory <- outputdir_path[i]
          new_command_line <- set_command_param(old_command_line, "outputdirectory", new_output_directory)

          # change parameter value
          new_command_line <- set_command_param(new_command_line, varying_parameter_name, maxent_settings[[varying_parameter_name]][i])

          # run the model
          mem <- if("memory_allocated" %in% maxent_settings) maxent_settings$memory_allocated else 512
          new_command_line <- set_command_arg(new_command_line, "density.MaxEnt", paste0('-mx',mem,'m',' -jar ',cv.model.args$path_to_maxent))
          system(new_command_line)


          if("ic" %in% eval_metrics){
            lbds <- list.files(new_output_directory, pattern="\\.lambdas",full.names=TRUE)
            if(length(lbds)==0L)
              warning("Unable to find lambdas file in directory ", new_output_directory)

            pred_raw <- try(rmaxent::project(lambdas=lbds, newdata=rst, quiet=TRUE)$prediction_raw, silent=TRUE)

            if(inherits(pred_raw,"try-error")){
              print(geterrmessage())
              warning("An error occurred when predicting maxent model located in ", new_output_directory)
            }

            IC <- try(rmaxent::ic(x=pred_raw, occ=occ, lambdas = lbds),silent=TRUE)

            if(inherits(IC,"try-error")){
              print(geterrmessage())
              IC<- tryCatch(data.frame(k=rmaxent::n_features(lbds),ll=NA, AIC=NA, AICc=NA, BIC=NA), error=function(err) data.frame(k=0,ll=NA, AIC=NA, AICc=NA, BIC=NA))
            }

            model_output %<>%
              dplyr::bind_rows(
                get_maxent_output(new_output_directory, eval_metrics=eval_metrics, index=j) %>%
                  dplyr::bind_cols(as.data.frame(subset(IC, select=c("k","ll","AIC","AICc","BIC"))))
              )

          }else{
            # get the output
            model_output %<>%
              dplyr::bind_rows(get_maxent_output(new_output_directory,eval_metrics=eval_metrics, index=j))
          }

        }
        # add varying parameter column
        model_output %<>%
          dplyr::bind_cols(data.frame(varying_param=paste(varying_parameter_name, maxent_settings[[varying_parameter_name]],sep="_")))
      }

      # remove temporary file if needed
      on.exit({
        tmp <- file.path(raster::tmpDir(), paste0(species_name,"_biasfile_",j,".asc"))
        if(file.exists(tmp))
          file.remove(tmp)
      })

      model_output

    }

  },finally = {
    if(do.parallel)
      doParallel::stopImplicitCluster()
    unregister()
  })

  return(cv_output)

}


#' utility function that performs a leave one out cross validation for geographic models
#' @export
loocv_geo <- function(pr, bg, algorithm="idw", tr=seq(from=1.E-3, to=0.99, by=1.E-2)){

  if(missing(pr))
    stop("Presence data locations are missing.")

  cv_eval_metrics <- vector(mode = "list", length = nrow(pr))

  for(i in 1:nrow(pr)){
    #======================
    #= 1. Train the model #
    #======================
    model <- switch(algorithm,
                      idw = try(dismo::geoIDW(p=pr[-i,], a=bg), silent=TRUE),
                      dist = try(dismo::geoDist(p=pr[-i,], lonlat=TRUE), silent=TRUE)
    )

    if(inherits(model,"try-error")){
      cat(model)
      stop("Unable to train the model")
    }

    #=========================
    #= 2. Evaluate the model #
    #=========================
    eval_geo_model <- dismo::evaluate(model=model, p=cbind(X=pr[i,1],Y=pr[i,2]), a=bg, tr=tr)

    cv_eval_metrics[[i]] <- data.frame(AUC=eval_geo_model@auc,
                                       TSS=max(eval_geo_model@TPR + eval_geo_model@TNR - 1),
                                       OR10=calc_omission_rate(model, occ_train=pr[-i,], occ_test=data.frame(X=pr[i,1],Y=pr[i,2]), threshold = 10))
  }

  colMeans(do.call(rbind,cv_eval_metrics))
}


#' utility function that calculates the omission rate from a model for different threshold values
#' @export
calc_omission_rate <- function(model,  occ_train, occ_test, threshold = 10) {

  suit_pred_train <- dismo::predict(model, occ_train)
  suit_pred_test <- dismo::predict(model, occ_test)

  om_rate <- vector("numeric", length = length(threshold))

  for (i in 1:length(threshold)) {
    val <- ceiling(length(occ_train[, 1]) * threshold[i] / 100) + 1

    quantval_or <- sort(suit_pred_train)[val]

    om_rate[i] <- as.numeric(length(suit_pred_test[suit_pred_test < quantval_or]) / length(suit_pred_test))
  }
  return(om_rate)
}


#' utility function: selects the best model based off a set of evaluation metrics.
#' @param model_output A dataframe object returned by the function \code{block_cv_maxent}.
#' @param eval_metrics A vector of evaluation metric names. Evaluation metrics available are: `auc`, `omission_rate`,`tss`, `ic`.
#' @param tolerance A numeric integer specifying the number of decimals to account for when comparing performance metrics. Default is 3 digits.
#' @export
get_best_maxent_model <- function(model_output, eval_metrics, tolerance=3){

  metric <- sortDirection <- c()
  for(j in eval_metrics){

    switch(j,

           ic = {
             metric <- append(metric,c("AICc","AIC","BIC"))
             sortDirection <- append(sortDirection,c(FALSE,FALSE,FALSE))
           },

           auc = {
             metric <- append(metric,"AUC_test")
             sortDirection <- append(sortDirection,TRUE)
           },

           omission_rate = {
             metric <- append(metric,"OR10_test")
             sortDirection <- append(sortDirection,FALSE)
            },

           tss = {
             metric <- append(metric,"TSS_test")
             sortDirection <- append(sortDirection,TRUE)
           }
    )
  }

  #-------------------------------------------------
  # average model performance across geographic bins
  #-------------------------------------------------
  if("varying_param" %in% names(model_output)){
    avg_metric <- model_output %>%
      dplyr::filter(k>0) %>%
      dplyr::select(c("varying_param",all_of(metric))) %>%
      dplyr::group_by(varying_param) %>%
      dplyr::summarise_at(metric, list(mean), na.rm=TRUE) %>%
      as.data.frame(row.names=NULL) %>%
      dplyr::mutate_at(.vars = metric,
                       .funs = ~round(.,tolerance))
    #-------------------------------------------------
    # sort from best to worst model
    #-------------------------------------------------
    sortOrder <- order(match(metric, c("AICc","AIC","OR10_test","AUC_test","TSS_test")))
    sortExpr <- sprintf("with(avg_metric, avg_metric[order(%s, decreasing=c(%s)), ])", paste(metric[sortOrder],collapse=", "), paste(sortDirection[sortOrder],collapse=", "))
    avg_metric <- eval(parse(text=sortExpr))

    # refine the selection if the two best models have similar performance in terms of AICc (i.e. difference in AICc <= 2)
    if("ic" %in% eval_metrics){

      first_two <- head(avg_metric$AICc, 2)

      if(length(first_two)>1){

        if(!any(is.na(first_two))){

          # for models performing similarly
          if(abs(diff(first_two))<=2){
            sortOrder <- order(match(metric, c("OR10_test","AUC_test","TSS_test"))) # drop off AIC criteria
            sortExpr <- sprintf("with(head(avg_metric,2), head(avg_metric,2)[order(%s, decreasing=c(%s)), ])", paste(metric[sortOrder],collapse=", "), paste(sortDirection[sortOrder],collapse=", "))
            avg_metric <- eval(parse(text=paste0("rbind(",sortExpr,",avg_metric[-c(1,2),])")))
          }

        }

      }

    }

  }else{
    avg_metric <- model_output %>%
      dplyr::filter(k>0) %>%
      dplyr::select(metric) %>%
      colMeans(na.rm=TRUE) %>%
      round(tolerance) %>%
      as.data.frame.list()
  }
  return(list(avg_performances = avg_metric, best_avg_performance=avg_metric[1,]))

}
