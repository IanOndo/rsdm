#' Sort species by number of occurrence records retrieve from occurrence data files.
#'
#' Assign species occurrence data to different folders based on minimum number of occurrence records
#'
#' @param base_directory A character string specifying the directory where to store the outputs of the workflow
#' @param species_directory A character string specifying the directory where are stored the species occurrence records.
#' @param env_directory A character string specifying the directory where are stored the environmental layers.
#' @param min_occ_envmodel A numeric integer specifying the minimum number of occurrence records required to assign a species to the environmental model directory
#' @param min_occ_geomodel A numeric integer specifying the minimum number of occurrence records required to assign a species to the geographic model directory
#' @param nchunk A numeric integer specifying the number of chunks in which the vector of species names must be divided. Must be > 1 to enable parallel processing.
#' @param ncores A numeric integer specifying the number of cores (or CPU's) to use in case of parallel processing. Default is `number of cores available` - 1.
#' @export
sortByAvailableRecords <- function(base_directory, species_directory, env_directory, algorithm, min_occ_envmodel=10, min_occ_geomodel=3, nchunk=20, tolerance=5, ncores=NULL, verbose=FALSE){

  cat("\n")
  if(verbose) cat(">...checking input directories...")
  # ensure all arguments are provided
  if(missing(base_directory)){
    cat("FAILED\n")
    stop("The base directory is missing.")
  }
  if(missing(species_directory)){
    cat("FAILED\n")
    stop("The directory with species occurrence data is missing.")
  }

  # ensure that the directories provided exist
  if(!dir.exists(base_directory)){
    cat("FAILED\n")
    stop("Base directory:", base_directory,"does not exist or could not be found.");flush.console()
  }

  if(verbose) cat("OK\n")
  if(verbose) cat("> ...ensure species directory has data...")
  # ensure species directory has data
  list.species <- list.files(species_directory, pattern="\\.csv$", full.names=TRUE, recursive=TRUE) # select species
  if(length(list.species)==0L)
    stop(paste0("\nUnable to find csv files in directory:",species_directory,". Please provide a directory with csv files."))
  empty.string <- nchar(list.species)==0L
  if(any(empty.string)){
    if(verbose)
      warning(paste("\nRemoving",sum(empty.string),"file(s) with no species names."))
    list.species <- list.species[!empty.string]
  }
  if(verbose) cat("OK\n")

  if(!is.numeric(nchunk)){
    stop("nchunk must be numeric")
  }
  if(nchunk<=0){
    stop("nchunk must be > 0")
  }

  ncores = if(is.null(ncores)) parallel::detectCores()-1 else min(ncores, parallel::detectCores()) # number of cores to use

  #if(verbose) cat("> ...split species files list in ",nchunk,"chunks and sort species by number of occurrence records...")

  #species_split_list <- chunk2(list.species, nchunk)

  if(length(list.species)>1)
    doParallel::registerDoParallel(ncores)
  else foreach::registerDoSEQ()

  assigned <- foreach::foreach(f=list.species, .packages=c("data.table","UsefulPlants","raster","dplyr"), .combine="c") %dopar% {

      dat <- try(data.table::fread(f), silent=TRUE)

      if(inherits(dat,"try-error"))
        return(-1)

      data.table::setDF(dat)

      #---------------------------------------------------------------------------
      #= 1. get coords (with a regular expression to look for longitude, latitude)
      #---------------------------------------------------------------------------
      x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
      y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
      coordHeaders <- c(x_lon,y_lat)

      #if the data file has bad headers
      if(all(!coordHeaders %in% names(dat))){
        dn <- file.path(base_directory,"problems/has_bad_headers")
        if(!dir.exists(dn))
          dir.create(dn, recursive =TRUE)
        file.copy(from=f,
                  to=file.path(dn, basename(f)))
        return(-1)
      }
      #---------------------------------------------------------------------------
      #= 2. remove duplicated coords
      #---------------------------------------------------------------------------
      xtmp  <- dat[, coordHeaders]
      repeated <- (1:nrow(xtmp)) >= which(coordHeaders[1] == xtmp[, coordHeaders[1]])
      if(any(repeated))
        xtmp <- xtmp[!repeated,]

      # make sure coordinates are numeric
      if(any(!sapply(xtmp[, coordHeaders], is.numeric)))
        xtmp[, coordHeaders][, !sapply(xtmp[, coordHeaders], is.numeric)] <-  sapply(xtmp[, coordHeaders][, !sapply(xtmp[, coordHeaders], is.numeric)],
                                                                                     function(f) if(is.factor(f)) as.numeric(levels(f))[f] else as.numeric(f))

      npts <- nrow(xtmp[!duplicated(round(xtmp, tolerance)),])

      #---------------------------------------------------------------------------
      #= 3. assign species to directories
      #---------------------------------------------------------------------------

      # if the data has no input
      has_no_input_data <- npts==0
      if(has_no_input_data){
        dn <- file.path(base_directory,"problems/has_no_input_data")
        if(!dir.exists(dn))
          dir.create(dn, recursive =TRUE)
        file.copy(from=f,
                  to=file.path(dn, basename(f)))
      }

      # if the species has enough points for modelling with environmental variables
      has_enough_data <- FALSE
      if(npts >= min_occ_envmodel){

        dirs_created <- init_model_directory(base_directory, type="environmental", algorithm="maxent")

        out_dn=file.path(base_directory,"envModels","maxent","occurrences",basename(f))

        if(!exists("env_directory")){
          warning("Cannot remove records with no environemental variable because the environmental layers directory is missing.")
        }
        env_layers <- try(raster::stack(list.files(env_directory,pattern="\\.tif$",full.names=TRUE)),silent=TRUE)

        if(inherits(env_layers,'try-error')){
          cat(env_layers)
          warning("Unable to read the environmental layers.
                  The occurrence data will not be write in environmental models folders.")
          # if(any(repeated))
          #   dat <- dat[!repeated,]
          #
          # if(any(!sapply(dat[, coordHeaders], is.numeric)))
          #   dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)] <-  sapply(dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)], function(f)as.numeric(levels(f))[f])
          #
          # write.csv(dat, out_dn, row.names=FALSE)
          #
          # has_enough_data <- TRUE # in order not to re-write data in geographic models folders

        }else{
          #= 3.a remove duplicated records and records with environmental values not available
          data_extract <- raster::extract(env_layers, xtmp, df=TRUE)[,-1]
          is_not_duplicated <-  data_extract %>%
            duplicated() %>% `!`

          is_complete <- data_extract %>%
            complete.cases()

          is_valid <- is_not_duplicated & is_complete

          if(sum(is_valid) >= min_occ_envmodel){

            if(any(repeated))
              dat <- dat[!repeated,]

            if(any(!sapply(dat[, coordHeaders], is.numeric)))
              dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)] <-  sapply(dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)],
                                                                                         function(f) if(is.factor(f)) as.numeric(levels(f))[f] else as.numeric(f))

            dat<- dat[is_valid,]

            write.csv(dat, out_dn, row.names=FALSE)

            # # copy environmental variables
            # files.copied <- file.copy(from=list.files(env_directory, pattern="\\.tif$"), to=file.path(directory,"envModels","modelsVariables", basename(list.files(env_directory, pattern="\\.tif$"))))
            #   if(any(!files.copied))
            #     stop("The following files could not be copied: ", list.files(env_directory, pattern="\\.tif$")[!files.copied])

            has_enough_data <- TRUE
          }

        }
      }

      # if the species does not have enough points for modelling based on environmental variables
      if( (npts >=min_occ_geomodel & npts < min_occ_envmodel) | (!has_enough_data & npts >=min_occ_geomodel) ){

        dirs_created <- init_model_directory(base_directory, type="geographic", algorithm="geo_dist")

        out_dn <- file.path(base_directory,"geoModels","geo_dist","occurrences",basename(f))

        if(!dirs_created){

          if(any(repeated))
            dat <- dat[!repeated,]

          if(any(!sapply(dat[, coordHeaders], is.numeric)))
            dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)] <-  sapply(dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)],
                                                                                       function(f) if(is.factor(f)) as.numeric(levels(f))[f] else as.numeric(f))

          dat<- dat[!duplicated(round(xtmp, tolerance)),]

          write.csv(dat, out_dn, row.names=FALSE)
        }
      }

      # if the species does not have enough points for modelling based on geographic distance
      if(npts < min_occ_geomodel){

        dirs_created <- init_model_directory(base_directory, type="point", algorithm="point")

        out_dn <- file.path(base_directory,"pointModels","point","occurrences",basename(f))

        if(!dirs_created){

          if(any(repeated))
            dat <- dat[!repeated,]

          if(any(!sapply(dat[, coordHeaders], is.numeric)))
            dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)] <-  sapply(dat[, coordHeaders][, !sapply(dat[, coordHeaders], is.numeric)],
                                                                                       function(f) if(is.factor(f)) as.numeric(levels(f))[f] else as.numeric(f))

          dat<- dat[!duplicated(round(xtmp, tolerance)),]

          write.csv(dat, out_dn, row.names=FALSE)
        }

      }

      gc()

      return(0)
  }

  if(length(list.species)>1){
    doParallel::stopImplicitCluster()
  }
  unregister()

  if(sum(assigned)){
    warning("The following file encountered an error:")
    sapply(paste0(list.species[assigned<0],"\n"), cat)
  }


}

#' helper function that split a vector x in n chunks of approximatively equal size
#' @export
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
