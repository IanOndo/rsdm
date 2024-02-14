## If the OS is Windows, set mclapply to the hackish version. Otherwise, leave the definition alone.
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {rsdm::mclapply.hack},
                    Linux   = {parallel::mclapply},
                    Darwin  = {parallel::mclapply})


#====================================================================
#== 0. Get parameters form command-line
#====================================================================
args = commandArgs(trailingOnly=TRUE)
for(arg in args) {

  arg_name = stringr::str_extract(arg,".+?(?=\\=)")

  switch(arg_name,

         'base_directory' 			= {BASE_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'species_directory' 		= {SPECIES_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'env_directory' 			= {ENV_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'geo_algorithm'  = {GEO_ALGO = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'force_sorting' = {FORCE_SORTING = as.logical(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'run_name'    = {RUN_NAME = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'points_proj' 		= {POINTS_PROJ = paste0("+",gsub(".+?(?<=\\=\\+)","",arg,perl=TRUE))},

         'mc_cores' 			= {MC_CORES = as.integer(gsub(".+?(?<=\\=)","",arg,perl=TRUE))}
  )
}
cat(paste0('#',paste(rep('-',times=100),collapse="")))
cat('\n')
cat(paste0("> Run ","'",RUN_NAME,"'",' started on: ',Sys.time()));started.at=proc.time()
cat('\n')
cat('#====================================================================\n')
cat('#== 1. Setup directories\n')
cat('#====================================================================\n\n')
# < 3 cells => points directory
# 10 > # of cells >= 3 cells => geographic distance-based directory
# >= 10 cells => SDM directory

xyOnly = FALSE
if(exists("ENV_DIR")){
  rsdm::init_main_directories(base_directory=BASE_DIR, species_directory=SPECIES_DIR, env_directory=ENV_DIR, force_sorting=FORCE_SORTING, ncores=MC_CORES)
}else{
  xyOnly = TRUE
  rsdm::init_main_directories(base_directory=BASE_DIR, species_directory=SPECIES_DIR, force_sorting=FORCE_SORTING, ncores=MC_CORES)
}
cat('\n\n')

if(!xyOnly){
  cat('#====================================================================\n')
  cat('#== 2. Prepare environmental layers\n')
  cat('#====================================================================\n')
  cat(">...copying environmental variables to ",file.path(BASE_DIR,"envModels","maxent","modelsVariables"),"...")
  files.copied <- file.copy(from=list.files(ENV_DIR, pattern="\\.tif$",full.names=TRUE),
                            to=file.path(BASE_DIR,"envModels","maxent","modelsVariables", basename(list.files(ENV_DIR, pattern="\\.tif$"))),
                            overwrite=TRUE)
  if(any(!files.copied)){
    cat("FAILED\n")
    stop("The following files could not be copied: \n", paste(list.files(ENV_DIR, pattern="\\.tif$",full.names=TRUE)[!files.copied],collapse="\n"))
  }else{
    cat("OK\n\n")
  }
}

cat('#====================================================================\n')
cat('#== 3. Prepare species occurrence data\n')
cat('#====================================================================\n\n')

# List of species to run with MaxEnt algorithm
if(!xyOnly){
  speciesListMaxEnt=list.files(file.path(BASE_DIR,"envModels","maxent","occurrences"), full.names=T)
}else
  speciesListMaxEnt=character(0)

# List of species to run with geographic distance-based algorithm
speciesListGeoDist=list.files(file.path(BASE_DIR,"geoModels","geo_dist","occurrences"), full.names=T)

# List of species to run as single points
speciesListPt=list.files(file.path(BASE_DIR,"pointModels","point","occurrences"), full.names=T)

# List of species to run
speciesList = c(speciesListPt, speciesListGeoDist, speciesListMaxEnt)
speciesAlgo = basename(dirname(dirname(speciesList)))

# List of algorithms to run
if(length(speciesList)==0L)
  stop("Oops !! Something probably went wrong with the species sorting...")

cat(paste(length(speciesList),'species will be run:\n'))
cat(paste('*',length(speciesListMaxEnt),"species with 'MaxEnt'\n"))
cat(paste('*',length(speciesListGeoDist),"species with a 'Geographic distance-based' model\n"))
cat(paste('*',length(speciesListPt),"species will be modelled as 'Points'\n"))

cat('\n\n')

if(!xyOnly){
  cat('#====================================================================\n')
  cat('#== 4. Model settings\n')
  cat('#====================================================================\n')

  cat('>...Setting up parameters for the models...')

  maxent_settings <- list(
    path_to_maxent=system.file('java/maxent.jar',package="dismo"),
    visible=FALSE,
    writemess=FALSE,
    writebackgroundpredictions=TRUE,
    maximumbackground=50000,
    betamultiplier=c(1,2,6,10),
    biasfile= system.file("exdata","sdm_biasfile.asc",package="UsefulPlants"),
    prefixes=FALSE,
    threshold=FALSE,
    hinge=TRUE,
    outputformat='raw',
    outputgrids=FALSE
  )

  geographic_distance_algorithm = ifelse(!exists("GEO_ALGO"), "idw", GEO_ALGO)
  geographic_grid_resolution = 0.1666667

  cat('OK\n')
  cat('\n\n')
}else{
  cat('#====================================================================\n')
  cat('#== 4. Model settings\n')
  cat('#====================================================================\n')

  cat('>...Setting up parameters for the models...')

  geographic_distance_algorithm = ifelse(!exists("GEO_ALGO"), "idw", GEO_ALGO)
  geographic_grid_resolution = 0.1666667

  cat('OK\n')
  cat('\n\n')
}

cat('#====================================================================\n')
cat('#== 5. Run Models\n')
cat('#====================================================================\n')

cat('...Starting computations...\n\n')
started.at = Sys.time()
runs <-mclapply(speciesList, function(j){

  # Get data for just that species
  if(file.exists(j)){
    species.data <- tryCatch(data.table::fread(j, showProgress=FALSE),error=function(err) return(NULL))
    k = gsub("\\.csv","",basename(j))
  }else{
    stop("Unable to locate the file :", j)
  }
  # if an error occurred during the reading returns an NULL
  if (is.null(species.data) || ncol(species.data) < 2){
    if(verbose) warning(paste0("Cannot read file from species",k))
    dn <- file.path(BASE_DIR,"problems")
    if(!dir.exists(dn))
      dir.create(dn, recursive =TRUE)
    fn <-file.path(dn,'Cannot_read_file.csv')
    if(!file.exists(fn))
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    else if(!k %in% read.csv(fn)$Species)
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    return(NULL)
  }
  # if the species has no points
  if(nrow(species.data)< 1){
    if(verbose) warning(paste0("Species '",k,"' has no points"))
    dn <- file.path(BASE_DIR,"problems")
    if(!dir.exists(dn))
      dir.create(dn, recursive =TRUE)
    fn <-file.path(dn,'Has_no_points.csv')
    if(!file.exists(fn))
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    else if(!k %in% read.csv(fn)$Species)
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    return(NULL)
  }

  ###########
  # OPTIONS #
  ###########
  memory.limit(90000)	# Increase size memory allocated
  terraOptions(chunksize = 2e+05, maxmemory = 2e+07, todisk=FALSE)

  algo = speciesAlgo[match(j, speciesList)]

  out = try({

    switch(algo,

     maxent = rsdm::run_maxentModel(loc_dat = j,
                                              species_name = k,
                                              outputdir = file.path(BASE_DIR,"envModels","maxent"),
                                              do.map=TRUE,
                                              maxent_settings = maxent_settings,
                                              grid_res = max(0.008333335, geographic_grid_resolution/2),
                                              varying_parameter_name="betamultiplier",
                                              verbose=TRUE),

     geo_dist = rsdm::run_geoModel(loc_dat = j,
                                             species_name = k,
                                             algorithm = geographic_distance_algorithm,
                                             outputdir = file.path(BASE_DIR,"geoModels","geo_dist"),
                                             train_res = geographic_grid_resolution,
                                             project_res = geographic_grid_resolution,
                                             verbose = TRUE),

     point = rsdm::run_pointModel(loc_dat = j,
                                            species_name = k,
                                            outputdir = file.path(BASE_DIR,"pointModels","point"),
                                            project_res = geographic_grid_resolution,
                                            bbox_ext = 2,
                                            verbose = TRUE)
  )},silent=TRUE)

  if(inherits(out,'try-error')){
    message(out)
    gc()
    rmv <- terra::tmpFiles(orphan=TRUE,remove=TRUE)
    j
  }
  gc()
  rmv <- terra::tmpFiles(orphan=TRUE,remove=TRUE)
}, mc.cores=MC_CORES, mc.preschedule = FALSE)
finished.at = Sys.time()
time.elapsed <- finished.at - started.at
cat('...End of computations...\n\n')
cat(paste0('> End of Run ',"'",RUN_NAME,"'",' on: ',Sys.time()),'\n\n');finished.at=proc.time()
cat(paste0('> Running time: ',as.numeric(time.elapsed),attr(time.elapsed,'units'),'.\n'))
cat(paste0('#',paste(rep('-',times=100),collapse="")))
