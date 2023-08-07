#'                                 Query Species Occurrences
#'
#' Queries species occurrences from different databases online and offline
#'
#' @param species_name Either a vector of species names or a csv file where species names are stored in the first column.
#' @param data_sources A (vector of) character string of the name of the database to search for occurrence data. Must be either `gbif` or `bien` or both for now.
#' @param download_dir A directory where to save species occurrence records data downloaded from online databases.
#' @param output_dir A directory where to save species occurrence records data after formatting.
#' @param mc_core A numeric integer specifying the number of cores to be used for parallel computing. Should be >1 for parallel processing.
#' @return None
#' @export
queryOcc <- function(species_name	= NULL,
                     data_sources = c('gbif','bien'),
                     download_dir = "",
                     output_dir   = "",
                     run_name			= 'Test',
                     mc_cores			= NULL)
{
  #--------------------
  #= 1. Check inputs
  #--------------------
  if(length(species_name)==0L || nchar(species_name)==0L)
    stop("Please provide a vector of species names or a file containing species names.")

  file_flag= tryCatch(file.exists(species_name), error=function(err) FALSE) && !tryCatch(dir.exists(species_name), error=function(err) FALSE)
  vector_flag = !file_flag & inherits(mode(species_name),"character")

  if(file_flag && length(readLines(species_name))==0L)
    stop(paste("The file provided:",species_name," is empty!"))

  if(length(data_sources)==0L)
      stop("Please provide the name of one or several data sources")
  if(sum(data_sources %in% c('gbif','bien','biotime','rainbio', 'genesys', 'spLink','cwr_gbif'))==0L)
    stop("Please select one or several data sources among: 'gbif','bien','biotime','rainbio', 'genesys', 'spLink','cwr_gbif' ")

  if(any(!data_sources %in% c('gbif','bien','biotime','rainbio', 'genesys', 'spLink','cwr_gbif')))
    stop(paste(data_sources[which(!data_sources %in% c('gbif','bien','biotime','rainbio', 'genesys', 'spLink', 'cwr_gbif'))]," is/are not valid data source(s)."))

  if( 'gbif' %in% data_sources & !dir.exists(download_dir)){
    warning("Argument 'download_dir' is missing, and will be set to the current working directory");flush.console()
    download_dir = getwd()
  }
  # make sure that the output directory provided exists
  if(!dir.exists(output_dir)){
    warning(paste("Output directory:", output_dir,"does not exist or could not be found, and will be set to the current working directory ./queryOcc"));flush.console()
    output_dir=file.path(getwd(), 'queryOcc')
    dir.create(output_dir, recursive=TRUE)
  }

  if(!is.null(mc_cores) & !is.numeric(mc_cores))
    stop("Argument 'mc_cores' must be a numeric integer")

  if(is.null(mc_cores)) mc_cores = parallel::detectCores()-1

  #---------------------
  #= 2. Set parameters
  #---------------------
  if(vector_flag){
    file_tmp <- tempfile("queryOcc_species_names",fileext=".csv")
    write.csv(data.frame(Species=species_name), file_tmp, row.names=FALSE)
    species_name <- file_tmp
  }

  params = c(

    'species_name' 	= species_name,		# path to species csv files directory

    'data_source'   = if(length(data_sources)>1) paste(data_sources,collapse=",") else data_sources,

    'download_dir' 	= if(!dir.exists(download_dir)) NULL else download_dir,

    'output_dir' 	  = output_dir,

    'run_name'			= run_name,

    'mc_cores' 			= round(mc_cores)  # number of cores to use
  )

  #====================================================================
  #== 3. Run species search
  #====================================================================
  commandArgs <- function(...) paste0(names(params),'=',params)
  assign('commandArgs',commandArgs,envir=.GlobalEnv)
  source(system.file("extdata/workflows/occQuery.R", package='rsdm'))
}
