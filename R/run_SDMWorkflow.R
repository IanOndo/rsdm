#' Running SDM Workflow
#'
#' Execute a modelling pipeline R script
#'
#' @param base_directory A character string specifying the path to the directory where to store the outputs of the workflow
#' @param species_directory A character string specifying  the path to the directory where are stored the species occurrence records.
#' @param env_directory A character string specifying the path to the directory where are stored the environmental layers.
#' @param force_sorting A logical. Should the sorting of occurrence records be repeated ? Only useful if `species_directory` contains new or updated occurrence data. Default is `FALSE`.
#' @param workflow A character string specifying the path to the modelling script to run
#' @return None
#' @export
#' @examples
#' run_SDMWorkflow()
run_SDMWorkflow <-function(base_directory,
                           species_directory,
                           env_directory,
                           force_sorting = FALSE,
                           workflow = system.file("extdata","workflows","1UsefulPlants_Workflow.R",package="UsefulPlants"),
                           run_name = "Test",
                           mc_cores = NULL){
  if(workflow == ""){
    stop("Couldn't find the modelling workflow script", call.=FALSE)
  }

  #=================
  #= 1. Check inputs
  #=================
  if(missing(base_directory)){
    warning("Base directory is missing, and will be set to the current working directory");flush.console()
    base_directory = getwd()
  }
  if(missing(species_directory))
    stop("Species directory is missing.")

  if(missing(env_directory))
    env_directory=NULL

  if(!is.null(mc_cores) & !is.numeric(mc_cores))
    stop("Argument 'mc_cores' must be a numeric integer")

  if(is.null(mc_cores)) mc_cores = parallel::detectCores()-1 else min(mc_cores, parallel::detectCores())

  #===================
  #= 2. Set parameters
  #===================
  params = c(
    'base_directory' 			= base_directory,

    'species_directory' 	= species_directory,		# path to species csv files directory

    'env_directory' 			= env_directory,			# path to environmental predictors directory

    'force_sorting' = force_sorting,

    'run_name' = run_name,

    'mc_cores' 			= round(mc_cores)  # number of cores to use
  )

  #==================
  #== 3. Run workflow
  #==================
  commandArgs <- function(...) paste0(names(params),'=',params)
  source(workflow, local=TRUE)
}
