#' Parallel version of Geometric binary predicates from the package `sf`
#'
#' Compute geometric binary predicates of intersection between 2 simple features geometry sets in parallel
#'
#' @param x An object of class `sf`, `sfc` or `sfg`
#' @param y An object of class `sf`, `sfc` or `sfg`
#' @param nchunks A numeric value specifying the number of subgeometries the x object should be split into.
#' @param ncores A numeric integer specifying the number of cpu cores to use for parallel processing.
#' @return A sparse list, with list element i an integer vector with all indices j for which intersection(x[i],y[i]) is `TRUE` if they intersect, hence integer(0) if none of them is `TRUE`.
#' @export
st_par_intersects <- function(x, y, nchunks=10, ncores=getOption("mc.cores")){

  if(missing(x))
    stop("x is missing.")

  if(missing(y))
    stop("y is missing.")

  if(!inherits(x, c("sf","sfc","sfg")))
    stop("x must be an sf or sfc object.")

  if(!inherits(y, c("sf","sfc","sfg")))
    stop("y must be an sf or sfc object.")

  #--------------------------------------------------
  #= 1. Split the geometric features in nchunks
  #--------------------------------------------------
  ngeom <- length(sf::st_geometry(x))
  iterations <- chunk2(1:ngeom, nchunks)

  #--------------------------------------------------
  #= 2. Apply st_intersects to each chunk in parallel
  #--------------------------------------------------
  on.exit({
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()}
  )
  ncores = if(is.null(ncores)) parallel::detectCores()-1 else min( max(1, ncores), parallel::detectCores()) # number of cores to use
  max_Rprocs = max(0, 15 - getNumRprocess())
  ncores = min(ncores, max(max_Rprocs, 1))
  # update maximum allowed number of additional R processes
  #options(mc.cores=ncores)

  if(ncores < 2)
    foreach::registerDoSEQ()
  else
    doParallel::registerDoParallel(ncores)

  out <- vector(length=length(iterations))

  out <- tryCatch({
    foreach::foreach(j=1:length(iterations), .packages="sf", .combine='c') %dopar% {
      sf::st_intersects(x[iterations[[j]], ], y)
    }
  },error=function(err){
    foreach::foreach(j=1:length(iterations), .combine='c') %do% {
      sf::st_intersects(x[iterations[[j]]], y)
    }
  }, finally = {
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()
  })

  out
}

#' Compute geometric binary predicates of intersection between 2 simple features geometry sets in parallel
#' then return the intersecting features
#'
#' @export
st_par_intersection <- function(x, y, nchunks=10, ncores=getOption("mc.cores")){

  is_intersected <- st_par_intersects(x, y, nchunks=nchunks, ncores=ncores) %>%
    purrr::map_lgl(function(x) length(x) > 0)

  x[is_intersected, ]
}


#' Parallelized version of geometric binary predicates function st_overlaps
#' @export
st_par_overlaps <- function(x, y, nchunks=10, ncores=getOption("mc.cores")){

  if(missing(x))
    stop("x is missing.")

  if(missing(y))
    stop("y is missing.")

  if(!inherits(x, c("sf","sfc","sfg")))
    stop("x must be an sf or sfc object.")

  if(!inherits(y, c("sf","sfc","sfg")))
    stop("y must be an sf or sfc object.")

  #--------------------------------------------------
  #= 1. Split the geometric features in nchunks
  #--------------------------------------------------
  ngeom <- length(sf::st_geometry(x))
  iterations <- chunk2(1:ngeom, nchunks)

  #--------------------------------------------------
  #= 2. Apply st_overlaps to each chunk in parallel
  #--------------------------------------------------
  on.exit({
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()}
  )
  ncores = if(is.null(ncores)) parallel::detectCores()-1 else min( min(1, ncores), parallel::detectCores()) # number of cores to use
  max_Rprocs = max(0, 15 - getNumRprocess())
  ncores = min(ncores, max(max_Rprocs, 1))
  # update maximum allowed number of additional R processes
  #options(mc.cores=ncores)

  if(ncores < 2)
    foreach::registerDoSEQ()
  else
    doParallel::registerDoParallel(ncores)

  out <- vector(length=length(iterations))

  out <- tryCatch({
    foreach::foreach(j=1:length(iterations), .packages="sf", .combine='c') %dopar% {
      sf::st_overlaps(x[iterations[[j]], ], y)
    }
  },error=function(err){
    foreach::foreach(j=1:length(iterations), .combine='c') %do% {
      sf::st_overlaps(x[iterations[[j]]], y)
    }
  }, finally = {
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()
  })

  out
}

#' Parallelized version of geometric binary predicates function st_covered_by
#' @export
st_par_covered_by <- function(x, y, nchunks=10, ncores=getOption("mc.cores")){

  if(missing(x))
    stop("x is missing.")

  if(missing(y))
    stop("y is missing.")

  if(!inherits(x, c("sf","sfc","sfg")))
    stop("x must be an sf or sfc object.")

  if(!inherits(y, c("sf","sfc","sfg")))
    stop("y must be an sf or sfc object.")

  #--------------------------------------------------
  #= 1. Split the geometric features in nchunks
  #--------------------------------------------------
  ngeom <- length(sf::st_geometry(x))
  iterations <- chunk2(1:ngeom, nchunks)

  #--------------------------------------------------
  #= 2. Apply st_covered_by to each chunk in parallel
  #--------------------------------------------------
  on.exit({
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()}
  )
  ncores = if(is.null(ncores)) parallel::detectCores()-1 else min( min(1, ncores), parallel::detectCores()) # number of cores to use
  max_Rprocs = max(0, 15 - getNumRprocess())
  ncores = min(ncores, max(max_Rprocs, 1))
  # update maximum allowed number of additional R processes
  #options(mc.cores=ncores)

  if(ncores < 2)
    foreach::registerDoSEQ()
  else
    doParallel::registerDoParallel(ncores)

  out <- vector(length=length(iterations))

  out <- tryCatch({
    foreach::foreach(j=1:length(iterations), .packages="sf", .combine='c') %dopar% {
      sf::st_covered_by(x[iterations[[j]], ], y)
    }
  },error=function(err){
    foreach::foreach(j=1:length(iterations), .combine='c') %do% {
      sf::st_covered_by(x[iterations[[j]]], y)
    }
  }, finally = {
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()
  })

  out
}

#' Parallelized version of geometric binary predicates function st_within
#' @export
st_par_within <- function(x, y, nchunks=10, ncores=getOption("mc.cores")){

  if(missing(x))
    stop("x is missing.")

  if(missing(y))
    stop("y is missing.")

  if(!inherits(x, c("sf","sfc","sfg")))
    stop("x must be an sf or sfc object.")

  if(!inherits(y, c("sf","sfc","sfg")))
    stop("y must be an sf or sfc object.")

  #--------------------------------------------------
  #= 1. Split the geometric features in nchunks
  #--------------------------------------------------
  ngeom <- length(sf::st_geometry(x))
  iterations <- chunk2(1:ngeom, nchunks)

  #--------------------------------------------------
  #= 2. Apply st_within to each chunk in parallel
  #--------------------------------------------------
  on.exit({
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()}
  )
  ncores = if(is.null(ncores)) parallel::detectCores()-1 else min(max(1, ncores), parallel::detectCores()) # number of cores to use
  max_Rprocs = max(0, 15 - getNumRprocess())
  ncores = min(ncores, max(max_Rprocs, 1))
  # update maximum allowed number of additional R processes
  #options(mc.cores=ncores)

  if(ncores < 2)
    foreach::registerDoSEQ()
  else
    doParallel::registerDoParallel(ncores)

  out <- vector(length=length(iterations))

  out <- tryCatch({
    foreach::foreach(j=1:length(iterations), .packages="sf", .combine='c') %dopar% {
      sf::st_within(x[iterations[[j]], ], y)
    }
  },error=function(err){
    foreach::foreach(j=1:length(iterations), .combine='c') %do% {
      sf::st_within(x[iterations[[j]]], y)
    }
  }, finally = {
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()
  })

  out
}

#' Parallelized version of the cross-join function from optiRum package
#' @export
CJ_par_dt <- function(X, Y, nThreads=getOption("mc.cores")){
  on.exit({
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()
  })
  if(!is.null(nThreads)){
    if(nThreads>1)
      doParallel::registerDoParallel(nThreads)
    else
      foreach::registerDoSEQ()
  }else{
    nThreads = min(nThreads, parallel::detectCores()) # number of cores to use
    doParallel::registerDoParallel(nThreads)
  }
  foreach::foreach(j=iterators::iter(X, by='row'), .packages=c("optiRum","data.table"), .combine='rbind') %dopar% {
    gc()
    optiRum::CJ.dt(X, Y)
  }

}
