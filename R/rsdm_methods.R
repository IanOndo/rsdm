## -------------------------------------------------------------------------- #
## 0. Generic Functions definition ------------------------------------------
## -------------------------------------------------------------------------- #
NULL

setGeneric(name="exdet",
           def= function(trg,ref,...)
             standardGeneric("exdet") )

setGeneric(name="micdet",
           def= function(trg,ref,...)
             standardGeneric("micdet") )

## extrapolation detection for matrix objects ------------------------------------------------
##'
##' @rdname exdetectors
##' @export
##'
setMethod("exdet",
          signature(trg="matrix", ref="matrix"),
          definition = function(trg, ref, ...){
              return(arma_exdet(trg, ref, ...))
          })

## extrapolation detection for RasterBrick objects ------------------------------------------------
##'
##' @rdname exdetectors
##' @export
##'
setMethod("exdet",
          signature(trg="RasterBrick", ref="RasterBrick"),
          definition = function(trg, ref, ...){
            return(arma_exdet_raster(trg, ref, ...))
          })

## detection of most influential variable for matrix objects ------------------------------------------------
##'
##' @rdname exdetectors
##' @export
##'
setMethod("micdet",
          signature(trg="matrix", ref="matrix"),
          definition = function(trg, ref, ...){
            return(arma_micdet(trg, ref, ...))
          })

## detection of most influential variable for RasterBrick objects ------------------------------------------------
##'
##' @rdname exdetectors
##' @export
##'
setMethod("micdet",
          signature(trg="RasterBrick", ref="RasterBrick"),
          definition = function(trg, ref, ...){
            return(arma_micdet_raster(trg, ref, ...))
          })
