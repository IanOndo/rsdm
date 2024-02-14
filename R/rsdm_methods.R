## -------------------------------------------------------------------------- #
## 0. Generic Functions definition ------------------------------------------
## -------------------------------------------------------------------------- #
NULL

setGeneric(name="nt1",
           def= function(trg,ref,...)
             standardGeneric("nt1") )

setGeneric(name="nt2",
           def= function(trg,ref,...)
             standardGeneric("nt2") )

setGeneric(name="exdet",
           def= function(trg,ref,...)
             standardGeneric("exdet") )

setGeneric(name="micdet",
           def= function(trg,ref,...)
             standardGeneric("micdet") )

## extrapolation detection of type 1 for matrix objects ------------------------------------------------
##'
##' @rdname nt1
##' @export
##'
setMethod("nt1",
          signature(trg="matrix", ref="matrix"),
          definition = function(trg, ref, ...){
            return(arma_calc_nt1(trg, ref, ...))
          })

## extrapolation detection of type 2 for matrix objects ------------------------------------------------
##'
##' @rdname nt2
##' @export
##'
setMethod("nt2",
          signature(trg="matrix", ref="matrix"),
          definition = function(trg, ref, ...){
            return(arma_calc_nt2(trg, ref, ...))
          })

## extrapolation detection for matrix objects ------------------------------------------------
##'
##' @rdname exdet
##' @export
##'
setMethod("exdet",
          signature(trg="matrix", ref="matrix"),
          definition = function(trg, ref, ...){
              return(arma_exdet(trg, ref, ...))
          })

## extrapolation detection for RasterBrick objects ------------------------------------------------
##'
##' @rdname exdet
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
