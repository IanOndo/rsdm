

.onUnload <- function (libpath) {
  library.dynam.unload("rsdm", libpath)
}
