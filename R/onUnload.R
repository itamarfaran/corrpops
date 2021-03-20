# https://stackoverflow.com/questions/27076732/dynamic-library-not-loading-in-r-binary-package-build
.onUnload <- function (libpath){
  library.dynam.unload("corrpops", libpath)
}
