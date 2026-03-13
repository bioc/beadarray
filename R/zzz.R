
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("Welcome to beadarray version ", packageDescription("beadarray", fields="Version"))
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.24")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))

 # packageStartupMessage("beadarray versions >= 2.0.0 are substantial updates from beadarray 1.16.0 and earlier. Please see package vignette for details")
  
      #setHook(packageEvent("ggplot2", "attach"),
            
        #function(...) {
        #    detach("package:plyr")
            # setGeneric("summarize", function(object) standardGeneric("summarize"))
            # setMethod("summarize", signature(object="beadLevelData"), function(object, channelList=list(greenChannel), probeIDs=NULL, useSampleFac = FALSE, sampleFac= NULL, weightNames = "wts", removeUnMappedProbes = TRUE) beadarray::summarize(object, "Detection"))
        #})

}

