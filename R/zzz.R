
.onAttach <- function(libname, pkgname) {

 packageStartupMessage("beadarray versions >= 2.0.0 are substantial updates from beadarray 1.16.0 and earlier. Please see package vignette for details")
  
      #setHook(packageEvent("ggplot2", "attach"),
            
        #function(...) {
        #    detach("package:plyr")
            # setGeneric("summarize", function(object) standardGeneric("summarize"))
            # setMethod("summarize", signature(object="beadLevelData"), function(object, channelList=list(greenChannel), probeIDs=NULL, useSampleFac = FALSE, sampleFac= NULL, weightNames = "wts", removeUnMappedProbes = TRUE) beadarray::summarize(object, "Detection"))
        #})

}

