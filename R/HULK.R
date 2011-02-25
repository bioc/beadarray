HULK <- function(BLData, array = 1, neighbours = NULL, invasions = 20, useLocs = TRUE, weightName ="wts", transFun = logGreenChannelTransform) {

  
    an <- sectionNames(BLData)
    tmp = BLData[[array]]
    probeIDs = tmp[,"ProbeID"]
    data <- transFun(BLData, array)

    out = NULL	

    for(i in array) {
    #print something encouraging
        cat("Array",i,":\n")

        if(is.null(neighbours)) {
            cat("Calculating Neighbourhood\n")
            neighbours <- generateNeighbours(BLData, i, useLocs=useLocs)
        }

        out[[i]] = data - HULKResids(BLData, i, transFun, useLocs, neighbours, invasions, weightName = weightName)
    }
    out
}

HULKResids <- function(BLData, array, transFun = logGreenChannelTransform, useLocs = TRUE, neighbours = NULL, invasions = 20, weightName = "wts") {

	
    tmp = BLData[[array]]

    probeIDs = tmp[,"ProbeID"]
            
    if(weightName %in% colnames(tmp)){
        weights <- tmp[,"wts"]
    }
    else {
        weights = rep(1, length(probeIDs))
    }
            
    an <- sectionNames(BLData)

    probeIDs = tmp[,"ProbeID"]

    data <- transFun(BLData, array)

   ###Remove outliers first	

    outliers = illuminaOutlierMethod(data, probeIDs)

    	

    beadTypeMeans = lapply(split(data[-outliers], probeIDs[-outliers]), mean, na.rm=TRUE)

    residuals = data - unlist(beadTypeMeans)[as.character(probeIDs)]
    residuals[which((is.na(residuals)) | (weights == 0))] = 0
    
    if(is.null(neighbours)) {
        cat("Calculating Neighbourhood\n")
        neighbours <- generateNeighbours(BLData, array, useLocs = useLocs)
    }
    
    cat("HULKing\n")      
    output <- .C("HULK", as.double(residuals), as.integer(t(neighbours)), as.integer(nrow(neighbours)), as.integer(invasions), results = as.double(residuals))
    output$results
}
