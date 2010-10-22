       
setMethod("combine", c(x = "beadLevelData", y = "beadLevelData"), function(x, y) {
    
    beadLevelData <- new(Class = "beadLevelData");
    
    ## process experiment data
    beadLevelData@experimentData <- x@experimentData;
    
    ## process the section data
    for(i in union(names(x@sectionData), names(y@sectionData))) {
        
        if(i == "SampleGroup") {
            if(any(y@sectionData$SampleGroup[,1] %in% x@sectionData$SampleGroup[,1])) {
                ySampleGroup <- paste(y@sectionData$SampleGroup[,1], 1, sep = ".")
            }
            else { 
                ySampleGroup <- y@sectionData$SampleGroup[,1] 
            }
            beadLevelData <- insertSectionData(beadLevelData, what = i, 
                                   data = c(as.character(x@sectionData$SampleGroup[,1]), as.character(ySampleGroup)) )
        }
        else {
            if( identical(colnames(x@sectionData[[i]]), colnames(y@sectionData[[i]])) ) {
                beadLevelData <- insertSectionData(beadLevelData, what = i, 
                                   data = rbind(x@sectionData[[i]], y@sectionData[[i]]) )
            }
            else {
                beadLevelData <- insertSectionData(beadLevelData, what = i, 
                                   data = merge(x@sectionData[[i]], y@sectionData[[i]], all = TRUE, sort = FALSE) )
            }
        }
    }
    
    ## process the bead data
    for(i in 1:dim(x)[1]) {
        sn = sectionNames(x)[i];
        tmp <- x[[i]]; ## assigning to a tmp variable make the function much faster
        for(j in 1:ncol(tmp)) {
            beadLevelData <- insertBeadData(beadLevelData, array = i, what = colnames(tmp)[j], data = tmp[,j])
        }
    }
    for(i in 1:dim(y)[1]) {
        sn = sectionNames(y)[i];
        tmp <- y[[i]]; ## assigning to a tmp variable make the function much faster
        for(j in 1:ncol(tmp)) {
            beadLevelData <- insertBeadData(beadLevelData, array = i + dim(x)[1], what = colnames(tmp)[j], data = tmp[,j])
        }
    }
    
    return(beadLevelData);
    } )
    