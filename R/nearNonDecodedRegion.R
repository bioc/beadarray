givenrow<-function(x, UndecNN){
    NDclass<-NULL
    tobeadded<-x
    while(length(tobeadded)>0){
        newx<-tobeadded[1]
        NDclass<-c(NDclass,newx)
        tobeadded<-tobeadded[-1]
        temp<-UndecNN[newx,2:7]
        temp<-temp[!is.na(temp)]
        temp<-temp[!(temp %in% NDclass)]
        temp<-temp[!(temp %in% tobeadded)]
        tobeadded<-c(tobeadded,temp)
    }
    return(NDclass)
}


takeclusterexpand <- function(probelist, newNN){
        probelist<-probelist[!is.na(probelist)]
        probelist<-as.vector(newNN[probelist,])
        probelist<-probelist[!is.na(probelist)]
        probelist<-unique(probelist)
        probelist
}

neighboursMatrixForAll <- function(nsegs = 9, ncols = 397, nrows = 326) {
## calculate the neighbours matrix assuming 
## all points in the grid are filled
## default values are for Expression Arrays at the moment
    polaritydown<-T
    NN<-matrix(NA,nrow=nsegs*nrows*ncols, ncol=6)
    i<-0
    storepolaritydown<-polaritydown
    for(segment in 1:nsegs){
        #cat(segment,"\n")
        for(col in 1:ncols){
            for(bead in 1:nrows){
                i<-i+1

                temp<-NULL
                if(bead>1){temp<-c(temp,i-1)}
                if(bead<nrows){temp<-c(temp,i+1)}

                if(polaritydown){
                        if(col<ncols){
                                temp<-c(temp,i+nrows)           
                                if(bead>1){temp<-c(temp,i+nrows-1)}
                        }
                        if(col>1){
                                temp<-c(temp,i-nrows)
                                if(bead>1){temp<-c(temp,i-nrows-1)}
                        }
                }
                else{
                        if(col<ncols){
                                temp<-c(temp,i+nrows)   
                                if(bead<nrows){temp<-c(temp,i+nrows+1)}
                        }
                        if(col>1){
                                temp<-c(temp,i-nrows)           
                                if(bead<nrows){temp<-c(temp,i-nrows+1)}
                        }
                }
                NN[i,1:length(temp)]<-temp
                }#bead
            polaritydown<-(!polaritydown)
            }#col
        polaritydown<-storepolaritydown
        }#segment
    return(NN)
}

beadsNearNonDecoded <- function(BLData, array = 1) {
    
    sdfFileName <- file.path(BLData@sectionData$Targets$directory[1], list.files(as.character(BLData@sectionData$Targets$directory[1]), pattern = ".sdf")[1]);
    if(file.exists(sdfFileName)) {
        sdf <- simpleXMLparse(readLines(sdfFileName, warn = FALSE))
    }
    else {
        stop("sdf file cannot be located.\nAborting registration check\n");
    }
    
    nSegs <- as.integer(sdf$RegistrationParameters$SizeBlockY[[1]]);
    nRows <- as.integer(sdf$RegistrationParameters$SizeGridX[[1]]);
    nCols <- as.integer(sdf$RegistrationParameters$SizeGridY[[1]]);
    
    
    locs <- obtainLocs(fileName = BLData@sectionData$Targets$greenLocs[array], filePath = BLData@sectionData$Targets$directory[array]);
       
#     ## create a neighbours matrix including all of the non-decoded beads
#     segmentNeighbours <- .Call("neighboursFromLocs", as.integer(c(nRows, nCols)), PACKAGE = "beadarray")   
#     locsNeighbours <- matrix(nrow = nSegs * nRows * nCols, ncol = 7)
#     for(i in 1:nSegs) {
#         locsNeighbours[(1 + (i-1) * (nRows * nCols) ):(i * (nRows * nCols) ),] <- segmentNeighbours + ( (i-1) * (nRows * nCols) )
#     }
    
    NN <- neighboursMatrixForAll(nsegs = nSegs, nrows = nRows, ncols =  nCols)
    #NN <- locsNeighbours[,2:6]
    
    output <- NULL
    
    
    for(i in array) {
    
        cat("Array",i,":\n")
    
        ## combine the data with the loc file
        combed <- BeadDataPackR:::combineFiles(BLData[[i]][,c("ProbeID", "Grn", "GrnX", "GrnY")], locs)
        
        newNN<-cbind(1:nrow(combed),NN)
        revind<-cbind(1:nrow(combed),combed[,5])[order(combed[,5]),]      
            
        for(k in 1:7){
            newNN[!is.na(newNN[,k]),k]<-revind[newNN[!is.na(newNN[,k]),k],1]
        }
        
        newNN<-newNN[order(newNN[,1]),]
        
        limUD<-sum(combed[,1]==0)
        UndecNN<-newNN[1:limUD,]

        UndecNN[UndecNN>limUD]<-NA

        targs<-which(apply(!is.na(UndecNN[,2:7]),1,sum)==6)

        toassign<-rep(F,limUD)
        toassign[targs]<-T
        nclass<-0
        maxsize<-0

        while(sum(toassign)>0){
            temp<-givenrow(min(which(toassign)), UndecNN)
            maxsize<-max(maxsize,length(temp))
            toassign[temp]<-F
            nclass<-nclass+1
        }

        NDclasses<-matrix(NA,nrow=maxsize,ncol=nclass)

        toassign<-rep(F,limUD)
        toassign[targs]<-T
        nclass<-0
        maxsize<-0

        while(sum(toassign)>0){
            temp<-givenrow(min(which(toassign)), UndecNN)
            toassign[temp]<-F
            nclass<-nclass+1
            NDclasses[1:length(temp),nclass]<-temp
        }

        if((is.null(ncol(NDclasses))) || (ncol(NDclasses) == 0) ) {
            message("No clusters")
            return(NA);
        }

        minclust<-50
        NDclasses <- NDclasses[,apply(!is.na(NDclasses),2,sum)>=minclust]
        if(length(NDclasses > 0) && (!is.matrix(NDclasses)))
            NDclasses <- matrix(NDclasses, ncol = 1)

        #### Now we have our clusters
        if((is.null(ncol(NDclasses))) || (ncol(NDclasses) == 0) ) {
            message("No clusters")
            return(NA);
        }

        idx <- NULL
        for(j in 1:ncol(NDclasses)) {
            
            temp1 <- NDclasses[,j]
            temp1 <- temp1[!is.na(temp1)]
            temp2 <- temp1
            for(k in 1:5) {
                temp2 <- takeclusterexpand(temp2, newNN);
            }
            idx <- c(idx, temp2[!(temp2 %in% temp1)])
        }
        idx <- idx - limUD
        idx <- idx[-which(idx < 0)]
        output$wts[[i]] <- !1:nrow(BLData[[i]]) %in% idx;
    }
    return(output);
}
