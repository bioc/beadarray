imageplot <- function(BLData, array = 1, transFun = logGreenChannelTransform, squareSize = NULL, useLocs = TRUE, horizontal = TRUE, low = NULL, high = NULL, ncolors = 100, zlim=NULL, legend=TRUE,...) {
    
    data = transFun(BLData, array = array)
    
    if(class(BLData)[1] %in% c("RGList", "MAList")) 
        stop("\nIt appears you are trying to use the imageplot() function on a Limma object, but imageplot() is currently masked by beadarray\n\nIf you wish to use the Limma function, you can either call it directly using:\n\t\"limma::imageplot()\"\nor detach the beadarray package using:\n\t\"detach(package:beadarray)\"\n")

    ## see if we can find the .locs file and use that
    locsFileName <- file.path(BLData@sectionData$Targets$directory[array], paste(BLData@sectionData$Targets$sectionName[array], "_Grn.locs", sep = ""))
       
    ## now see if this is a Matrix or BeadChip.  Locs should only be used on chips;
    if( (!is.null(BLData@experimentData$platformClass)) && (useLocs) ) {
        useLocs <- ifelse(grepl("Matrix", BLData@experimentData$platformClass), FALSE, TRUE);
    }
    else { ## don't use locs if we have no information about platform type
        useLocs <- FALSE;
    }       
    
    ## Adjust the squareSize for SAM or BeadChip (only if it hasn't been set manually)
    if(is.null(squareSize)) {
        if( !is.null(BLData@experimentData$platformClass) && grepl("Matrix", BLData@experimentData$platformClass) ) 
            squareSize = 20
        else 
            squareSize = 50
    }
        
    
    if( (file.exists(locsFileName)) && (useLocs == TRUE) ) {
        
        locs <- readLocsFile(locsFileName);

        sdf <- beadarray:::simpleXMLparse(readLines(file.path(BLData@sectionData$Targets$directory[1], list.files(as.character(BLData@sectionData$Targets$directory[1]), pattern = ".sdf")[1]), warn = FALSE))
        
        nSegs <- as.integer(sdf$RegistrationParameters$SizeBlockY[[1]]);
        nRows <- as.integer(sdf$RegistrationParameters$SizeGridX[[1]]);
        nCols <- as.integer(sdf$RegistrationParameters$SizeGridY[[1]]);
        beadsPerSeg <- nRows * nCols;

        comb <- BeadDataPackR:::combineFiles(cbind(BLData[[array]][,"ProbeID"], data, BLData[[array]][,c("GrnX", "GrnY")]), locs)
        comb[which(comb[,1] == 0),2] <- NA;
        comb <- comb[order(comb[,5]),]

        resList <- segList <- list()
        
        ## process each segment individually
        for(j in 1:nSegs) {
            seg <- comb[(((j-1) * beadsPerSeg) + 1):(j * beadsPerSeg),];
            seg[,3] <- seg[,3] - min(seg[,3]);
            seg[,4] <- seg[,4] - min(seg[,4]);
            segList[[j]] <- seg;
        }
        
       # nRow <- floor(min(unlist(lapply(segList, FUN = function(seg) { return(max(seg[,3])) }))) %/% squareSize);
        nRow <- max(segList[[1]][,3]) %/% squareSize;
       
        for(j in 1:length(segList)) {
            seg <- segList[[j]];
            tmpXsize <- max(seg[,3]) / nRow;
            
            splitVals <- list(as.character(seg[,3] %/% tmpXsize), as.character(seg[,4] %/% squareSize));           
            tmp1 <- split(seg[,2], splitVals)
            tmp2 <- split(seg[,3], splitVals)
            tmp3 <- split(seg[,4], splitVals)
            
            resList[[j]] <- res <- matrix(NA, ncol = max(seg[,3] %/% tmpXsize) + 1, nrow = max(seg[,4] %/% squareSize) + 1)
            
            tmp.mean <- unlist(lapply(tmp1, mean, na.rm = TRUE))
            for(i in 1:length(tmp.mean)) {
                resList[[j]][(tmp3[[i]][1] %/% squareSize) + 1, (tmp2[[i]][1] %/% tmpXsize) + 1] <- tmp.mean[i];
            }
        }
        res <- matrix(NA, ncol = ncol(resList[[1]]), nrow = sum(unlist(lapply(resList, nrow))) + length(resList) - 1);
        rowIdx <- 1;
        for(i in 1:length(resList)) {
            for(j in 1:nrow(resList[[i]])) {
                res[rowIdx,] <- resList[[i]][j,];
                rowIdx <- rowIdx + 1;
            }
            rowIdx <- rowIdx + 1;
        }    
    }
    ## if the locs file can't be found, plot without knowledge of array segments
    else {
        
        xdat <- BLData[[array]][,"GrnX"]
        ydat <- BLData[[array]][,"GrnY"]
        #shift coordinates for plot
        xdat <- xdat - min(xdat)
        ydat <- ydat - min(ydat)
        
        splitVals <- list(as.character(xdat %/% squareSize), as.character(ydat %/% squareSize));
        tmp1 <- split(data, splitVals);
        tmp2 <- split(xdat, splitVals);
        tmp3 <- split(ydat, splitVals);

        res <- matrix(NA, ncol = max(xdat %/% squareSize) + 1, nrow = max(ydat %/% squareSize) + 1)

        tmp.mean <- unlist(lapply(tmp1, mean, na.rm = TRUE))
        for(i in 1:length(tmp.mean)) {
            res[(tmp3[[i]][1] %/% squareSize) + 1, (tmp2[[i]][1] %/% squareSize) + 1] <- tmp.mean[i]
        }
    }

    ## if we want the long axis to be vertical do so
    if(!horizontal)
        res <- t(res[nrow(res):1,])

    ## set colours
    if (is.character(low)) 
        low = col2rgb(low)/255
    if (is.character(high)) 
        high = col2rgb(high)/255
    if (!is.null(low) && is.null(high)) 
        high = c(1, 1, 1) - low
    if (is.null(low) && !is.null(high)) 
        low = c(1, 1, 1) - high
    if (is.null(low)) 
        low = c(1, 0.84, 0)
    if (is.null(high)) 
        high = c(0, 0, 1)

    col = rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
          high[2], len = ncolors), seq(low[3], high[3], len = ncolors))

    zr = NULL	
    zr[1] = min(zr[1], res[!is.na(res)], na.rm=TRUE)
    zr[2] = max(zr[2], res[!is.na(res)], na.rm=TRUE)

    if(!is.null(zlim)) {
        res[!is.na(res)] = pmax(zlim[1], res[!is.na(res)], na.rm=TRUE)
        res[!is.na(res)] = pmin(zlim[2], res[!is.na(res)], na.rm=TRUE)
    }
    else
        zlim=range(res, na.rm=TRUE)	

    ## create the plot
    image(res, col = col, xaxt = "n", yaxt = "n", zlim=zlim,...)
    if(legend)
        mtext(paste("z-range ",round(zr[1],1)," to ",round(zr[2],1)," (saturation ",round(zlim[1],1),", ",round(zlim[2],1),")",sep=""),side=1,cex=0.6)
}
