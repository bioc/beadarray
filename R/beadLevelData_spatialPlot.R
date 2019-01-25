
showArrayMask2 <- function(BLData,array = 0, override = FALSE, wtsName = "wts", transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, horizontal=TRUE)
{
	if(class(BLData) != "beadLevelData")
	{stop("Must be performed on a beadLevelData object.")}

	if(array == 0)
	{
		stop(paste("Please specify an array."))
	}

	an = sectionNames(BLData)

	if(wtsName %in% colnames(BLData[[array]]))


	{	wtCol = which(colnames(BLData[[array]]) == wtsName)
	
		sel <- which(BLData[[array]][,wtCol]==0)

		##get out now if there are too many masked beads
		if(!override & length(sel) > 200000)
		{
			stop(paste("There are over 200 000 beads in the mask, thus plotting the mask may cause R to freeze. (You can override this error by setting override = TRUE.)\nNumber of masked beads:",length(sel),"\n"))
		}

		transInten = transFun(BLData, array=array)	

		pidCol = which(colnames(BLData[[array]]) == "ProbeID")

		o <- outlierFun(transInten, BLData[[array]][,pidCol])

		plotBeadLocations2(BLData, array = array, BeadIDs = c(o,sel), main = an[array], pch = ".",horizontal=horizontal,ptsCol=c(rep("outlier", length(o)), rep("masked", length(sel))))



		#if(elim)
		#{
		#	eliminated <- listEliminatedProbes(BLData, 1)
		#	sel <- which(BLData[[an[array]]]$ProbeID %in% eliminated)
		#	x.cds <- BLData[[an[array]]]$GrnX[sel]
		#	y.cds <- BLData[[an[array]]]$GrnY[sel]
		#	points(x.cds, y.max - y.cds, pch = "X", col = "Blue")			
		#}
	}
}



plotBeadLocations2 = function(BLData, ProbeIDs=NULL, BeadIDs=NULL, array=1, SAM=FALSE, xCol = "GrnX", yCol="GrnY", xlab="x-coordinate", ylab="y-coordinate", horizontal = TRUE, main=paste("Bead", ProbeIDs, "locations"),ptsCol=NULL,...){


xs.orig = getBeadData(BLData, array=array, what="GrnX")
ys.orig = getBeadData(BLData, array=array, what="GrnY")

##If plotting with the longest edge going along the screen, flip the input coordinates

if(horizontal){

	tmp=ys.orig

	ys = xs.orig

	xs = tmp
	rm(tmp)

}

##If keeping the same orientation as the image, flip the y coordinates so that y=0 is at bottom of screen

else{

	xs = xs.orig
	ys = max(ys.orig) - ys.orig 

}


##xs = max(xs) -xs



##Put the origin at (0,0)

xs = xs - min(xs)
ys = ys - min(ys)

xmax = max(xs)
ymax = max(ys)



#plot(1, xlim=range(0:xmax), ylim=range(0:ymax) ,type="n", new=TRUE, axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i",main=main,...)

if(!(is.null(ProbeIDs))){

	pIDs = getBeadData(BLData, what="ProbeID",array=array)


	xs = xs[which(pIDs %in% ProbeIDs)]
	ys = ys[which(pIDs %in% ProbeIDs)]

}
else{

	xs = xs[BeadIDs]
	ys = ys[BeadIDs]
 
}

p <- qplot(x=xs, y=ys, size=I(0.3),colour=ptsCol)

p+opts(axis.line = theme_blank(), axis.text.x = theme_blank(), 
            axis.text.y = theme_blank(), axis.ticks = theme_blank(), 
            axis.title.x = theme_blank(), axis.title.y = theme_blank(), 
            legend.position = "none", panel.background = theme_blank(), 
            panel.border = theme_blank(), panel.grid.major = theme_blank(), 
            panel.grid.minor = theme_blank(), plot.background = theme_blank())
#points(xs, ys,...)
#box()

}


outlierplot2 <- function(BLData, array=array, transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, wtsname=NULL,horizontal = TRUE, nSegments = NULL, lowOutlierCol = "blue", highOutlierCol = "pink", outlierPch = ".", main="",...){

    locsFileName <- file.path(BLData@sectionData$Targets$directory[array], paste(BLData@sectionData$Targets$sectionName[array], "_Grn.locs", sep = ""))

    ##Find all outliers on the array
        
    wts<-1
    if(!is.null(wtsname)){wts<-getBeadData(BLData,array=array,what=wtsname)}
    
    oList = outlierFun(transFun(BLData, array=array), probeList = getBeadData(BLData,array=array,what="ProbeID"),wts = wts,...)

    cat(length(oList), " outliers found on the section\n")

    beadMeans = lapply(split(getBeadData(BLData, what="Grn", array=array), getBeadData(BLData, what="ProbeID", array=array)), mean,na.rm=TRUE)

    resids = getBeadData(BLData, what="Grn",array=array) - unlist(beadMeans[match(getBeadData(BLData, what="ProbeID", array=array), names(beadMeans))])

    oCols = NULL

    oCols[which(resids > 0 )] = "high"
    oCols[which(resids < 0 )] = "low"
    oCols[which(resids == 0 )] = "low"


    p <- plotBeadLocations2(BLData, array=array, BeadIDs = oList, horizontal = horizontal, col=oCols,pch=outlierPch,main=main, ptsCol=oCols[oList])

    ##Work out the position of segments

    ## can we read the sdf file?
    sdfFileName <- file.path(BLData@sectionData$Targets$directory[1], list.files(as.character(BLData@sectionData$Targets$directory[1]), pattern = ".sdf")[1]);
    if( file.exists(sdfFileName) && (is.null(nSegments)) ){
        sdf <- simpleXMLparse(readLines(sdfFileName, warn = FALSE))
        nSegments <- as.integer(sdf$RegistrationParameters$SizeBlockY[[1]]);
    }
    
    ## only plot the segments if we have a value and this is a BeadChip rather than a SAM
    if( !is.null(nSegments) && !is.null(BLData@experimentData$platformClass) && !grepl("Matrix", BLData@experimentData$platformClass) ) {
        ys = getBeadData(BLData, what="GrnY", array=array)
        ys = ys - min(ys)
                
        segEnds = seq(from=0, to = max(ys), by = max(ys)/(nSegments + 1))

        if(horizontal) 
            p <- p +  geom_vline(xintercept=segEnds, colour="green")
        else 
           p <- p + geom_hline(yintercept=segEnds, colour="green") 
    }

  p

}



	



