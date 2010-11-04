
showArrayMask <- function(BLData,array = 0, override = FALSE, wtsName = "wts", transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, horizontal=TRUE)
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

		plotBeadLocations(BLData, array = array, BeadIDs = c(o,sel), main = an[array], SAM = SAM, pch = ".",horizontal=horizontal,col=c(rep("black", length(o)), rep("red", length(sel))))



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
