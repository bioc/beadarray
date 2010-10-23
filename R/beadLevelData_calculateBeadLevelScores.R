
	

plotBeadIntensities = function(BLData, array=1, BeadIDs,transFun=logGreenChannelTransform,cols=NULL,...){


	pIDs = getBeadData(BLData, array=array, what="ProbeID")

	values = transFun(BLData, array = array)

	subset = values[which(pIDs %in% BeadIDs)]		
	
	pIDs = pIDs[which(pIDs %in% BeadIDs)]


	intenList = split(subset, pIDs)
	
	if(!is.null(cols)){

		for(i in 1:length(BeadIDs)){
	
		colList = rep(cols[i], length(intenList[[i]]))
		}

	}
	
	else colList = rainbow(n=length(intenList))
	
	genericBeadIntensityPlot(intenList, colList,...)
	

}

genericBeadIntensityPlot = function(intenList,colList="black",...){

	if(length(intenList) == 0){

		stop("No intensities to be plotted\n")
	}
		
	yvals = unlist(intenList)
	
	nItems = unlist(lapply(intenList, length))

	xvals = rep(1, nItems[1]) + runif(nItems[1], 0, 0.1)
	
	cols = rep(colList[1], nItems[1])

	if(length(intenList) > 1){
	
		for(i in 2:length(intenList)){

			xvals = c(xvals, rep(i, nItems[i]) + runif(nItems[i], 0, 0.1))
			cols = c(cols, rep(colList[i], nItems[i]))
		}
	}	

	plot(xvals, yvals, col=cols,pch=16,axes=FALSE,xlab="",ylab="log2 intensity", las=2,...)
	box()
	axis(2)
	axis(1, at = 1:length(intenList), labels =names(intenList), las=2)
}




controlProbeDetection = function(BLData, transFun = logGreenChannelTransform, array = 1, controlProfile = NULL, tagsToDetect = list(housekeeping = "housekeeping", Biotin = "phage_lambda_genome", Hybridisation = "phage_lambda_genome:high"), negativeTag = "permuted_negative", detThresh=0.05){

	
	detect= function(x) 1 - (sum(x>negVals)/(length(negVals)))

##If control profile is not specified we will try and use the annotation of the beadLevelData object

	if(is.null(controlProfile)){
		anno = getAnnotation(BLData)

		controlProfile = makeControlProfile(anno)

	}

	##Remove any duplicated IDs

	controlProfile = controlProfile[-which(duplicated(controlProfile[,1])),]

	negIDs = controlProfile[which(controlProfile[,2] == negativeTag),1]

	if(length(negIDs) == 0) stop("Could not find any IDs with Tag ", negativeTag)

	pIDs = getBeadData(BLData, array=array ,what="ProbeID")

	transInten = transFun(BLData,array=array)

	negVals = transInten[which(pIDs %in% negIDs)]


	##No control tags were specified; assume we want to test all
	
	if(is.null(tagsToDetect)){
	
		transInten = transInten[which(pIDs %in% controlProfile[,1])]

		transIDs = pIDs[which(pIDs %in% controlProfile[,1])]

		transFac = controlProfile[match(transIDs, controlProfile[,1]),2]

		transInten = split(transInten, transFac)
	}

	else{
		selIDs = controlProfile[which(controlProfile[,2] %in% unlist(tagsToDetect)),1]
		
		transInten = transInten[which(pIDs %in% selIDs)]
			
		transIDs = pIDs[which(pIDs %in% selIDs)]

		transFac = as.character(controlProfile[match(transIDs, controlProfile[,1]),2])

		transInten = split(transInten, transFac)
		names(transInten) = names(tagsToDetect)[match(names(transInten), tagsToDetect)]
	}

	M = median(negVals, na.rm=TRUE)
	MAD = mad(negVals, na.rm=TRUE)

	negVals = negVals[negVals < M+3*MAD]
	negVals = negVals[!is.na(negVals)]

	resList = NULL


	for(i in 1:length(transInten)){
	
		resList[[i]] = sum(sapply(transInten[[i]],detect) < detThresh,na.rm=TRUE)/length(transInten[[i]])*100

	}	

	names(resList)= names(transInten)

	unlist(resList)

}


	
poscontPlot = function(BLData, array=1, transFun = logGreenChannelTransform, positiveControlTags = c("housekeeping", "biotin"), colList=c("red","blue"), controlProfile = NULL,...){

	data(ExpressionControlData)

	if(is.null(controlProfile)){
			
		controlInfo = ExpressionControlData[[getAnnotation(BLData)]]

	}	
	
	else controlInfo = controlProfile

	if(is.null(colList)) colList = rainbow(length(positiveControlTags))
	


	posInten = NULL

	pIDs = getBeadData(BLData, array=array, what="ProbeID")

	transInten = transFun(BLData, array=array)

	cols = NULL

	for(i in 1:length(positiveControlTags)){

		Ids = controlInfo[controlInfo[,2] == positiveControlTags[i],1]
	
		if(length(Ids) == 0){
			warning("Could not find any IDs matching the description", positiveControlTags[i])
		}

		selBeads = which(pIDs %in% Ids)

		subset = cbind(pIDs[selBeads], transInten[selBeads])
	
		posInten = append(posInten, split(subset[,2], subset[,1]))

		cols = c(cols, rep(colList[i], length(Ids)))

	}
	
	
	genericBeadIntensityPlot(posInten, colList = cols,...)


}



quickSummary = function(BLData, array=1, transFun = logGreenChannelTransform, reporterIDs = NULL, reporterTags = NULL,reporterFun = function(x) mean(x, na.rm=TRUE)){

	if(is.null(reporterIDs)){
		controlInfo = ExpressionControlData[["Humanv3"]]

		reporterIDs = controlInfo[,1]
		reporterTags = controlInfo[,3]
	}

	tmp = BLData[[array]]

	inten = transFun(BLData, array=array)

	if(any(is.infinite(inten))){
		
		probeIDs = tmp[-which(is.infinite(inten)),1]

		inten = inten[-which(is.infinite(inten))]
	}
 	
	else probeIDs = tmp[,1]

	tagFac = reporterTags[match(probeIDs, reporterIDs)]

	lapply(split(inten, tagFac), reporterFun)

	
}



makeControlProfile = function(annoName){

	if(annoName %in% names(ExpressionControlData)){

		controlProfile = data.frame(ArrayAddress = ExpressionControlData[[as.character(annoName)]][,1], Tag = 	ExpressionControlData[[as.character(annoName)]][,3])
		controlProfile
	}

	else cat("Could not make control profile for annotation:", annoName, "\n")



}


makeQCTable = function(BLData, transFun = logGreenChannelTransform, controlProfile = NULL, summaryFns = list(Mean = function(x) mean(x, na.rm=TRUE), Sd = function(x) sd(x, na.rm=TRUE)), channelSuffix = NULL){

	##If no profile was specified, use the annotation of the BLData object
	if(is.null(controlProfile)){
		anno = getAnnotation(BLData)
		
		controlProfile = makeControlProfile(anno)
	
	}

	an = sectionNames(BLData)
	
	uIDs = unique(controlProfile[,2])
	
	firstCol = 1 #The first column to be filled in

	qcTable = matrix(nrow = length(an), ncol = length(uIDs)*length(summaryFns))
	colnames(qcTable) = paste(rep(names(summaryFns), each = length(uIDs)), rep(sort(uIDs),length(names(summaryFns))),sep=":")

	if(!is.null(channelSuffix)) colnames(qcTable) = paste(colnames(qcTable), channelSuffix, sep=":")		


	
	rownames(qcTable) = an

	for(j in 1:length(summaryFns)){

		lastCol = firstCol + length(uIDs) - 1
			
		for(i in 1:length(an)){
			qcTable[i,firstCol:lastCol] = unlist(quickSummary(BLData, array=i, transFun = transFun, reporterIDs = controlProfile[,1], reporterTags = controlProfile[,2], reporterFun = summaryFns[[j]]))

																
		}

		firstCol = lastCol + 1

	
	}

	qcTable
}


