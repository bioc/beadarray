setMethod("channel",
    signature(object = "ExpressionSetIllumina"),
    function (object, name, ...) 
    {
        allNames = channelNames(object)

	if(!name %in% allNames) stop("No channel of name ", name, " was not found. Use channelNames function to see what names are valid")

	selArray = which(object@channelData[[1]] == name)	
		
	BSData = new("ExpressionSetIllumina")
	assayData(BSData)=assayDataNew(exprs = exprs(object)[,selArray], se.exprs = se.exprs(object)[,selArray],nObservations=nObservations(object)[,selArray],storage.mode="list")
	
	#assayData(BSData)=assayDataNew(exprs = exprs(object)[,selArray], se.exprs = se.exprs(object)[,selArray],storage.mode="list")

	##Create new channel information

        ##Remove channel name from all assayData matrices 

	for(el in assayDataElementNames(BSData)){

		colnames(assayDataElement(BSData, el)) = gsub(paste(name, ":",sep=""), "", colnames(assayDataElement(BSData, el)))

	}
	

	
	cData = NULL

	cData[[1]] = object@channelData[[1]][which(object@channelData[[1]] == name)]
	cData[[2]] = object@channelData[[2]][[which(channelNames(object) == name)]]

	BSData@channelData = cData


		


	BSData@phenoData = object@phenoData
	BSData@featureData = object@featureData	
	BSData@QC = object@QC
	BSData
	}
    
)

