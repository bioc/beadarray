logGreenChannelTransform = function(BLData, array){

	x = getBeadData(BLData, array=array,what="Grn")

        log2.na(x)
}

greenChannelTransform = function(BLData, array){

	x = getBeadData(BLData, array=array,what="Grn")
	return(x)
}



logRedChannelTransform = function(BLData,array){

	x = getBeadData(BLData, array=array,what="Red")

        log2.na(x)


}


redChannelTransform = function(BLData,array){

	x = getBeadData(BLData, array=array,what="Red")
	return(x)
}


logRatioTransform = function(BLData, array=array){

	x = getBeadData(BLData, array=array,what="Grn")
	
	y = getBeadData(BLData, array=array,what="Red")	

        log2.na(x) - log2.na(y)
	
}




illuminaOutlierMethod= function(inten, probeList,n= 3)
 {

    probes = sort(unique(probeList[probeList > 0]))

    nasinf = is.na(inten) | !is.finite(inten)

    inten = inten[!nasinf]
    probeList = probeList[!nasinf]
    nbeads = length(inten)
    start = 0
    foo <- .C("findAllOutliers", as.double(inten), binStatus = integer(length = nbeads), 
        as.integer(probeList), as.integer(probes), as.integer(length(probes)), 
        as.integer(nbeads), as.integer(start), as.double(n), 
        PACKAGE = "beadarray")
    sel = which((probeList > 0) & (foo$binStatus == 0))
    which(!nasinf)[sel]
}



