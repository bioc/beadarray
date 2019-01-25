##Use normalize.invariant set to find a list of genes which do not change rank across arrays and then normalise to the 
##previously defined target distrbution


rankInvariantNormalise = function(exprs, T=NULL){
requireNamespace("affy")
if(is.null(T)){

T = apply(exprs, 1, mean, na.rm=TRUE)

}


for (i in 1:ncol(exprs)){

curve = normalize.invariantset(exprs[,i], T)

exprs[,i] = predict(curve$n.curve, exprs[,i])$y
}

exprs

}





medianNormalise = function(exprs, log=TRUE){

exprs = as.matrix(exprs)

narrays = ncol(exprs)

if(log){
                          
exprs = log2(exprs)

}
                          
med = median(exprs,na.rm=TRUE)

for(i in 1:narrays){

exprs[,i]  = exprs[,i] - median(exprs[,i], na.rm=TRUE)

}


exprs = exprs + med


exprs

}


normaliseIllumina = function(BSData, method="quantile", transform="none", T=NULL, status=fData(BSData)$Status,negctrl="negative",regular="regular",...) {
  rownms = rownames(exprs(BSData))
  colnms = colnames(exprs(BSData))
  transform = match.arg(transform, c("none", "vst", "log2"))
  method = match.arg(method, c("quantile", "qspline", "vsn", "rankInvariant", "median", "none", "neqc", "rsn"))
  if(method=="vsn" && transform!="none"){
          cat(paste("\nmethod =", method, "not compatible with transform =", transform, "\nResetting transform = \"none\"\n\n"))
          transform="none"
  }
  if(transform=="log2" && method != "neqc") {
          BSData = assayDataElementReplace(BSData, "exprs", as.matrix(log2(exprs(BSData))))
  }
  else if(transform=="vst") {
          requireNamespace("lumi")
          x = new("LumiBatch")
          x = assayDataElementReplace(x, "exprs", exprs(BSData))
          x = assayDataElementReplace(x, "se.exprs", se.exprs(BSData))
          x = assayDataElementReplace(x, "beadNum", nObservations(BSData))
          if(!all(is.na(Detection(BSData))))
             x = assayDataElementReplace(x, "detection", Detection(BSData))
	     else try(x <- assayDataElementReplace(x, "detection", calculateDetection(BSData)))	
		
          BSData = assayDataElementReplace(BSData, "exprs", exprs(lumiT(x, method="vst", ...)))
          rm(x)
  }

  switch(method,
         quantile={
            BSData = assayDataElementReplace(BSData, "exprs", normalizeQuantiles(as.matrix(exprs(BSData))))
            rownames(exprs(BSData)) = rownms
            colnames(exprs(BSData)) = colnms
         },
         qspline={
            requireNamespace("affy")
            BSData = assayDataElementReplace(BSData, "exprs", normalize.qspline(as.matrix(exprs(BSData))))
            rownames(exprs(BSData)) = rownms
            colnames(exprs(BSData)) = colnms
         },
         vsn={
            requireNamespace("vsn")
            BSData = assayDataElementReplace(BSData, "exprs", exprs(vsn2(exprs(BSData))))
         },
         rankInvariant={
            BSData = assayDataElementReplace(BSData, "exprs", rankInvariantNormalise(exprs(BSData), T=T))

         },
         median={
            BSData = assayDataElementReplace(BSData, "exprs", medianNormalise(exprs(BSData), log=FALSE))
         },

	rsn={
		requireNamespace("lumi")
		##Need to check there are no NAs in the data 

		noNA <- apply(exprs(BSData), 1, function(x) !any(is.na(x)))

		newObj <- BSData

		exprs(newObj)[noNA,] = lumi::rsn(exprs(BSData)[noNA,],...)

		BSData = newObj
		

	},

	neqc={
		
	      #noNA <- apply(exprs(BSData), 1, function(x) !any(is.na(x)))

	      ##note that neqc removes any control probes when it returns it output 

	      newObj <- BSData[which(status == regular),]

#	      tmp <- neqc(exprs(BSData)[noNA,],status=status[noNA], negctrl=negctrl, regular = regular,...)
	      tmp <- neqc(exprs(BSData),status=status, negctrl=negctrl, regular = regular,...)
	      
#	      exprs(newObj) <- tmp[match(featureNames(newObj), rownames(tmp)),]
        exprs(newObj) <- tmp

	      BSData = newObj
	

	  })
	

  BSData
}
