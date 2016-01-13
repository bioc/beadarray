combinedControlPlot <- function(data, array = 1, controlProfile=NULL, wtsName = NULL, negativeLabel = "negative", excludeString = "ERCC"){

#require("ggplot2")

if(is.null(controlProfile)){ 
	annoName = annotation(data)
	cProf <- makeControlProfile(annoName)


}

else cProf <- controlProfile

if(is.null(cProf)){
	  
  message("ControlProfile could not be created\n")
}

else{

  inten <- getBeadData(data, array=array, what="Grn")
  pIDs <- getBeadData(data, array=array, what="ProbeID")


  if(!is.null(wtsName)){

    if(wtsName %in% colnames(data[[array]])){

    wts <- getBeadData(data, array, what=wtsName)

  }

  }

  else wts = rep(1, length(inten))



  bsv <-  identifyControlBeads(data, array=array, cProf)

  if(!is.null(bsv)){


  excludeInd <- grep(excludeString, bsv)


  if(length(excludeInd) > 0){

	  bsv <- bsv[-excludeInd]
	  pIDs <- pIDs[-excludeInd]
	  inten <- inten[-excludeInd]
	  wts <- wts[-excludeInd]

  }

  bigGrps <-  names(which(table(bsv)>5000))



  if(length(bigGrps) >0){

	  removeInd<- NULL
	  
	  for(i in 1:length(bigGrps)){

		  len <- length(which(bsv == bigGrps[i]))

		  removeInd <- c(removeInd, sample(which(bsv == bigGrps[i]),len - 5000))

	  
	  }



	  bsv <- bsv[-removeInd]
	  pIDs <- pIDs[-removeInd]
	  inten <- inten[-removeInd]
	  wts <- wts[-removeInd]

  }


  if(any(inten < 0  | inten == 0 | is.na(inten))){

    removeInd = which(inten <0 | inten ==0 | is.na(inten))
	  bsv <- bsv[-removeInd]
	  pIDs <- pIDs[-removeInd]
	  inten <- inten[-removeInd]
	  wts <- wts[-removeInd]

    

  }
      



  df <- data.frame(ControlType = bsv, ID = pIDs, Log2Intensity = log2(inten), Masked = sapply(wts, function(x) x == 0), Negative = sapply(bsv, function(x) x == negativeLabel), Regular = sapply(bsv, function(x) x == "regular"), Control = sapply(bsv, function(x) !(x %in% c("regular", negativeLabel))))

  df.controls <- subset(df,Control)

  df.negative <- subset(df, Negative)

  d1 <- density(df.negative[,3])

  qs <- quantile(subset(df,Negative)[,3])

  p1 <- ggplot(data=subset(df, Control), aes(x = factor(ID), y = Log2Intensity, fill=factor(ControlType))) + geom_boxplot() + geom_hline(yintercept = qs[4],color="green") + geom_hline(yintercept = qs[3],color="green") + geom_hline(yintercept = qs[2],color="green") + geom_hline(yintercept = qs[1],color="green") + facet_wrap(~ControlType)

  #p1 <- ggplot(df, aes(x = factor(ID[Control]), y = Log2Intensity[Control], fill=factor(ControlType[Control]))) + geom_boxplot(data = subset(df, Control)) + geom_hline(y = qs[4],color="green") + geom_hline(y = qs[3],color="green") + geom_hline(y = qs[2],color="green") + geom_hline(y = qs[1],color="green")

  if(any(df$Masked)){

  p1 <- ggplot(data=subset(df, Control), aes(x = factor(ID), y = Log2Intensity, fill=factor(ControlType))) + geom_boxplot() + geom_point(data = subset(df.controls, Masked), color="red") + geom_hline(yintercept = qs[4],color="green") + geom_hline(yintercept = qs[3],color="green") + geom_hline(yintercept = qs[2],color="green") + geom_hline(yintercept = qs[1],color="green") + facet_wrap(~ControlType,scales="free_x")

  }

  else p1 <- ggplot(data=subset(df, Control), aes(x = factor(ID), y = Log2Intensity, fill=factor(ControlType))) + geom_boxplot() + geom_hline(yintercept = qs[4],color="green") + geom_hline(yintercept = qs[3],color="green") + geom_hline(yintercept = qs[2],color="green") + geom_hline(yintercept = qs[1],color="green") + facet_wrap(~ControlType,scales="free_x")


  p1


  }


}

}