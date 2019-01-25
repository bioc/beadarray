expressionQCPipeline = function(BLData, transFun = logGreenChannelTransform, qcDir = "QC", plotType = ".jpeg", horizontal = TRUE,controlProfile=NULL,overWrite=FALSE,nSegments=9,outlierFun=illuminaOutlierMethod,tagsToDetect = list(housekeeping = "housekeeping", Biotin = "biotin", Hybridisation = "cy3_hyb"),zlim=c(5,7),positiveControlTags = c("housekeeping", "biotin"), hybridisationTags =  c("cy3_hyb"), negativeTag= "negative", boxplotFun = logGreenChannelTransform, imageplotFun = logGreenChannelTransform){


an = sectionNames(BLData)

###Use the controlProfile if specified

if(is.null(controlProfile)){

	controlProfile = makeControlProfile(annotation(BLData))



}


if(is.null(controlProfile)){
      
  message("ControlProfile could not be created\n")
}


else{
  ##First step doing per-array plots

  dir.create(qcDir, showWarnings=F)
  dir.create(paste(qcDir, "/outliers",sep=""), showWarnings=F)
  dir.create(paste(qcDir, "/imageplot",sep=""), showWarnings=F)
  dir.create(paste(qcDir, "/controls",sep=""), showWarnings = F)

  for(i in 1:length(an)){


	  cat("Making per-array plots for section", i, "\n")
	  
	  ##Create a HTML page to write to

	  ###Plots of control probes

	  ##Positive controls
	  
	  cat("Positive controls\n")


	  if(plotType == ".jpeg"){

		  fname.pc = paste(qcDir, "/controls/", an[i], ".jpeg",sep="")

	  }

	  else if(plotType == ".png"){
		  fname.pc = paste(qcDir, "/controls/", an[i], ".png",sep="")
	  
	  }

	  if(file.exists(fname.pc)){

		  if(overWrite){
		  p <- combinedControlPlot(BLData, array=i,controlProfile = controlProfile)
		  ggsave(p, filename = fname.pc,width=8,height=8, dpi=100)
		  }
		  
		  else cat("Positive control plot exists. Skipping to next plot\n")	

	  }

	  else {
	      p <- combinedControlPlot(BLData, array=i,controlProfile = controlProfile)
		  ggsave(p, filename = fname.pc,width=8,height=8, dpi=100)
	  }
	  
	  


	  

	  cat("Outliers\n")

	  if(plotType == ".jpeg") {
		  fname.out = paste(qcDir, "/outliers/", an[i], ".jpeg",sep="")
		  
		  jpeg(fname.out, width=1200, height=300)
	  }


	  else if(plotType == ".pdf"){
		  fname.out = paste(qcDir, "/outliers/", an[i], ".pdf",sep="")
		  pdf(fname.out, width=12, height=3)
	  }
	  else if(plotType == ".png"){
		  fname.out = paste(qcDir, "/outliers/", an[i], ".png",sep="")	
		  png(fname.out, width=1200, height=300)
	  
	  }


	  if(file.exists(fname.out)){
		  if(overWrite){
			  outlierplot(BLData, array=i, nSegments = nSegments, horizontal = horizontal, outlierFun=outlierFun)
	  
		  }

		  else cat("Outlier plot exists. Skipping to next plot\n")	

	  }

	  else{
		  outlierplot(BLData, array=i, nSegments = nSegments, horizontal = horizontal, outlierFun = outlierFun)
	  }

	  dev.off()



	  cat("imageplot\n")
	  
	  
	  if(plotType == ".jpeg") {
		  fname.im = paste(qcDir, "/imageplot/", an[i], ".jpeg",sep="")
		  jpeg(fname.im, width=1200, height=300)
	  }


	  else if(plotType == ".pdf"){
		  fname.im = paste(qcDir, "/imageplot/", an[i], ".pdf",sep="")
		  pdf(fname.im, width=12, height=3)
	  }
	  else if(plotType == ".png"){
		  fname.im = paste(qcDir, "/imageplot/", an[i], ".png",sep="")	
		  png(fname.im, width=1200, height=300)
	  
	  }


	  if(file.exists(fname.im)){
		  if(overWrite){

			  im <- imageplot(BLData, array=i, useLocs=TRUE,zlim=zlim, horizontal = horizontal, transFun = imageplotFun)	
			  ggsave(im, filename = fname.im,width=4,height=1)
		  }

		  else cat("Positive control plot exists. Skipping to next plot\n")	

	  }				
	  
	  else{
		  im <- imageplot(BLData, array=i, useLocs=TRUE,zlim=zlim, horizontal = horizontal, transFun = imageplotFun)	
		  ggsave(im, filename = fname.im,width=4,height=1)

		  
	  }




	  if(requireNamespace("hwriter")) {

	      ##Make the HTML page
	      outfile = openPage(filename = paste(qcDir, "/",an[i], ".htm", sep=""))
	      
	      hwrite(paste("Quality assessment for ", an[i]), heading=1,outfile)

	      hwrite("Imageplot", heading=2,outfile)

	      hwrite("Imageplot created from the log2 transformed green intensiites. White space indicates beads that could not be decoded after array manufacture", outfile)

	      hwriteImage(gsub(paste(qcDir, "/", sep=""), "", fname.im),outfile)

	      hwrite("Outlier locations", heading=2,outfile)

	      hwrite("Locations of beads that are flagged as outliers using Illumina's outlier detection procedure on log2 intensities", outfile)

	      hwriteImage(gsub(paste(qcDir, "/", sep=""), "",fname.out),outfile)

	      hwrite("Positive Controls", heading=2,outfile)

	      hwriteImage(gsub(paste(qcDir, "/", sep=""), "",fname.pc), outfile)

	      closePage(outfile)
	  }


	  else warning("Could not create HTML page. Make sure that 'hwriter' package is installed\n")

  }

	  if(requireNamespace("hwriter")) {
	      outfile = openPage(filename = paste(qcDir, "/Summary.htm", sep=""))


	      hwrite("Quality assessment summary", heading=1, outfile)


	      ##Create boxplot using defined functions

	      if(plotType == ".jpeg") {jpeg(paste(qcDir, "/Boxplot.jpeg",sep=""), width = 1200, height = 300);hwriteImage("Boxplot.jpeg", outfile)}
	      if(plotType == ".png") {png(paste(qcDir, "/Boxplot.png",sep=""), width = 1200, height = 300);hwriteImage("Boxplot.png", outfile)}
	      if(plotType == ".pdf") pdf(paste(qcDir, "/Boxplot.pdf",sep=""), width = 12, height = 3);

	      boxplot(BLData, transFun = boxplotFun, outline=FALSE)
		      
	      dev.off()
		      


	      hwrite("Scan Metrics", heading=2,outfile)

	      if("Metrics" %in% colnames(BLData@sectionData)) hwrite(BLData@sectionData$Metrics, outfile)

	      hwrite("Bead-level control summary", heading=2, outfile)


	      cat("Creating probe metrics\n")

	      beadLevelQC = makeQCTable(BLData, transFun = transFun, controlProfile = controlProfile)


	      hwrite(beadLevelQC, outfile)

	      cat("Calculating outlier Metrics\n")

	      outlierTable = matrix(nrow = length(an), ncol = nSegments)

	      colnames(outlierTable) = paste("Segment", 1:nSegments)
	      rownames(outlierTable) = an

	      for(i in 1:length(an)){

		      outlierTable[i,] = calculateOutlierStats(BLData, transFun = transFun, array=i,nSegments=nSegments, outlierFun=outlierFun)

	      }
	      

		  hwrite("Outlier Metrics", outfile,heading=2)


		  hwrite(round(outlierTable,2), outfile)


			  detectionTable = matrix(nrow = length(an), ncol=length(tagsToDetect))
			  colnames(detectionTable) = names(tagsToDetect)
			  rownames(detectionTable) = an
			  for(i in 1:length(an)){
				  detectionTable[i,] = controlProbeDetection(BLData, transFun = transFun, array=i, tagsToDetect = tagsToDetect, negativeTag = negativeTag, controlProfile=controlProfile)
			  }

				  

		  hwrite("Detection Metrics", outfile,heading=2)

		  hwrite(round(detectionTable, 2),outfile)




		  closePage(outfile)




		  ##Write to csv


		  write.csv(beadLevelQC, file=paste(qcDir,"/probeMetrics.csv",sep=""), quote=FALSE)

		  write.csv(metrics(BLData), file=paste(qcDir, "/scanMetrics.csv",sep=""), quote=FALSE)	

		  write.csv(outlierTable, file=paste(qcDir, "/outlierMetrics.csv",sep=""), quote=FALSE)

		  write.csv(detectionTable, file=paste(qcDir, "/detectionMetrics.csv",sep=""), quote=FALSE)
	  }	

	  else warning("Could not create HTML page. Make sure that 'hwriter' package is installed\n")


  }

}
