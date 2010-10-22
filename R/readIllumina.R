readIllumina <- function(dir= ".", useImages = FALSE, illuminaAnnotation=NULL, sectionNames = NULL, metricsFile = NULL,...) 
{
	
    dir <- normalizePath(dir);
    
	if(!is.null(sectionNames)){
	###User has specified which section names to read
		dirFiles = dir(dir)	
		
		
		txtNames = paste(sectionNames, ".txt", sep="")
		locsNames = paste(sectionNames,"_Grn.locs", sep="")
		xmlNames = paste(sectionNames, "_Grn.xml", sep="")
		tifNames = paste(sectionNames, "_Grn.tif", sep="")

		txtNames[which(!txtNames %in% dirFiles)] = NA
		locsNames[which(!locsNames %in% dirFiles)] = NA
		xmlNames[which(!xmlNames %in% dirFiles)] = NA
		tifNames[which(!tifNames %in% dirFiles)] = NA
			
		validNames = which(!is.na(txtNames))

		targets = data.frame(directory = dir, sectionName = sectionNames[validNames], textFile =txtNames[validNames], greenImage = tifNames[validNames], locs = locsNames[validNames], xmlNames = xmlNames[validNames])
		##Try to read the metrics file
			
		if(!is.null(metricsFile)){
			metrics = read.table(metricsFile, sep="\t", header=TRUE)

			###Try and match up the metrics to those we have read in 

			metricsNames = paste(metrics[,2], metrics[,3], sep="_")

			if(any(sectionNames %in% metricsNames)){
				metMat = match(sectionNames, metricsNames)
				metMat[!is.na(metMat)]
				metrics = metrics[metMat,]
			}

		}
		else metrics = NULL
	}


	else{
		##infer targets from contents of directory

    		targets <- createTargetsFile(dir, nochannels = NULL)
    		metrics = targets$metrics 	
    		targets = targets$targets
    
	}
  
    ## if there's an .sdf file, read it 
	sdf = NULL
	sdfName = list.files(dir, pattern=".sdf")	
    if(length(sdfName)){ 
        sdf <- simpleXMLparse(readLines(paste(dir, sdfName, sep = .Platform$file.sep), warn = FALSE))	
    }
    nSections <- nrow(targets);

    ## report how many channels there are
    nChannels <- numberOfChannels(paste(targets$directory[1], targets$textFile[1], sep = .Platform$file.sep), sep = "\t");
    
    BLData <- new(Class = "beadLevelData");

    BLData = insertSectionData(BLData, what = "Targets", data=targets)
    if(!is.null(metrics)) BLData = insertSectionData(BLData, what="Metrics", data = metrics)


    if(!is.null(sdf)){
        BLData@experimentData$sdfFile <- paste(dir, sdfName, sep= .Platform$file.sep)
        BLData@experimentData$platformClass <- sdf$Class[[1]];
    }



##    BLData@sectionData <- targets[,1:2];        
    nBeads <- vector(length = nSections);

    for(i in 1:nSections) {
        
        message(paste("Processing section ", targets$sectionName[i], sep = ""));
        
        data <- readBeadLevelTextFile(file.path(targets$directory[i], targets$textFile[i]),...);
    
        ##record the ProbeIDs, X and Y coords
	BLData <- insertBeadData(BLData, array = i, what = "ProbeID", data = data[,1])
        BLData <- insertBeadData(BLData, array = i, what = "GrnX", data = data[,3])
        BLData <- insertBeadData(BLData, array = i, what = "GrnY", data = data[,4])

        ## record the number of decoded beads
        nBeads[i] <- nrow(data);
        
        ## read the green images
        if(useImages && !is.null(targets$greenImage[i])) {
            greenImage <- readTIFF(fileName = as.character(targets$greenImage[i]), path = as.character(targets$directory[i]));
            ## there are wrapper functions for these, but using .Call doesn't require
            ## copying the data in the function call
            bg <- .Call("illuminaBackground", greenImage, data[,3:4], PACKAGE = "beadarray")
            greenImage <- .Call("illuminaSharpen", greenImage, PACKAGE = "beadarray");
            fg <- .Call("illuminaForeground", greenImage, data[,3:4], PACKAGE = "beadarray");
            rm(greenImage);
            
            BLData <- insertBeadData(BLData, array = i, what = "Grn", data = fg - bg)
            BLData <- insertBeadData(BLData, array = i, what = "GrnF", data = fg)
            BLData <- insertBeadData(BLData, array = i, what = "GrnB", data = bg)
              
        }
        ## or extract the data from the .txt file
        else {
            BLData <- insertBeadData(BLData, array = i, what = "Grn", data = data[,2])
        }
        
        ## if this is two channel, read the red data too
        if(nChannels == 2) {
        
            BLData <- insertBeadData(BLData, array = i, what = "RedX", data = data[,3])
            BLData <- insertBeadData(BLData, array = i, what = "RedY", data = data[,4])
            
            ## read the images
            if(useImages && !is.null(targets$redImage[i])) {
                image <- readTIFF(fileName = as.character(targets$redImage[i]), path = as.character(targets$directory[i]));
                ## there are wrapper functions for these, but using .Call doesn't require
                ## copying the data in the function call
                bg <- .Call("illuminaBackground", image, data[,6:7], PACKAGE = "beadarray")
                image <- .Call("illuminaSharpen", image, PACKAGE = "beadarray");
                fg <- .Call("illuminaForeground", image, data[,6:7], PACKAGE = "beadarray");
                rm(image);
                
                BLData <- insertBeadData(BLData, array = i, what = "Red", data = fg - bg)
                BLData <- insertBeadData(BLData, array = i, what = "RedF", data = fg)
                BLData <- insertBeadData(BLData, array = i, what = "RedB", data = bg)
                
            }
            ## or extract the data from the .txt file
            else {
                BLData <- insertBeadData(BLData, array = i, what = "Red", data = data[,5])
            }
        }
    }
    
    ## incorporate things into sectionData slot
    ## number of beads
   		
## sample groupings
    sampleGroup <- vector(length = nrow(targets))

  	if(!is.null(sdf)){	
    		tmp <- sapply(sdf$SampleLabels$string[[1]], grep, x = targets$sectionName)


    	if(is.matrix(tmp)){
		###Each sample is one more than one section
    		for(i in 1:ncol(tmp)){
        		sampleGroup[tmp[,i]] <- colnames(tmp)[i]
		}

    	}

		else{
		##One sample per section
			for(i in 1:length(tmp)){
				sampleGroup[as.numeric(tmp[i])] <- names(tmp)[i] 
    			    
			}
		}

	
	}

	BLData = insertSectionData(BLData, what="SampleGroup", data = data.frame(SampleGroup = sampleGroup))
	BLData = insertSectionData(BLData, what="numBeads", data=data.frame(numBeads = nBeads))

 	if(!is.null(illuminaAnnotation)){
		BLData = setAnnotation(BLData, illuminaAnnotation)
	}   	
	else warning("No Illumina annotation was specified. Try setting manually using setAnnotation..\n")


    return(BLData);

}
