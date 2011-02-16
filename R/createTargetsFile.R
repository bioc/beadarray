createTargetsFile <-
function(dir=NULL, nochannels=1, channel1="Grn", channel2="Red", txtsuff="txt", imgsuff="tif", locssuff="locs", xmlsuff="xml", verbose=FALSE,special=c("sdf","fiducial"),ColourConfusionStop=T,metricsflag="Metrics",metsep="\t",metricsection="Section",metricchip="Matrix"){

    ##check number of channels makes sense
    if(!is.null(nochannels)){
        if((nochannels!=1)&(nochannels!=2)){stop("the number of channels should be NULL, 1 or 2.")}
    }


    ## has a directory been specified?
    ## if not, assume working directory
    #if(is.null(dir)){
    #dir<-getwd()
    #if(verbose){cat("No directory supplied, using working directory.\n")}}

    ## check that directory exists
    if(!file.exists(dir)){ stop("Directory does not exist.\n") }

    ## report back on directory choice
    if(verbose){ cat("Directory is",dir,".\n") }

    ## are there any files in the directory?
    filelist<-list.files(dir)
    nofiles<-length(filelist)
    if(nofiles==0){ stop("Directory is empty.\n") }
    if(verbose){ cat("Found",nofiles,"files in the directory.\n") }

    ## look for metrics files

    fmet <- grep(metricsflag,filelist)
    if(verbose){ cat("Found",length(fmet),"metrics files in the directory.\n") }

    storemet<-NULL
    for(i in fmet){
        storemet<-rbind(storemet,read.table(file.path(dir,filelist[i],fsep = .Platform$file.sep),header=T,as.is=T,sep=metsep))
    }

    if(length(fmet)) { filelist<-filelist[-fmet] }
    nofiles<-length(filelist)

    ## Remove fiducial, sdf files etc.
    for(i in special){
        remspec<-grep(i,filelist)
        if(length(remspec)>0){
            if(verbose){
                cat("Ignoring",length(remspec),"files matching",i,".\n")
            }	
            filelist<-filelist[-remspec]
        }
    }

    nofiles<-length(filelist)
    if(nofiles==0){stop("Directory is empty.\n")}
    if(verbose){cat(nofiles,"remaining for consideration.\n")}

    ## sort out channel colours
    filegreen<-rep(F,nofiles)
    filered<-rep(F,nofiles)
    filegreen[grep(channel1,filelist)] <- TRUE
    filered[grep(channel2,filelist)] <- TRUE

    ## Check the number of channels

    if(!is.null(nochannels)){
        ## if the number of channels is meant to be 1
        if(nochannels==1){
            if(ColourConfusionStop&(sum(filered)>0)){
                stop("One channel was specified, but this appears to be a two channel array.")
            }
        }

        ## if the number of channels is meant to be 2
        if(nochannels==2){
            if(ColourConfusionStop&(sum(filered)==0)){
                stop("Two channels were specified, but this appears to be a single channel array.")
            }
        }
    }

    ## if the number of channels is to be inferred
    if(is.null(nochannels)){
        if(sum(filered)>0){
            nochannels<-2
            if(verbose){cat("Number of channels was not specified, going with two.")}
        }
        if(sum(filered)==0){
            nochannels<-1
            if(verbose){cat("Number of channels was not specified, going with one.")}
        }
    }

    ##sort out file types

    filetype<-rep(NA,nofiles)
    for(i in 1:nofiles){
    temp<-strsplit(filelist[i],".",fixed=T)[[1]]
    filetype[i]<-temp[length(temp)]}

    if(verbose){
    for(j in unique(filetype)){
    cat("found",sum(filetype==j,na.rm=T),"files of type",j,"\n")
    }
    }


    ## ultimately there has to be a textfile, so we use these as the key

    ##first check there is a textfile
    keys<-filelist[which(filetype==txtsuff)]
    if(length(keys)==0){stop("There were no text files with the required suffix")}


    for(i in 1:length(keys)){
    keys[i]<-substr(keys[i],1,nchar(keys[i])-nchar(txtsuff)-1)
    }


    ## one channel 

    if(nochannels==1){
        sectionName <- rep(NA, length(keys) );
        sectiontext<-rep(NA,length(keys))
        greenimage<-rep(NA,length(keys))
        greenlocs<-rep(NA,length(keys))
        greenxml<-rep(NA,length(keys))
        for(i in 1:length(keys)){
            temp<-grep(keys[i],filelist)
            if(length(temp[filetype[temp]==txtsuff])==1){

                sectionName[i] <- strsplit(filelist[temp[filetype[temp]==txtsuff]], "\\.")[[1]][1];
		####In order to make iScan files available
		sectionName[i] = gsub("_perBeadFile", "", sectionName[i])	


                sectiontext[i]<-filelist[temp[filetype[temp]==txtsuff]]
            }
            if(length(temp[(filetype[temp]==imgsuff)&(filered[temp]==F)])==1){
                greenimage[i]<-filelist[temp[(filetype[temp]==imgsuff)&(filered[temp]==F)]]
            }
            if(length(temp[(filetype[temp]==locssuff)&(filered[temp]==F)])==1){
                greenlocs[i]<-filelist[temp[(filetype[temp]==locssuff)&(filered[temp]==F)]]
            }
            if(length(temp[(filetype[temp]==xmlsuff)&(filered[temp]==F)])==1){
                greenxml[i]<-filelist[temp[(filetype[temp]==xmlsuff)&(filered[temp]==F)]]
            }
        }
        targets<-cbind(rep(dir,length(keys)),sectionName, sectiontext,greenimage,greenlocs,greenxml)
        colnames(targets)<-c("directory","sectionName", "textFile","greenImage","locs","xml")
    }




    ## twochannel 

    if(nochannels==2){

        sectionName <- rep(NA, length(keys) );
        sectiontext<-rep(NA,length(keys))
        greenimage<-rep(NA,length(keys))
        greenlocs<-rep(NA,length(keys))
        greenxml<-rep(NA,length(keys))
        redimage<-rep(NA,length(keys))
        redlocs<-rep(NA,length(keys))
        redxml<-rep(NA,length(keys))

        for(i in 1:length(keys)){
            temp<-grep(paste(keys[i],"[_\\.]", sep = ""),filelist)

            if(length(temp[filetype[temp]==txtsuff])==1){
                sectionName[i] <- strsplit(filelist[temp[filetype[temp]==txtsuff]], "\\.")[[1]][1];
                sectiontext[i]<-filelist[temp[filetype[temp]==txtsuff]]
            }
            if(length(temp[(filetype[temp]==imgsuff)&(filegreen[temp]==T)])==1){
                greenimage[i]<-filelist[temp[(filetype[temp]==imgsuff)&(filegreen[temp]==T)]]
            }
            if(length(temp[(filetype[temp]==locssuff)&(filegreen[temp]==T)])==1){
                greenlocs[i]<-filelist[temp[(filetype[temp]==locssuff)&(filegreen[temp]==T)]]
            }
            if(length(temp[(filetype[temp]==xmlsuff)&(filegreen[temp]==T)])==1){
                greenxml[i]<-filelist[temp[(filetype[temp]==xmlsuff)&(filegreen[temp]==T)]]
            }
            if(length(temp[(filetype[temp]==imgsuff)&(filered[temp]==T)])==1){
                redimage[i]<-filelist[temp[(filetype[temp]==imgsuff)&(filered[temp]==F)]]
            }
            if(length(temp[(filetype[temp]==locssuff)&(filered[temp]==T)])==1){
                redlocs[i]<-filelist[temp[(filetype[temp]==locssuff)&(filered[temp]==F)]]
            }
            if(length(temp[(filetype[temp]==xmlsuff)&(filered[temp]==T)])==1){
                redxml[i]<-filelist[temp[(filetype[temp]==xmlsuff)&(filered[temp]==F)]]
            }

            targets<-cbind(rep(dir,length(keys)), sectionName, sectiontext,greenimage,greenlocs,greenxml,redimage,redlocs,redxml)
            colnames(targets)<-c("directory", "sectionName", "textFile","greenImage","greenLocs","greenxml","redImage","redLocs","redxml")

        }
    }


                metrow<-rep(NA,length(keys))
	
    ##match up metrics
    if(length(fmet)>0){

        if(verbose){cat("Matching up metrics file.\n")}
        if(verbose){cat("Will now look for Section and Chip columns of metrics file.\n")}

        if(!all(colnames(storemet)!=metricsection)){
            if(verbose){cat("Found section column.\n")}

            if(!all(colnames(storemet)!=metricchip)){
                if(verbose){cat("Found chip column.\n")}


		keys = gsub("_perBeadFile", "", keys)		

                for(i in 1:(dim(storemet)[1])){

                    matchsect<-grep(paste(storemet[[metricsection]][i], "$", sep = ""), keys)
                    matchchip<-grep(storemet[[metricchip]][i],keys)

                    if(sum(duplicated(c(matchchip,matchsect)))==1){

                        mymatch<-c(matchchip,matchsect)[duplicated(c(matchchip,matchsect))]

                        metrow[mymatch]<-i
                    }
                }

            }
        }
    return(list(targets=as.data.frame(targets), metrics = as.data.frame(storemet[metrow,])))
    }

    ## remove reordering for now
    ## this is too simplistic and breaks a lot!
    #targets <- targets[order(as.integer(rownames(targets))),]
    return(list(targets = as.data.frame(targets)))
}

