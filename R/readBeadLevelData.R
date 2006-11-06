"readBeadLevelData" <-
  function(targets, path=NULL, imageManipulation = "sharpen", backgroundSize=17, columns=list(ProbeID="Code", x="x", y="y"), sep = "\t", numrow = NULL){
	
   if(is.null(targets$Image1)) stop("Error: beadTargets object must contain Image1 column")
   if(is.null(targets$xyInfo)) stop("Error: beadTargets object must contain xyInfo column")

   manip = switch(imageManipulation, none = 0, sharpen = 1, sasik = 2, sasikFaster = 3, 4)

   if(manip == 4){
     stop("The imageManipulation arguement must be one of: \"none\", \"sharpen\" or \"sasik\"")
   }
   
   foregroundCalc = "Illumina"

   if(foregroundCalc == "Illumina"){
     fground = 0
   }
   else if(foregroundCalc == "sasik"){
     fground = 1
   }
     
	
   tifFiles = as.character(targets$Image1)
   csv_files = as.character(targets$xyInfo)

     if(!is.null(targets$Image1))
    tifFiles2 = as.character(targets$Image2)
  	

  csvNcol = ncol(read.table(csv_files[1], sep = sep, header = T, 
    nrows = 1))

  if (is.null(numrow)) {
    if (csvNcol == 4)
      WHAT = list(integer(0), NULL, NULL, NULL)

    else if(csvNcol == 7)
      WHAT = list(integer(0),NULL, NULL, NULL, NULL, NULL, NULL)
    
    else {
      stop("Unknown format for xyFile")
    }
    
    ffun <- function(x)
      length(scan(x, what = WHAT, sep = sep, skip = 1, quiet = TRUE)[[1]])
    
    r = sapply(as.list(csv_files),ffun)
        
    numrow = max(r)
  }
  else
    {
      if(length(numrow) != csvNcol)
        stop(paste("numrow must be a numeric vector of length ",csvNcol,", on value for each sample."))

      r = numrow
    }
  
 

   

   if(csvNcol == 4){
     BLData <- new("BeadLevelList",R=matrix(nrow=0, ncol=0), Rb=matrix(nrow=0, ncol=0),G = matrix(nrow = numrow, ncol=k), Gb = matrix(nrow = numrow, ncol=k),
                    GrnX = matrix(nrow = numrow, ncol=k), GrnY = matrix(nrow = numrow, ncol=k),
                    ProbeID = matrix(nrow = numrow, ncol=k))
   }
   else if(csvNcol == 7){
     BLData <- new("BeadLevelList",G = matrix(nrow = numrow, ncol=k), Gb = matrix(nrow = numrow, ncol=k),
                    R = matrix(nrow = numrow, ncol=k), Rb = matrix(nrow = numrow, ncol=k),
                    GrnX = matrix(nrow = numrow, ncol=k), GrnY = matrix(nrow = numrow, ncol=k),
                    RedX = matrix(nrow = numrow, ncol=k), RedY = matrix(nrow = numrow, ncol=k),
                    ProbeID = matrix(nrow = numrow, ncol=k,0))
   }
  


   for(i in 1:k){

     file=csv_files[i]
     
     if(!is.null(path)) file=file.path(path, file) 

     if(csvNcol == 4){
       dat1 <- scan(file, what = list(ProbeID = integer(0), NULL, GrnX = numeric(0),
                            GrnY = numeric(0)), sep = sep, skip = 1, quiet = TRUE)
     }
     else if(csvNcol == 7){
       dat1 <- scan(file, what = list(ProbeID = integer(0), NULL, GrnX = numeric(0),
                            GrnY = numeric(0), NULL, RedX = numeric(0), RedY = numeric(0)),
                    sep = sep, skip = 1, quiet = TRUE)
     }
     ord <- order(dat1$ProbeID)
     


    BLData@ProbeID[1:r[i], i] <- dat1$ProbeID[ord]

     numBeads = length(dat1$GrnX)
     
     greenIntensities <- .C("readBeadImage", as.character(tifFiles[i]), as.double(dat1$GrnX[ord]),
                       as.double(dat1$GrnY[ord]), as.integer(numBeads), foreGround = double(length = numBeads),
                       backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
                       as.integer(fground), PACKAGE = "beadarray")

    BLData@G[1:r[i], i] <- greenIntensities[[5]]
    BLData@Gb[1:r[i], i] <- greenIntensities[[6]]
    BLData@GrnX[1:r[i], i] <- (dat1$GrnX[ord] - min(dat1$GrnX))
    BLData@GrnY[1:r[i], i] <- (dat1$GrnY[ord] - min(dat1$GrnY))


     rm(greenIntensities)
     gc()

     if(csvNcol == 7){
       redIntensities <- .C("readBeadImage", as.character(tifFiles2[i]), as.double(dat1$RedX[ord]),
                            as.double(dat1$RedY[ord]), as.integer(numBeads), foreGround = double(length = numBeads),
                            backGround = double(length = numBeads), as.integer(backgroundSize), as.integer(manip),
                            as.integer(fground), PACKAGE = "beadarray")

       BLData@R[1:r[i], i] <- redIntensities[[5]]
      BLData@Rb[1:r[i], i] <- redIntensities[[6]]
      BLData@RedX[1:r[i], i] <- (dat1$RedX[ord] - min(dat1$RedX))
      BLData@RedY[1:r[i], i] <- (dat1$RedY[ord] - min(dat1$RedY))
      rm(redIntensities)


     }
     rm(dat1)
     gc()
     
   }

   BLData@targets = targets

   BLData
   
 }                                                              
