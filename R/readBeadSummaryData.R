readBeadSummaryData<- function(targets=NULL, header=T, sep=",",path=NULL,
                                   columns = list(ProbeID = "TargetID",
                                     AvgSig = "AVG_Signal", NoBeads = "Avg_NBEADS",
                                     Detection="Detection", BeadStDev="BEAD_STDEV"),
                                   other.columns = NULL, skip = 7)
{


if(!is.null(targets)) {

if(!is.null(path)) targets = file.path(path, targets)

}


else{

if(!is.null(path)){

setwd(path)
}

targets=dir()


}

  filecounter = arraycounter = 0

  r = try(read.table(as.character(targets[1,1]), header = header, sep = sep, skip = skip, nrows = 1))
  cat("Reading file ", as.character(targets[1,1]), "\n")

  temp <- scan(file = as.character(targets[1,1]), skip = skip+1, quiet=TRUE,sep = sep, what = as.list(rep("character",  ncol(r))))


  samples = read.table(as.character(targets[1,2]), sep=",", header=T, skip=7)

  QC = readQC(as.character(targets[1,3]))

  #cat("First line", "\n")
  #print(r[1,])

  #calculate how many probes are on each array
  nrows = length(temp[[1]])
  
   #calculate how many arrays are in these files
  
#Find how many column names match the names we expect

   
#How many columns match the column heading for the Signal. ie how many arrays are there

#By default we expect each array to be listed in columns going across the file
readAcross = TRUE
#If only one column match was found then array are listed under each other in the file

avg.cols = which(match(strtrim(colnames(r), nchar(columns$AvgSig)), columns$AvgSig)==1)

matches = length(avg.cols)

if(length(which(match(strtrim(colnames(r), nchar(columns$ProbeID)), columns$ProbeID)==1))==0) stop(paste("Could not find any columns called", columns$ProbeID))



  #We assume that arrays are   
  if (matches >1) {
    ArraysPerFile = matches
  }

  else{
    readAcross=FALSE
    ArraysPerFile = nrow(r) / length(unique(r[,columns$ProbeID]))
  }




    for(i in 1:nrow(targets)){
      r = read.table(as.character(targets[i,1]), header = header, sep = sep, skip = skip)
       cat("Reading file ", as.character(targets[i,1]), "\n")   

      if(i==1){
        if(readAcross){
#First item in columns is the ProbeID columnd, we only read this once 
        a = readCols(r, columns=columns[2:5])
        

    #	Other columns
	if(!is.null(other.columns)) {

        BSData$other = readCols(r, columns=other.columns, colnames=other.columns) 
	}

          
      }

   
     else{
       for(j in 1:ArraysPerFile){
         if(j==1){

           
         
         r= r[(((j-1)*nrows)+1):(j*nrows),]

         BSData = readCols(r, columns=columns[2:5])
         
         if(!is.null(other.columns)) {

           BSData$other = readCols(r, columns=other.columns, colnames=other.columns) 
         }

       }

         else{
          r= r[(((j-1)*nrows)+1):(j*nrows),]

         BSData = cbind(BSData, readCols(r, columns=columns[2:5]))
         
         if(!is.null(other.columns)) {

           BSData$other = cbind(BSData$other,readCols(r, columns=other.columns, colnames=other.columns)) 
         }
           
       }
         
     }
     }
    }   
    else{
      
      if(readAcross){
      
      temp = readCols(r, columns=columns[2:5])
         
      
        if(!is.null(other.columns)) {
          temp$other = readCols(r, columns=other.columns, colnames=other.columns)
	}


    }
      else{

        
       for(j in 1:ArraysPerFile){

         
         r= r[(((j-1)*nrows)+1):(j*nrows),]

         temp= readCols(r, columns=columns[2:5])
         if(!is.null(other.columns)) {

           temp$other = readCols(r, columns=other.columns, colnames=other.columns) 
         }
       }
         
       }

      
    BSData = cbind(BSData, temp)
      
        
    }

    }


names(a$exprs) = names(a$BeadStDev) = names(a$NoBeads)=names(a$Detection) <- gsub("AVG_Signal.(\.+)","\\1",names(a$exprs))
rownames(a$exprs) = rownames(a$BeadStDev) = rownames(a$NoBeads) = rownames(a$Detection) = as.character(r[,1])

BSData = new("ExpressionSetIllumina")
assayData(BSData) = a
f = annotatedDataFrameFrom(a, byrow=TRUE)
featureData(BSData)=f

rownames(samples) = names(a$exprs)
p=new("AnnotatedDataFrame", samples,data.frame(labelDescription=colnames(samples), row.names=colnames(samples)))
phenoData(BSData) = p
BSData@QC = QC

BSData    

}



readCols <- function(table, columns = list(AvgSig = "AVG_Signal", NoBeads = "Avg_NBEADS",
                              Detection="Detection", BeadStDev="BEAD_STDEV"),
                     colnames = c("G", "NoBeads", "Detection", "BeadStDev")){
  
  cols.to.read=list(length=length(columns))

  nrows=nrow(table)
  
for(i in 1:length(columns)){

cols.to.read[[i]] = which(match(strtrim(colnames(table), nchar(columns[i])), columns[i])==1)

if(length(cols.to.read[[i]])==0) stop(paste("Could not find any columns called", columns[i]))

}



  #How many columns match the column heading for the Signal. ie how many arrays are there


  output=list()
  


  a = assayDataNew(exprs = table[,cols.to.read[[1]]], BeadStDev=table[,cols.to.read[[4]]], NoBeads=table[,cols.to.read[[2]]], Detection=table[,cols.to.read[[3]]], storage.mode="list")

  a

}


