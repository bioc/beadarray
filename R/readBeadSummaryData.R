readBeadSummaryData = function(dataFile, qcFile=NULL, sampleSheet=NULL, sep="\t", skip=8, ProbeID="ProbeID", columns = list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR", nObservations = "Avg_NBEADS", Detection="Detection Pval"), qc.sep="\t", qc.skip=8, controlID="ProbeID", qc.columns = list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", nObservations="Avg_NBEADS", Detection="Detection Pval"), annoPkg=NULL, dec=".", quote="", annoCols = c("TargetID", "PROBE_ID","SYMBOL"))
{

if(!(is.null(sampleSheet))){ 
samples = read.table(sampleSheet, sep=",", header=TRUE, skip=7, as.is=TRUE)
}

r = read.table(as.character(dataFile), sep=sep, header=TRUE, skip=skip, dec=dec, quote=quote, as.is=TRUE, row.names=NULL, check.names=FALSE, strip.white=TRUE, comment.char="", fill=TRUE)

#foundColumns = NULL

index = grep(ProbeID, colnames(r))

##Try to match the annotation columns

annoCols = match(annoCols, colnames(r))

annoCols = annoCols[!is.na(annoCols)]

annoMat = r[,annoCols]


if(length(index)!=0){

ProbeID = r[,index]
}

else {
  stop("Could not find a column called ", ProbeID, " to use as bead identifiers.  Check your file and try changing the \'ProbeID\' and/or \'skip\' arguments.")
}

#check for non-unique probe names
if(length(ProbeID) != length(unique(ProbeID))){
      notdup = !duplicated(ProbeID)
      warning("ProbeIDs non-unique: consider setting 'ProbeID' to another column containing unique identifiers. ",   sum(!notdup), " repeated entries have been removed.\n")
      ProbeID = ProbeID[notdup]
      r = r[notdup,]

      if(!is.na(annoCols)) annoMat = r[notdup,annoCols]

#    warning("ProbeIDs non-unique: consider setting 'ProbeID' to another column containing unique values.  Adding extension '.repX' to Xth replicate ID to enforce uniqueness.\n")
#    dups = unique(ProbeID[duplicated(ProbeID)])
#    for(j in 1:length(dups)) {
#      sel = ProbeID==dups[j]
#      ProbeID[sel] = paste(ProbeID[sel], ".rep", seq(1:sum(sel)), sep="")
#    }
}

data = index = list()
ncols = NULL
nrows = nrow(r)

for(i in 1:length(columns)){

  index[[i]] = grep(columns[[i]], colnames(r))
  ncols[i] = length(index[[i]])
  if(ncols[i] == 0){
    cat("Could not find a column called: ", columns[[i]], "\n")
  }
}

if(sum(ncols)==0) {
  stop("No data found, check your file or try changing the \'skip\' argument")
}

i = seq(1:length(ncols))[ncols==max(ncols)][1]
defColNames = sub(paste("(.|)", columns[[i]], "(.|)", sep=""), "", colnames(r)[index[[i]]])

for(i in 1:length(columns)){
  if(ncols[i]==max(ncols)) {
    data[[i]] = r[,index[[i]]]
    colNames = sub(paste("(.|)", columns[[i]], "(.|)", sep=""), "", colnames(r)[index[[i]]])
    dupColNames = unique(colNames[duplicated(colNames)])
    if(length(dupColNames)!=0) {
      for(j in 1:length(dupColNames)) {
        sel = colNames==dupColNames[j]
        colNames[sel] = paste(colNames[sel], ".rep", seq(1:sum(sel)), sep="")
      }
    }
    colnames(data[[i]]) = colNames
  }
  else { # data missing
    data[[i]] = matrix(NA, nrows, max(ncols))
    colnames(data[[i]]) = defColNames
  }
  rownames(data[[i]]) = ProbeID
#  if(!(is.null(sampleSheet))){
#    colnames(data[[i]]) = as.character(samples[,4])
#  }
}

names(data) = names(columns) #foundColumns

BSData = new("ExpressionSetIllumina")

if(!is.null(annoPkg) && is.character(annoPkg))
  BSData@annotation = annoPkg

for(i in 1:length(data)){
  index = which(names(assayData(BSData))== names(data)[i])

  if(ncols[i]==0) {
    cat("Missing data - NAs stored in slot", names(data)[i], "\n")
  }

  assayData(BSData)[[index]] = as.matrix(data[[i]])
}

if(!(is.null(qcFile))){
  QC = readQC(file=qcFile, sep=qc.sep, skip=qc.skip, columns=qc.columns, controlID=controlID, dec=dec, quote=quote)
  if(ncol(QC$exprs)!=ncol(exprs(BSData))) {
     warning("Number of arrays doesn't agree: ", ncol(exprs(BSData)), " in dataFile, versus ", ncol(QC$exprs), " in qcFile.  qcFile ignored.")
  }
  else {
    reorder = match(colnames(QC[[i]]), colnames(exprs(BSData)))
    notagree = colnames(QC$exprs)!=colnames(exprs(BSData))
    if(sum(notagree)==0) {
       	
	for(i in 1:length(QC)){
		slotMatch = match(names(QC)[i], names(assayData(BSData)))	
		##Only bind the QC column to the existing assayData if it contains data
		if(!is.na(slotMatch)){
			if(ncol(QC[[i]]) == ncol(assayData(BSData)[[slotMatch]])){
			assayData(BSData)[[slotMatch]] = rbind(assayData(BSData)[[slotMatch]], QC[[i]])
			
			dupIDs = which(duplicated(rownames(assayData(BSData)[[slotMatch]])))
			
			if(length(dupIDs) > 0){
				rownames(assayData(BSData)[[slotMatch]])[dupIDs] = paste(rownames(assayData(BSData)[[slotMatch]])[dupIDs], ".2", sep="")
			}

			}
		}
	
	}
    }
    else {
      if(length(reorder)!=0) {
          for(i in 1:length(BSData@QC)) {
            reorder = sapply(colnames(QC[[i]]), FUN="grep",  colnames(exprs(BSData)), fixed=TRUE)
            if(length(reorder)>0) {
              QC[[i]] = QC[[i]][, reorder]
		##Only bind the QC column to the existing assayData if it contains data
		slotMatch = match(names(QC)[i], names(assayData(BSData))[i])
		if(!is.na(slotMatch)){			
			if(ncol(QC[[i]]) == ncol(assayData(BSData)[[slotMatch]])){
			assayData(BSData)[[slotMatch]] = rbind(assayData(BSData)[[slotMatch]], QC[[i]])

			dupIDs = which(duplicated(rownames(assayData(BSData)[[slotMatch]])))
			
			if(length(dupIDs) > 0){
				rownames(assayData(BSData)[[slotMatch]])[dupIDs] = paste(rownames(assayData(BSData)[[slotMatch]])[dupIDs], ".2", sep="")
			}

			}
		}
            }
          }
      }
      else {
        warning("Could not match array names used in dataFile with those in qcFile.  qcFile ignored.")
      }
    }
  }
}

else QC = NULL

if(!(is.null(sampleSheet))){
  colmatch = grep(colnames(exprs(BSData)), samples, fixed=TRUE)
  ord = match(colnames(exprs(BSData)), samples[,colmatch])
  if(length(colmatch)==1 && sum(is.na(ord))==0) {
    samples = samples[ord,]
    rownames(samples) = colnames(exprs(BSData))
    p = new("AnnotatedDataFrame", samples, data.frame(labelDescription=colnames(samples), row.names=colnames(samples)))
  }
  else { 
    warning("Could not reconcile dataFile with sampleSheet information. sampleSheet ignored.")
    p = new("AnnotatedDataFrame", data.frame(sampleID=colnames(exprs(BSData)),row.names=colnames(exprs(BSData))))
  }
}

else {
  p = new("AnnotatedDataFrame", data.frame(sampleID=colnames(exprs(BSData)),row.names=colnames(exprs(BSData))))
}


##Have we added any more IDs to BSData

Status = rep("Gene", length(ProbeID))

if(!is.null(QC)){
	Type = QC$Type
	newIDs = rownames(QC[[1]])

	##Some IDs might already be in the QC and the Gene data


	ProbeID = c(ProbeID, newIDs)

	dupIDs = which(duplicated(ProbeID))


	newAnno = matrix(nrow = length(newIDs), ncol=ncol(annoMat), NA)
	colnames(newAnno) = colnames(annoMat)	


	annoMat = rbind(annoMat, newAnno)
	Status = c(Status, Type)	


	##Some 
	if(length(dupIDs) > 0){ 
		ProbeID[dupIDs]	= paste(ProbeID[dupIDs], ".2", sep="")	
	}	

}

featureData = new("AnnotatedDataFrame", data=data.frame(ProbeID,annoMat, row.names=ProbeID, Status=Status))

phenoData(BSData) = p
featureData(BSData) = featureData

BSData

}

readQC=function(file, sep="\t", skip=8, controlID = "ProbeID", columns=list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", nObservations="Avg_NBEADS", Detection="Detection Pval"),  dec=".", quote=""){

##
  
  r=read.table(as.character(file), sep=sep, header=TRUE, skip=skip, quote=quote, as.is=TRUE, check.names=FALSE, strip.white=TRUE, comment.char="", fill=TRUE)


##If there is an ArrayID column, the QC file is BeadStudio version 1 and each row in the file is an array
#  read.table(file = fileName, header = TRUE, sep = sep, 
#        skip = nMetaDataLines, row.names = NULL, quote = "", 
#        as.is = TRUE, check.names = FALSE, strip.white = TRUE, 
#        comment.char = "", fill = TRUE)

 ArrayID = grep("ArrayID", colnames(r))

type = NULL

 typeCol = grep("TargetID", colnames(r))

 if (!is.na(typeCol)){
	type = r[,typeCol]
}

 if(length(ArrayID) != 0){
    BeadStudioVersion=1
    index = grep(columns$exprs, colnames(r))
    if(length(index)!=0){
       ProbeID = sub(paste("(.|)", columns$exprs, "(.|)", sep=""), "", colnames(r))[index]
       nrows = length(ProbeID)
    }
    else { # can't find controlIDs - report warning
       stop("Could not find a column called ", controlID, " to use as control identifiers.  Check your file and try changing the \'controlID\' and/or \'skip\' arguments.")
    }
}

 else {
    BeadStudioVersion=2 
    nrows = nrow(r)
    index = grep(controlID, colnames(r))
  
    if(length(index)!=0){
       ProbeID = r[,index]
    }
    else { # can't find controlIDs - report warning
       stop("Could not find a column called ", controlID, " to use as control identifiers.  Check your file and try changing the \'controlID\' and/or \'skip\' arguments.")
    }
}

# check for non-unique probe names
if(length(ProbeID) != length(unique(ProbeID))){
  notdup = !duplicated(ProbeID)
  warning("controlIDs non-unique: ", sum(!notdup), " repeated entries have been removed.\n")
  ProbeID = ProbeID[notdup]
  r = r[notdup,]
  type = type[notdup]	
}

data = index = list()
ncols = NULL

for(i in 1:length(columns)){
  index[[i]] = grep(columns[[i]], colnames(r))
  ncols[i] = length(index[[i]])
  if(ncols[i] == 0){
    cat("[readQC] Could not find a column called: ", columns[[i]], "\n")
  }
}

if(sum(ncols)==0) {
  stop("No data found, check your file or try changing the \'skip\' argument")
}

i = seq(1:length(ncols))[ncols==max(ncols)][1]
defColNames = sub(paste("(.|)", columns[[i]], "(.|)", sep=""), "", colnames(r)[index[[i]]])

for(i in 1:length(columns)){
  if(ncols[i]==max(ncols)) {
    if(BeadStudioVersion == 2) { 
      data[[i]] = as.matrix(r[,index[[i]]])
      colnames(data[[i]]) = sub(paste("(.|)", columns[[i]], "(.|)", sep=""), "", colnames(r)[index[[i]]])
    }
    else {
      data[[i]] = as.matrix(t(r[,index[[i]]]))
      colnames(data[[i]]) = r[,ArrayID]
    }
  }
  else { # data missing
    if(BeadStudioVersion == 2) { 
      data[[i]] = matrix(NA, nrows, max(ncols))
      colnames(data[[i]]) = defColNames
    }
    else {
      data[[i]] = matrix(NA, nrows, nrow(r))
      colnames(data[[i]]) = r[,ArrayID]
    }
  }
  rownames(data[[i]]) = ProbeID
}

names(data) = names(columns) # foundColumns
QC = list(exprs=new("matrix"), se.exprs=new("matrix"), Detection=new("matrix"), nObservations=new("matrix"),Type = new("character"))


for(i in 1:length(data)){
  index = which(names(QC)== names(data)[i])

  if(ncols[i]==0) {
    cat("[readQC] Missing data - NAs stored in slot", names(data)[i], "\n")
  }

  QC[[index]] = data[[i]]

}

if(!is.null(type)){

	QC$Type = type

}
 QC

}
