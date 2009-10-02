boxplotBeads = function(BLData, whatToPlot="G", arrays=NULL,
                                  log=TRUE, varwidth=TRUE, method="illumina",
                                  n = 3, trim=0.05,...) {
  tmp = list()
  arraynms = arrayNames(BLData)
  narrays = length(arraynms)
  if(is.null(arrays))  # plot all arrays
    arrays = 1:narrays
  for(i in arrays)
      tmp[[arraynms[i]]] = getArrayData(BLData, array=i, what=whatToPlot, log=log, method=method, n=n, trim=trim)
  boxplot(tmp,varwidth=varwidth,...)
}

plotRG = function(BLData, ProbeIDs=NULL, BeadIDs=NULL, log=TRUE, arrays=1,
                   xlim=c(8,16), ylim=c(8,16), xlab="G intensities",
                   ylab="R intensities", main=arrayNames(BLData)[arrays],
                   smooth=TRUE, cols=NULL, ...) {
  arraynms = arrayNames(BLData)
  narrays = length(arrays)
  if(length(arrays)==1 && is.null(ProbeIDs) && is.null(BeadIDs)) {
     if(smooth)
       smoothScatter(getArrayData(BLData, what="G", array=arrays, log=log),
            getArrayData(BLData, what="R", array=arrays, log=log),
            xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
            main=main, ...)
     else
       plot(getArrayData(BLData, what="G", array=arrays, log=log),
            getArrayData(BLData, what="R", array=arrays, log=log),
            xlim=xlim, ylim=ylim, new=TRUE, xlab=xlab, ylab=ylab,
            main=main, ...)
  }
  else{
    if(is.null(cols))
      cols = rainbow(narrays)
    plot(1, xlim=xlim, ylim=ylim, type="n", new=TRUE, xlab=xlab, ylab=ylab, main=main, ...)
     for(i in arrays) {
           if(!is.null(ProbeIDs))
             for(j in 1:length(ProbeIDs)) {
               points(getArrayData(BLData, what="G", array=i, log=log)[which(BLData[[i]]$ProbeID %in% ProbeIDs[j])],
                      getArrayData(BLData, what="R", array=i, log=log)[which(BLData[[i]]$ProbeID %in% ProbeIDs[j])], col = cols[which(arrays %in% i)], ...)
                }
            if (!is.null(BeadIDs)) 
                points(getArrayData(BLData, what = "G", array = i, 
                  log = log)[BeadIDs], getArrayData(BLData, what = "R", 
                  array = i, log = log)[BeadIDs], col = cols[which(arrays %in% i)],
                  ...)
        }
    }
}

plotBeadLocations = function(BLData, ProbeIDs=NULL, BeadIDs=NULL, array=1, SAM=FALSE, xlab="x-coordinate", ylab="y-coordinate", main=paste("Bead", ProbeIDs, "locations"),...){

xmax = max(BLData[[array]]$GrnX)
ymax = max(BLData[[array]]$GrnY)

plot(1, xlim=range(0:xmax), ylim=range(0:ymax) ,type="n", new=TRUE, xlab=xlab, ylab=ylab,main=main,...)

if(!(is.null(ProbeIDs))){

xs = BLData[[array]]$GrnX[which(BLData[[array]]$ProbeID %in% ProbeIDs)]
ys = ymax-BLData[[array]]$GrnY[which(BLData[[array]]$ProbeID %in% ProbeIDs)]

}
else{

xs = BLData[[array]]$GrnX[BeadIDs]
ys = ymax-BLData[[array]]$GrnY[BeadIDs]
 
}

if(SAM) {

yy = c(ymax/2, 0, 0, ymax/2, ymax, ymax)

xx = c(0,xmax/4, 0.75*xmax, xmax, 0.75*xmax, xmax/4)

polygon(xx, yy)

}

else polygon(x=c(0,0,xmax,xmax), y=c(0, ymax, ymax,0))


points(xs, ys,...)


}


plotBeadIntensities = function(BLData, ProbeIDs, arrays, log=FALSE, whatToPlot="G", ProbeCols=NULL, ylim=NULL,...){

nplots  = length(ProbeIDs)*length(arrays)

if(is.null(ProbeCols)){

ProbeCols = rainbow(n=length(ProbeIDs), start=0, end=5/6)

}

if(is.null(ylim)) {
  if(log==TRUE)
    ylim=c(0,16)
  else
    ylim=c(1,2^16)
}

plot(4, xlim=range(0,nplots), ylim=ylim, type="n", axes=FALSE, xlab="", ylab="intensities",...)
count=1


for(i in 1:length(arrays)){
for(j in 1:length(ProbeIDs)){
intens = getProbeIntensities(BLData, ProbeIDs=ProbeIDs[j], array=arrays[i], log=log, what=whatToPlot)
  if (length(intens) >= 1) {
      sel = is.finite(intens) & !is.na(intens)
      boxplot(intens[sel], at = count - 0.5, add = TRUE, axes = FALSE, col = ProbeCols[j])
  }
  else {
    cat("\nNo bead with ID", ProbeIDs[j], "on array", i, "\n")
  }
#points(x=rep(count-0.5, length(o)), y = o, pch=16, col=ProbeCols[j])

count = count+1

}

if(i!=length(arrays)) abline(v=count-1)

}
axis(2)
box()

}




imageplot = function(BLData, array = 1, nrow = 100, ncol = 100,
                     low = NULL, high = NULL, ncolors = 123,
                     whatToPlot ="G", log=TRUE, zlim=NULL,
                     main=whatToPlot, method="illumina",
                     n = 3, trim=0.05, legend=TRUE, SAM=FALSE, ...){
    
    if(class(BLData)[1] %in% c("RGList", "MAList")) 
        stop("\nIt appears you are trying to use the imageplot() function on a Limma object, but imageplot() is currently masked by beadarray\n\nIf you wish to use the Limma function, you can either call it directly using:\n\t\"limma::imageplot()\"\nor detach the beadarray package using:\n\t\"detach(package:beadarray)\"\n")

  par(mar = c(2,1,1,1), xaxs = "i", yaxs = "i")
  
#Not needed since the co-ords are automatically scaled to zero now  
#  xs = floor(BLData@GrnX[,array] - min(BLData@GrnX[,array]))
#  ys = floor(BLData@GrnY[,array] - min(BLData@GrnY[,array]))

  whatToPlot = match.arg(whatToPlot, choices=c("G", "Gb", "R", "Rb", "wtsG", "wtsR", "residG", "residR", "M", "residM", "A", "beta"))
  if((whatToPlot=="R" | whatToPlot=="residR" | whatToPlot=="M" | whatToPlot=="residM" | whatToPlot=="A" | whatToPlot=="beta") & BLData@arrayInfo$channels!="two")
    stop(paste("Need two-channel data to plot", whatToPlot, "values"))
                                          
  data = getArrayData(BLData, what=whatToPlot, array=array, log=log, method=method, n=n, trim=trim) 
  ind = is.na(data) | is.infinite(data)
  if(sum(ind)>0) {
    cat(paste("Warning:", sum(ind), "NA, NaN or Inf values, which will be ignored.\nCheck your data or try setting log=\"FALSE\"\n"))
    data = data[!ind] # dat[ind]= 0
    #rm(ind)
  }

  if (is.character(low)) 
    low = col2rgb(low)/255
  if (is.character(high)) 
    high = col2rgb(high)/255
  if (!is.null(low) && is.null(high)) 
    high = c(1, 1, 1) - low
  if (is.null(low) && !is.null(high)) 
    low = c(1, 1, 1) - high

  if (is.null(low)) 
    low = c(1, 1, 1)
  if (is.null(high)) 
    high = c(0, 0, 1)
#  if(whatToPlot=="G" | whatToPlot=="Gb")
    col = rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
          high[2], len = ncolors), seq(low[3], high[3], len = ncolors))

#  else  # plot in Red colour scheme
#    col = rgb(seq(low[2], high[2], len = ncolors), seq(low[1], 
#          high[1], len = ncolors), seq(low[3], high[3], len = ncolors))
  if(SAM) {
    xs = floor(BLData[[array]]$GrnX[!ind])
    ys = floor(BLData[[array]]$GrnY[!ind])
    rm(ind)
  }
  else { # for BeadChip - switch X and Y
    xs = floor(BLData[[array]]$GrnY[!ind])
    ys = floor(BLData[[array]]$GrnX[!ind])
  }
  xgrid = floor(seq(0, max(xs), by = max(xs)/ncol))
  ygrid = floor(seq(0, max(ys), by = max(ys)/nrow))

  imageMatrix = matrix(ncol = ncol, nrow = nrow)
  zr = NULL
  for(i in 1:ncol){
    idx = which((xs > xgrid[i]) & (xs < xgrid[i+1]))
    if(length(idx)>0) {
      fground = data[idx]
      yvalues = ys[idx]
#      yvalues = BLData@GrnY[idx,array]

      out = .C("BLImagePlot", length(fground), as.double(fground), as.double(yvalues), as.integer(ygrid),
              result = double(length = nrow), as.integer(nrow), PACKAGE = "beadarray")
      zr[1] = min(zr[1], out$result[!is.na(out$result)], na.rm=TRUE)
      zr[2] = max(zr[2], out$result[!is.na(out$result)], na.rm=TRUE)
      if(!is.null(zlim)) {
         out$result[!is.na(out$result)] = pmax(zlim[1], out$result[!is.na(out$result)], na.rm=TRUE)
         out$result[!is.na(out$result)] = pmin(zlim[2], out$result[!is.na(out$result)], na.rm=TRUE)
      }
      if(SAM)
         imageMatrix[,i] = rev(out$result)
      else
         imageMatrix[,i] = out$result
    }
  }
#  if(!is.null(zlim)) {
#     imageMatrix = apply(imageMatrix, 1, FUN="pmin", zlim[2], na.rm=TRUE)
#     imageMatrix = apply(imageMatrix, 1, FUN="pmax", zlim[1], na.rm=TRUE)
#  }
#  zr =  range(imageMatrix, na.rm=TRUE)

  imageMatrix = t((imageMatrix))

  if(is.null(zlim)) zlim=range(imageMatrix, na.rm=TRUE)
  image(x = c(0:ncol), y = c(0:nrow), z = imageMatrix,  xaxt = "n", yaxt = "n", col = col, main=main, zlim=zlim,...)

  if(legend)
    mtext(paste("z-range ",round(zr[1],1)," to ",round(zr[2],1)," (saturation ",round(zlim[1],1),", ",round(zlim[2],1),")",sep=""),side=1,cex=0.6)
}

SAMSummary =
function(BLData, mode="outliers", whatToPlot="G", samID=NULL, log=TRUE, n=3, colour=TRUE, scale = NULL, low="yellow", high="red",...) {

if(is.null(samID))
  samID=BLData@arrayInfo$chip[1]

split.screen(c(1,2))
screen(1)

sel = BLData@arrayInfo$chip==samID

if(sum(sel)>96)
   stop("SAMSummary cannot plot more than 96 arrays - please check your data")

arraynms = arrayNames(BLData)[sel]
rows = BLData@arrayInfo$row[sel]
cols = BLData@arrayInfo$col[sel]

roword = c("R001", "R002", "R003", "R004", "R005", "R006", "R007", "R008")
colord = c("C001", "C002", "C003", "C004", "C005", "C006", "C007", "C008", "C009", "C010", "C011", "C012")

len = 96
allnames = rep(NA, len)
values = vector(length=len)
arrayseq = 1:len

if(mode == "outliers"){

  o = list(length=len)

  for(i in 1:8) {
     for(j in 1:12){
        tmparray = arraynms[rows==roword[i] & cols==colord[j]]
        tmpind = (i-1)*12+j
        if(length(tmparray)>0) {
           allnames[tmpind] = tmparray
           o[[tmpind]] = findAllOutliers(BLData, array=tmparray, log=log, what=whatToPlot, n=n,...)
           values[tmpind] = length(o[[tmpind]])
        }
        else {
          o[[tmpind]] = NA
          values[tmpind] = NA
        }
     }
  }
}

if(mode == "intensities") {
  for(i in 1:8) {
     for(j in 1:12){
        tmparray = arraynms[rows==roword[i] & cols==colord[j]]
        tmpind = (i-1)*12+j
        if(length(tmparray)>0) {
           allnames[tmpind] = tmparray
           values[tmpind] =  median(getArrayData(BLData, array=tmparray, log=log, what=whatToPlot, ...), na.rm=TRUE)
        }
        else
          values[tmpind] = NA
      }
   }
  values = values-min(values, na.rm=TRUE)
}

doLoop = TRUE
counter = 0
while(doLoop){
counter = counter+1
x=1/13
y=1/9

ys = c(y/2, 0, 0, y/2, y, y)

for(i in 1:8){

xs = c(0,x/4, 0.75*x, x, 0.75*x, x/4)

if(is.null(scale)) scale=max(values, na.rm=TRUE)

for(j in 1:12){

array_index =  (12 * (8-i)) + j

control_intensity = values[array_index]

if(is.na(control_intensity)){
polygon(xs,ys, col="white") # rgb(1,1,1))
}
else{

if(colour){ polygon(xs,ys, col=rgb(control_intensity/scale,0,0))}
else{ polygon(xs,ys,col=gray(control_intensity/scale))}

}


#text(xs[3], ys[4], BLData@other@SAMPLE[1,array_index])

xs = xs + 1/12

}

ys = ys + 1/8


}
if(counter>1) {
  answer = readline("Type Q to quit or anything else to continue, followed by <ENTER>:  ")
  if(answer=="Q")
    break# doLoop=FALSE
}
cat(paste("\nClick on plot device to select an array to show", mode, "from\n\n"))
l=locator(n=1)

x.clicked = strtrim(as.character(l$x),5)
#print(l)

y.clicked = strtrim(as.character(l$y),5)

y.clicked = 8-(as.double(y.clicked)*8)

x.clicked = as.double(x.clicked)*12

#print(x.clicked)
#print(y.clicked)

ArrayClickedOn = ceiling(y.clicked)*12 + ceiling(x.clicked) - 12

if(!is.na(allnames[ArrayClickedOn])) { #print(ArrayClickedOn)
  cat(paste("Displaying", mode, "from array", allnames[ArrayClickedOn], "\n\n"))

  screen(2, new=TRUE)

  plot(1:10, type="n", axes=FALSE, xlab="", ylab="")

  if(mode == "outliers"){


    os= o[[ArrayClickedOn]]

    plotBeadLocations(BLData, BeadIDs=os, array=allnames[ArrayClickedOn], SAM=TRUE, new=TRUE, main=paste(allnames[ArrayClickedOn], "outliers"))
    screen(1)
    }

  if(mode == "intensities") {

    imageplot(BLData, array=allnames[ArrayClickedOn], whatToPlot=whatToPlot, log=log, high=high, low=low,main=paste(allnames[ArrayClickedOn], whatToPlot), SAM=TRUE, ...)
    screen(1)
  }
}
else 
  cat("No data to plot, try again\n\n")
  
}
}

BeadChipSummary =
function(BLData, mode="outliers", whatToPlot="G", chipID=NULL, stripsPerChip=12, log=TRUE, n=3, colour=TRUE, scale = NULL,low="yellow", high="red", ...){

if(is.null(chipID))
  chipID = BLData@arrayInfo$chip[1]

split.screen(c(1,2))
screen(1)

#if(mode=="outliers") title=paste("Number of Outliers:", whatToPlot)
#if(mode=="intensities") title=paste("Median Intensities:", whatToPlot)

sel = BLData@arrayInfo$chip==chipID

if(sum(sel)>stripsPerChip)
   stop(paste("BeadChipSummary cannot plot more than", stripsPerChip, "arrays - please check your data"))

arraynms = arrayNames(BLData)[sel]
sub = BLData@arrayInfo$row[sel]
strp = BLData@arrayInfo$col[sel]

subarray = c("A", "B", "C", "D", "E", "F")
if(stripsPerChip==12)
  strip = 1:2
else
  strip = 1

#arraynms = arrayNames(BLData)[sel]

#narrays = length(arraynms)
len = stripsPerChip # min( , narrays)
allnames = rep(NA, len)
values = vector(length=len)

if(mode == "outliers"){

  o = list(length=len)

  for(i in 1:length(subarray)) {
     for(j in strip){
        tmparray = arraynms[sub==subarray[i] & strp==j]
        tmpind = (i-1)*2+j
        if(length(tmparray)>0) {
           allnames[tmpind] = tmparray
           o[[tmpind]] = findAllOutliers(BLData, array=tmparray, log=log, what=whatToPlot, n=n,...)
           values[tmpind] = length(o[[tmpind]])
        }
        else {
          o[[tmpind]] = NA
          values[tmpind] = NA
         }
     }
  } 
}
if(mode == "intensities"){
 for(i in 1:length(subarray)) {
    for(j in strip){
        tmparray = arraynms[sub==subarray[i] & strp==j]
        tmpind = (i-1)*2+j
        if(length(tmparray)>0) {
           allnames[tmpind] = tmparray
           values[tmpind] =  median(getArrayData(BLData, array=tmparray, log=log, what=whatToPlot,...), na.rm=TRUE)
        }
        else
          values[tmpind] = NA
    }
  }
  values = values-min(values, na.rm=TRUE)
}
tol=0.015
doLoop = TRUE
counter = 0
while(doLoop){
counter = counter + 1
  y=1
ys = c(y-tol, y-tol, y+tol-y/stripsPerChip, y+tol-y/stripsPerChip)

  

for(i in 1:len){
control_intensity = values[i]


xs=c(0.1,1,1,0.1)
if(is.null(scale)) scale= max(values, na.rm=TRUE)
if(is.na(control_intensity)){
polygon(xs,ys, col="white") # rgb(1,1,1))
}
else{

if(colour){ polygon(xs,ys, col=rgb(control_intensity/scale,0,0))}
else{ polygon(xs,ys,col=gray(control_intensity/scale))}

}

ys = ys - 1/stripsPerChip

}

if(counter>1) {
  answer = readline("Type Q to quit or anything else to continue, followed by <ENTER>:  ")
  if(answer=="Q")
    break# doLoop=FALSE
}
cat(paste("\nClick on plot device to select an array to show", mode, "from\n\n"))
l=locator(n=1)


y.clicked = strtrim(as.character(l$y),5)

y.clicked = as.double(y.clicked)*stripsPerChip


ArrayClickedOn = stripsPerChip - ceiling(y.clicked) + 1
if(!is.na(allnames[ArrayClickedOn])) {
  cat(paste("Displaying", mode, "from array", arraynms[ArrayClickedOn], "\n\n"))
  screen(2, new=TRUE)

  plot(1:10, type="n", axes=FALSE, xlab="", ylab="")

  if(mode == "outliers"){

    os= o[[ArrayClickedOn]]

    plotBeadLocations(BLData, BeadIDs=os, array=allnames[ArrayClickedOn], new=TRUE, main=paste(allnames[ArrayClickedOn], "outliers"))

    screen(1)
  }

  if(mode == "intensities"){

    imageplot(BLData, array=allnames[ArrayClickedOn], whatToPlot=whatToPlot, log=log, high=high, low=low, main=paste(allnames[ArrayClickedOn], whatToPlot),...)

    screen(1)
  }
}

else 
  cat("No data to plot, try again\n\n")
}

}

####### taken from the sma package #####
plot.smooth.line  <- function(x, M, f = 0.1, ...)
{
#  A <- x
  ind <- !(is.na(x) | is.na(M) | is.infinite(x) | is.infinite(M))
  #lines(lowess(A[ind], M[ind], f = f), ...)
  lines(approx(lowess(x[ind], M[ind], f = f)), ...)  
}