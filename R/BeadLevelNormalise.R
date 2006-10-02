
BeadLevelNormalise = function(targets, imageManipulation="sharpen", backgroundCorrect="none", probes=NULL,imagesPerArray=2){

  
print("Calculating Target distribution...")

if(imagesPerArray==2){

  
BLData = readBeadLevelData(targets[1:2,], imageManipulation=imageManipulation)

BLData = backgroundCorrect(BLData, method=backgroundCorrect)

G = rbind(BLData@G[,1], BLData@G[,2])

plot(density(log2(G),na.rm=TRUE))

o = order(BLData@G[,1], decreasing=TRUE)

T = G[o]

arraysToRead = c(1,2)

}

else if (imagesPerArray==1){

BLData = readBeadLevelData(targets[1,], imageManipulation=imageManipulation)
  
BLData = backgroundCorrect(BLData, method=backgroundCorrect)

o = order(BLData@G[,1], decreasing=TRUE)

T = BLData@G[o,1]

plot(density(log2(BLData@G[,1]),na.rm=TRUE))
}


narrays = nrow(targets) / imagesPerArray


for(i in 2:narrays){

if(imagesPerArray==1){

BLData = readBeadLevelData(targets[i,], imageManipulation=imageManipulation)

BLData = backgroundCorrect(BLData, method=backgroundCorrect)

o = order(BLData@G[,1], decreasing=TRUE)

lines(density(log2(BLData@G[,1]),na.rm=TRUE))

T = T + BLData@G[o,1]

}

else if(imagesPerArray==2){

arraysToRead = arraysToRead+2  

BLData = readBeadLevelData(targets[arraysToRead,], imageManipulation=imageManipulation)  
  
BLData = backgroundCorrect(BLData, method=backgroundCorrect)

G = rbind(BLData@G[,1], BLData@G[,2])

o = order(G, decreasing=TRUE)

lines(density(log2(G),na.rm=TRUE))

T = T + G[o]

rm(G)

}

rm(BLData)

gc()

}



T = T / narrays

lines(density(log2(T),na.rm=TRUE), col="red")

print("Now reading and normalising...")

if(imagesPerArray==1){

  arraysToRead=1
}
else{
  arraysToRead =c(1,2)

}


BLData = readBeadLevelData(targets[arraysToRead,], imageManipulation=imageManipulation)

BLData = backgroundCorrect(BLData, method=backgroundCorrect)

BLData@G = normalize.qspline(as.matrix(BLData@G), target=T)

BSData = createBeadSummaryData(BLData, imagesPerArray=imagesPerArray, probes=probes)


rm(BLData)

gc()





for(i in 2:(narrays)){

  if(imagesPerArray==1){
    arraysToRead = i
  }
  else{
    arraysToRead = arraysToRead + 2
  }

  

BLData = readBeadLevelData(targets[arraysToRead,], imageManipulation=imageManipulation)

BLData = backgroundCorrect(BLData, method=backgroundCorrect)

BLData@G = normalize.qspline(as.matrix(BLData@G), target=T)


BSData2 = createBeadSummaryData(BLData, imagesPerArray=imagesPerArray, probes=probes)

e = new("ExpressionSetIllumina")

assayData(e) = combine(assayData(BSData), assayData(BSData2))

BSData = e


rm(BLData)

gc()

}



BSData

}
