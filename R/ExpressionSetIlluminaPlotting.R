


maplots <- function(object, SampleGroup = NULL,do.log=T){
  
  e<- exprs(object)
  
  if(do.log) e <- log2(exprs(object))
    
  if(is.null(SampleGroup)){
    
    message("No sample factor specified. Comparing to reference array")
    
    refdata <- rowMeans(e)
    
    plts <- alist()
    pCount <- 1
    
    #  for(i in 1:ncol(e)){
    
    #   refArray <- i
    #message("Computing M and A values with reference", colnames(e)[i])
    
    otherarrays <- e
    
    M <- lapply(1:ncol(otherarrays), function(x) refdata - otherarrays[,x])
    A <- lapply(1:ncol(otherarrays), function(x) 0.5*(refdata + otherarrays[,x]))
    
    mvals <- do.call("cbind", M)
    avals <- do.call("cbind", A)
    colnames(mvals) <- colnames(avals) <- colnames(otherarrays)
    
    
    df <- data.frame(melt(mvals),melt(avals))
    
    
#    plts[[1]] <- ggplot(df,aes(x=value.1,y=value))+
 #     stat_density2d(aes(alpha=..level..), geom="polygon") +
  #    scale_alpha_continuous(limits=c(0,0.2),breaks=seq(0,0.2,by=0.025))+
   #   geom_point(colour="steelblue",alpha=0.02)+ theme_bw()+geom_smooth(col="red",method="loess")+xlab("A") + ylab("M") + facet_wrap(~Var2) + theme(legend.position="none")
    
    plts <- ggplot(df, aes(x = value.1,y=value)) + stat_binhex(na.rm=T) + theme_bw()+xlab("A") + ylab("M") + facet_wrap(~Var2) + theme(legend.position="none") + ggtitle("Comparisons with Average array intensities")
    
  }
  
  else{
    
    if(is.null(SampleGroup)) stop("You must define a SampleGroup for the differential expression\n")
    
    if (SampleGroup %in% colnames(pData(object)))  esets <- split(sampleNames(object), pData(object)[,SampleGroup])
    else {
      print(paste(colnames(pData(object)),collapse=" "))
      stop("The SampleGroup argument must be a column name in the phenoData slot. See above for list of valid strings")
    }
    

    
    
    plts <- alist()
    
    for(i in 1:length(esets)){
      
      df <- alist()
      
      for(j in 1:length(esets[[i]])){
        
        refdata <- e[,esets[[i]][j]]
        otherarrays <- e[,esets[[i]][-j]]             
        
        M <- lapply(1:ncol(otherarrays), function(x) refdata - otherarrays[,x])
        A <- lapply(1:ncol(otherarrays), function(x) 0.5*(refdata + otherarrays[,x]))
        
        mvals <- do.call("cbind", M)
        avals <- do.call("cbind", A)
        colnames(mvals) <- colnames(avals) <- colnames(otherarrays)
        
        
        df[[j]] <- data.frame(melt(mvals),melt(avals), RefArray = esets[[i]][j])
        
      }
      
      df <- do.call("rbind",df)
    #  plts[[i]] <- #ggplot(df,aes(x=value.1,y=value))+
        #stat_density2d(aes(alpha=..level..), geom="polygon") +
        #scale_alpha_continuous(limits=c(0,0.2),breaks=seq(0,0.2,by=0.025))+
        #geom_point(colour="steelblue",alpha=0.02)+ theme_bw()+geom_smooth(col="red",method="loess")+xlab("A") + ylab("M") + facet_wrap(RefArray~Var2,ncol=length(esets[[i]])-1) + theme(legend.position="none")
      plts[[i]] <- ggplot(df, aes(x = value.1,y=value)) +
          stat_binhex(na.rm=T) + 
          theme_bw()+xlab("A") + 
          ylab("M") + 
          facet_wrap(RefArray~Var2) + 
          theme(legend.position="none") + 
          ggtitle(names(esets)[[i]])
      

      
    }
    names(plts) <- names(esets)

  }
  
  plts
}

setMethod("plotMA", signature(object="ExpressionSetIllumina"),maplots)

##DEPRECATED.....

plotMAXY <- function(exprs, arrays, log = TRUE, genesToLabel=NULL,labels=colnames(exprs)[arrays],labelCol="red", labelpch=16,foldLine=2,sampleSize=NULL,...){
  
  .Deprecated("maplots")
  
  mat <- matrix(c(0,1,0,0.04, 0,1,0.96,1, 0,0.04,0.04,0.96,
                  0.96,1,0.04,0.96, 0.04,0.96,0.04,0.96), byrow = T, ncol= 4)
  
  close.screen(all.screens = TRUE)
  
  split.screen(mat)
  
  split.screen(figs = c(length(arrays), length(arrays)), screen = 5)
  
  for(i in 1:length(arrays)){
    for(j in 1:length(arrays)){
      screen(((i-1)*length(arrays))+j+5)
      
      par(mar = c(0.3,0.3,0.3,0.3), cex.axis = 0.7)
      if(i == j){
        #        plot(0, col.axis = "white", cex = 0, col.lab = "white", tcl = -0, xlab = "", ylab = "")
        plot(0, axes = TRUE, type = "n", tcl = -0, col.axis = "white")
        text(1.0,0, labels = labels[i], cex=1)
      }
      else if(j < i){
        plotXY(exprs, array1 = arrays[i], array2 = arrays[j], xaxt = "n", yaxt = "n", log=log,genesToLabel=genesToLabel,foldLine=foldLine, labelCol=labelCol,labelpch=labelpch,sampleSize=sampleSize)
        if(i == length(arrays)){
          axis(1)
        }
        if(j == 1){
          axis(2)
        }
      }
      else{
        plotMA(exprs, array1 = arrays[i], array2 = arrays[j], xaxt = "n", yaxt = "n", log=log,genesToLabel=genesToLabel,foldLine=foldLine, labelCol=labelCol,labelpch=labelpch,sampleSize=sampleSize)
        if(i == 1){
          axis(3)
        }
        if(j == length(arrays)){
          axis(4)   
        }
      }
    }
  }
}


"plotMA" <- function(exprs, array1=1, array2=2, genesToLabel=NULL, labelCol="red", foldLine=2, log=TRUE,labelpch=16,ma.ylim=2,sampleSize=NULL,...){
  
  #try to catch any Limma objects and tell the user.
  if(class(exprs)[1] %in% c("RGList", "MAList")) 
    stop("\nIt appears you are trying to use the plotMA() function on a Limma object, but plotMA() is currently masked by beadarray\n\nIf you wish to use the Limma function, you can either call it directly using:\n\t\"limma::plotMA()\"\nor detach the beadarray package using:\n\t\"detach(package:beadarray)\"\n")
  
  exprs=as.matrix(exprs)
  
  if(log) exprs=log2(exprs)
  
  if(array2!=0 && array2!=array1){
    x = 0.5*(exprs[,array1] + exprs[,array2])
    y = exprs[,array1]- exprs[,array2]
  }
  else{
    #x = log2(exprs[,array1])
    #y = log2(exprs[,array1])
    stop("\'array1\' and \'array2\' must be different")
  }
  
  if(!is.null(sampleSize)){
    s = sample(1:length(x), sampleSize)
    x=x[s]
    y=y[s]
  }  
  
  naInd = intersect(which(!is.na(x)),which(!is.na(y)))
  naInd = intersect(naInd, which(!(is.infinite(x))))
  naInd = intersect(naInd, which(!(is.infinite(y))))
  
  smoothScatter(x[naInd],y[naInd], pch=16,cex=0.4, ylim=range(ma.ylim,-ma.ylim), xlab = "", ylab = "", ...) 
  
  abline(h=c(-log2(foldLine),0,log2(foldLine)),lty=c(2,1,2)) 
  
  if(!is.null(genesToLabel)){
    index = which(rownames(exprs) %in% genesToLabel)
    points(x[index], y[index], col=labelCol, pch=labelpch)
  }
}



"plotXY" <-
  function(exprs,array1=1, array2=2, genesToLabel=NULL, labelCol="red", log=TRUE,labelpch=16,foldLine=2,sampleSize=NULL,...){
    
    #XY plot of either two samples against each other, or red and green channels of one channel
    
    if(log) exprs=log2(exprs)
    
    exprs=as.matrix(exprs)
    
    if (array2!=0 && array2!=array1){
      
      
      x = exprs[,array1]
      y = exprs[,array2]
      
      
      xmax = 16
      xbox=18
      yspacing=0.3
      
      
      
      
    }
    else{
      #      x = log2(exprs[,array1])
      #      y = log2(exprs[,array1])
      #
      #      xmax = 16
      #      xbox=18
      #      yspacing=0.3
      stop("\'array1\' and \'array2\' must be different")
    }
    
    
    if(!is.null(sampleSize)){
      s = sample(1:length(x), sampleSize)
      x=x[s]
      y=y[s]
    }
    
    
    
    naInd = intersect(which(!is.na(x)),which(!is.na(y)))
    naInd = intersect(naInd, which(!(is.infinite(x))))
    naInd = intersect(naInd, which(!(is.infinite(y))))
    
    
    
    smoothScatter(x[naInd],y[naInd], xlim=range((max(0,min(x),na.rm=TRUE)),16), xlab = "", ylab = "", pch = 16, cex = 0.4, ...)
    abline(log2(foldLine), 1, lty=2)
    abline(-log2(foldLine),1,lty=2)
    
    if(!is.null(genesToLabel)){
      
      index = which(rownames(exprs) %in% genesToLabel)
      
      points(x[index], y[index], col=labelCol, pch=labelpch)
    }
    
  }
