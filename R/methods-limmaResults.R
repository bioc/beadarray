setMethod("initialize", "limmaResults",
          function(.Object,
                   assayData = assayDataNew(LogFC=LogFC,LogOdds=LogOdds, PValue=PValue, storage.mode="list"),
                   
                   phenoData = new("AnnotatedDataFrame"),
                   LogFC=new("matrix"),
                   LogOdds=new("matrix"),
                   PValue=new("matrix"),
                   annotation = character(),
                   featureData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME")
          )
{
            .Object<-callNextMethod(.Object,
                                    LogFC = LogFC,
                                    LogOdds = LogOdds,
                                    PValue = PValue,  
                                    experimentData = experimentData,
                                    annotation = annotation,
                                    featureData = featureData
            )
            
            .Object
          })


setAs("limmaResults", "GRanges",
      function(from)
      {
        annoName <- annotation(from)
        
        annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)
        
        if(annoLoaded){
          
          
          mapEnv <-  as.name(paste("illumina", annoName, "GENOMICLOCATION",sep=""))
          fn <- featureNames(from)
          fn <- fn[which(fn %in% mappedkeys(eval(mapEnv)))]
          
          locs <- mget(fn,eval(mapEnv),ifnotfound=NA)
          
          locs <- lapply(locs, function(x) gsub(" ", ",", x,fixed=T))
          
          asLocMatrix <- function(str){
            x<- do.call("rbind",sapply(strsplit(as.character(str), ",",fixed=T)[[1]], function(x) as.vector(strsplit(x, ":",fixed=T))))
          }
          
          locMat <- lapply(locs, asLocMatrix)
          
          rn <- rep(names(locs), unlist(lapply(locMat, nrow)))
          
          locMat <- do.call("rbind", locMat)
          
          grL <- GRangesList()
          for(i in 1:ncol(from)){
          rng <- GRanges(locMat[,1], IRanges(as.numeric(locMat[,2]), as.numeric(locMat[,3]),names=rn),strand=locMat[,4])
          #mcols(rng) <- df[match(names(rng), rownames(df)),]
          
          mcols(rng)$LogFC <- as.numeric(as.vector(LogFC(from[rn,i])))
          mcols(rng)$LogOdds <- as.numeric(as.vector(LogOdds(from[rn,i])))
          mcols(rng)$PValue <- as.numeric(as.vector(PValue(from[rn,i])))
        
          grL[[i]] <- sort(rng)
          }
          names(grL) <- sampleNames(from)
          grL
        }
        
      }
      
)


setMethod("dim", "limmaResults", function(x) {
  
  nFeatures = nrow(fData(x))
  nSamps = length(sampleNames(x))
  
  c("Features"=nFeatures, "Contrasts"=nSamps)
} )



setGeneric("LogFC", function(object) standardGeneric("LogFC"))

setMethod("LogFC", signature(object="limmaResults"), function(object) assayDataElement(object, "LogFC"))

setGeneric("LogFC<-", function(object, value) standardGeneric("LogFC<-"))


setGeneric("LogOdds", function(object) standardGeneric("LogOdds"))

setMethod("LogOdds", signature(object="limmaResults"), function(object) assayDataElement(object, "LogOdds"))

setGeneric("LogOdds<-", function(object, value) standardGeneric("LogOdds<-"))



setGeneric("PValue", function(object) standardGeneric("PValue"))

setMethod("PValue", signature(object="limmaResults"), function(object) assayDataElement(object, "PValue"))

setGeneric("PValue<-", function(object, value) standardGeneric("PValue<-"))




setGeneric("ArrayWeights", function(object) standardGeneric("ArrayWeights"))

setMethod("ArrayWeights", signature(object="limmaResults"), function(object) object@ArrayWeights)

setGeneric("ArrayWeights<-", function(object, value) standardGeneric("ArrayWeights<-"))
setReplaceMethod("ArrayWeights",
                 signature=signature(
                   object="limmaResults",
                   value="numeric"),
                 function(object, value) {
                   object@ArrayWeights <- value
                   object
                 })




setGeneric("DesignMatrix", function(object) standardGeneric("DesignMatrix"))

setMethod("DesignMatrix", signature(object = "limmaResults"), function(object) object@DesignMatrix)
setGeneric("DesignMatrix<-", function(object, value) standardGeneric("DesignMatrix<-"))

setReplaceMethod("DesignMatrix",
                 signature=signature(
                   object="limmaResults",
                   value="matrix"),
                 function(object, value) {
                   object@DesignMatrix <- value
                   object
                 })


setGeneric("ContrastMatrix", function(object) standardGeneric("ContrastMatrix"))

setMethod("ContrastMatrix", signature(object = "limmaResults"), function(object) object@ContrastMatrix)
setGeneric("ContrastMatrix<-", function(object, value) standardGeneric("ContrastMatrix<-"))

setGeneric("SampleGroup", function(object) standardGeneric("SampleGroup"))

setMethod("SampleGroup", signature(object = "limmaResults"), function(object) object@SampleGroup)



setReplaceMethod("ContrastMatrix",
                 signature=signature(
                   object="limmaResults",
                   value="matrix"),
                 function(object, value) {
                   object@ContrastMatrix <- value
                   object
                 })

setMethod("show", signature(object="limmaResults"), function(object) {
  cat("Results of limma analysis\n")
  cat("Design Matrix used...\n")
  print(head(DesignMatrix(object)))
  cat(".....\n\n")
  cat("Array Weights....\n")
  cat("Contrast Matrix used...\n")
  print(head(ContrastMatrix(object)))
  cat("Array Weights....\n")
  cat(selectSome(round(ArrayWeights(object),digits=3)))
  cat("\nTop Table\n")
  
  for(i in 1:ncol(object)){
    cat(paste("Top 10 probes for contrast", sampleNames(object)[i], "\n"))
    topN <- order(LogOdds(object)[,i],decreasing=T)[1:10]
    sub <- object[topN,]
    df <- data.frame(fData(sub), LogFC=LogFC(sub)[,i], LogOdds=LogOdds(sub)[,i], pvalue = PValue(sub)[,i])
    print(head(df,4))
    cat("\n\n")
    Direction <- rep(0, nrow(object))
    sig <- which(p.adjust(PValue(object)[,i]) <0.05)
    if(length(sig)>0){
      Direction[sig] <- sapply(LogFC(object)[sig,i], function(x) ifelse(x>0,1,-1))
    cat("Significant probes with adjusted p-value < 0.05\n")
    }
  
    print(table(Direction))                    
    
    cat("\n\n")
  }
  
})

setMethod("plot",
          signature(x = "limmaResults"),
          function (x) 
          {
            df <- NULL
            for(i in 1:ncol(x)){
              df[[i]] <-  data.frame(LogFC = LogFC(x)[,i], LogOdds = LogOdds(x)[,i],Contrast=sampleNames(x)[i])
              
            }
            df<-do.call("rbind",df)
            
            ggplot(df, aes(x = LogFC, y = LogOdds)) + geom_point(color="steelblue",alpha=0.3) + facet_wrap(~Contrast)
            
          })
          



