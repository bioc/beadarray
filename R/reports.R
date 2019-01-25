#library(biovizBase)
#require(beadarrayExampleData)
#data(exampleSummaryData)

#exampleSummaryData.log2 <- channel(exampleSummaryData, "G")
#exampleSummaryData.norm <- normaliseIllumina(exampleSummaryData.log2)
#library(illuminaHumanv3.db)
#ids <- as.character(featureNames(exampleSummaryData.norm))
#qual <- unlist(mget(ids, illuminaHumanv3PROBEQUALITY, ifnotfound=NA))
#table(qual)

#rem <- qual == "No match" | qual == "Bad" | is.na(qual)


#exampleSummaryData.filt <- exampleSummaryData.norm[!rem,]

#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#tx <- TxDb.Hsapiens.UCSC.hg19.knownGene

#limmaRes <- limmaDE(exampleSummaryData.filt, SampleGroup="SampleFac")

#genes <- c("ERBB2", "HBE1", "HBG1", "HBZ", "ALB")
#data(genesymbol)
#geneList <- genesymbol[genes]

makeReport <- function(geneList, summaryData, limmaRes, tx,genome="hg19"){
  requireNamespace(Nozzle.R1)
  requireNamespace(ggbio)

  
  lgr <- as(limmaRes, "GRanges")
  prb.ranges <- lgr[[1]]
  r <- newCustomReport("My Report")
  
  ReportSum <- newSection("ReportSummary")
  
  regionz <- newTable(as.data.frame(geneList))
  
  p <- newParagraph("The following regions were used in the report..")
  ReportSum <- addTo(ReportSum, p,regionz)
  r <- addTo(r, ReportSum)
  
  for(i in 1:length(geneList)){
    
    s <- newSection(paste0("Report for ", names(geneList)[i]))
    
    gr <- geneList[i]
    
    pbs <- suppressWarnings(names(prb.ranges[prb.ranges %over% gr]))
    
    tab <- newTable(fData(summaryData)[pbs,])
    
    s <- addTo(s, tab)
    
    ss1 <- newSection("Probe Locations")
    
    

    
    p1 <- autoplot(tx, which = gr)
    p2 <- suppressWarnings(autoplot(prb.ranges[prb.ranges %over% gr]))
    outfile <- paste0("reports/",names(geneList)[i], "-locations.png")
    if(!is.null(genome)) {
      p.ideo <- plotIdeogram(genome="hg19",subchr = as.character(seqnames(gr)))
        
      t <- tracks(p.ideo,p1,p2)
    } else t <- tracks(p1,p2)
      
      ggsave(t, filename=outfile)
    
    
  
    locs <- newFigure(basename(outfile),
                   type = IMAGE.TYPE.RASTER, exportId = NULL,
                   protection = PROTECTION.PUBLIC)
    ss1 <- addTo(ss1, locs)
    
    ss2 <- newSection( "Boxplot" )
    
    df <- data.frame(melt(exprs(summaryData)[pbs,]))
    
    if(length(pbs)==1){
      df <- data.frame(df, SampleGroup)
      gg <- ggplot(df, aes(x= SampleGroup, y = value,fill=SampleGroup)) + geom_boxplot()
    } else{
      df <- data.frame(df, SampleGroup = SampleGroup[match(df$Var2, sampleNames(summaryData))])
      gg <- ggplot(df, aes(x= SampleGroup, y = value,fill=SampleGroup)) + geom_boxplot() + facet_wrap(~Var1)
    } 
    
    outfile <- paste0("reports/",names(geneList)[i], "-boxplot.png")
    ggsave(gg, filename=outfile)
    
    k <- newFigure(basename(outfile),
                   type = IMAGE.TYPE.RASTER, exportId = NULL,
                   protection = PROTECTION.PUBLIC)
    
    
    
    ss2 <- addTo(ss2, k)
    
    
    ss3 <- newSection("Differential Expression Stats")
    
    df <- NULL
    for(i in 1:ncol(limmaRes)){
      
    df[[i]] <- data.frame(Probe = pbs, LogFC = LogFC(limmaRes)[pbs,i], LogOdds = LogOdds(limmaRes)[pbs,i], PValue = PValue(limmaRes)[pbs,i], Contrast = sampleNames(limmaRes)[i])
    
    }
    df <- do.call("rbind",df)
    detab <- newTable(df)
    
    ss3 <- addTo(ss3, detab)
    
    s <- addTo(s, ss1)
    s <- addTo(s, ss2)
    s<- addTo(s,ss3)
    

    r <- addTo(r, s)
  }
  
  
  writeReport(r, filename="reports/my_report")
  
}


