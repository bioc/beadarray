\name{expressionQCPipeline}
\alias{expressionQCPipeline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Flexible bead-level QC pipeline
}
\description{
Function to produce various QC plots and HTML summary pages for bead-level data.
}
\usage{
expressionQCPipeline(BLData, transFun = logGreenChannelTransform, qcDir = "QC", plotType = ".jpeg", horizontal = TRUE,controlProfile=NULL,beadLevelQC = NULL, bashMetrics = NULL, outlierMetrics = NULL, detectionMetrics=NULL,overWrite=FALSE,nSegments=9,outlierFun=illuminaOutlierMethod,tagsToDetect = list(housekeeping = "housekeeping", Biotin = "phage_lambda_genome", Hybridisation = "phage_lambda_genome:high"),zlim=c(5,7),positiveControlTags = c("housekeeping", "phage_lambda_genome"), hybridisationTags =  c("phage_lambda_genome:low", "phage_lambda_genome:med","phage_lambda_genome:high"), negativeTag= "permuted_negative", boxplotFun = logGreenChannelTransform, imageplotFun = logGreenChannelTransform)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{beadLevelData} object
}
  \item{transFun}{
what transformation function to apply
}
  \item{qcDir}{
a directory to write output to 
}
  \item{plotType}{
desired file extension for plots (jpeg or png)
}
  \item{horizontal}{
if TRUE imageplots and outlier plots are produced with longest edge on x axis
}
  \item{controlProfile}{
a data frame defining all control types 
}
  \item{beadLevelQC}{
(optinal) a character string giving the name under which bead-level qc is stored in sectionData
}
  \item{bashMetrics}{
(optinal) a character string giving the name under which BASH qc is stored in sectionData
}
  \item{outlierMetrics}{
(optinal) a character string giving the name under which outlier qc is stored in sectionData
}
  \item{detectionMetrics}{
(optinal) a character string giving the name under which detection qc is stored in sectionData
}
  \item{overWrite}{
if FALSE any plots that exist in the directory will not be recreated
}
  \item{nSegments}{
how many segments each section is divided into
}
  \item{outlierFun}{
a function to removed outlier
}
  \item{tagsToDetect}{
which control types to used in the detection metrics
}
  \item{zlim}{
the range of the imageplots
}
  \item{boxplotFun}{
what transformation function to be used in boxplots
}
  \item{imageplotFun}{
what transformation function to be used for imageplots
}
  \item{positiveControlTags}{
character strings defining which positive controls to plot
}
  \item{hybridisationTags}{
additional control types to be plotted 
}

  \item{negativeTag}{
character string to identify which control type in the control profile corresponds to negative controls
}
  \item{\dots}{
other plot arguments
}
}
\details{

This function is a convient way of automatically generating QC plots for each section within a \code{beadLevelData} object. The following plots are produced for each section. i) scatter plots of all bead observation of the positive controls. See \code{\link{poscontPlot}}. ii) Further scatter plots of other controls of interest using \code{\link{poscontPlot}}. iii) imageplot (\code{\link{imageplot}}) of section data after applying transformation function iv) plot of outlier locations using specified outlier function. A HTML page displaying all the plots is produced. 

After plots have been produced for each section, \code{\link{makeQCTable}} is run to make a table of mean and standard deviations for the defined control types, followed by the results of \code{\link{calculateOutlierStats}} and \code{\link{controlProbeDetection}} for each section and written to a HTML page in the requested directory.

The function should be able to run automatically for expression data that has its annotation stored using \code{\link{setAnnotation}} or using \code{\link{readIllumina}}. Otherwise the \code{controlProfile} data frame can be used to define the control types on the array and their associated ArrayAddressIDs. Similarly, the function assumes single-channel data but a transformation function can be passed.

}

\value{

}
\references{
}
\author{
Mark Dunning
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{poscontPlot}}
\code{\link{imageplot}}
\code{\link{outlierplot}}
\code{\link{controlProbeDetection}}

}
\examples{

data(BLData)
data(controlProfile)

#expressionQCPipeline(BLData, controlProfile=controlProfile, positiveControlTags=c("housekeeping", "biotin"), hybridisationTags = c("cy3_hyb", "high_stringency_hyb"), zlim=c(10,12), horizontal=F, negativeTag = "negative", tagsToDetect=list(Housekeeping = "housekeeping", Biotin = "biotin", hybridisation = "cy3_hyb"))


}
