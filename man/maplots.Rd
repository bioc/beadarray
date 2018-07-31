\name{plotMA-methods}
\docType{methods}

\alias{plotMA,ExpressionSetIllumina-method}


\title{
Function to construct the classic MA plots from a dataset
}
\description{
Note that a ggplot2 object is returned, which can be modified programatically.

The default operation is to construct a reference array be averaging all probes across the whole datasets, and then making an MA plot by comparing each array to this reference. A log2 transformation is performed, but this can be turned-off by setting \code{do.log=FALSE}. 

Alternatively, a SampleGroup may be specified. This is a character string that should match a column in the phenoData for the object. Then, the resulting MA plot will contain all pairwise comparisons between samples belonging to the same group. A list of ggplot objects is returned, each item in the list being a different group.
}

%- maybe also 'usage' for other objects documented here.

\details{
Densities are computed using the 'hexbin' package
}
\value{
A list of ggplot2 objects is returned.

}

\author{
Mark Dunning
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

if(require(beadarrayExampleData)){
data(exampleSummaryData)
rawdata <- channel(exampleSummaryData, "G")
 mas <- plotMA(rawdata,do.log=FALSE)
mas
##Smoothed lines can be added
mas + geom_smooth(col="red")
##Added lines on the y axis
mas + geom_hline(yintercept=c(-1.5,1.5),col="red")
##Changing the color scale
mas + scale_fill_gradient2(low="yellow",mid="orange",high="red")

mas <- plotMA(rawdata,do.log=FALSE,SampleGroup="SampleFac")

\dontrun{
mas[[1]]
mas[[2]]
}

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
