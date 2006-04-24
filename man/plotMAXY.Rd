\name{plotMAXY}
\alias{plotMAXY}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatter plots and MA-plots for all specified arrays }
\description{

Produces scatter and MA-plots for all pairwise comparisons that can be made from a set of arrays.
}
\usage{
plotMAXY(BSData, vec, log = TRUE, labels=vec, label=FALSE,ma.ylim = 2, sampleSize = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BSData}{a BeadSummaryList containing bead summary data}
  \item{vec}{numeric vector giving the indices of the arrays in BSData we want to plot}
  \item{log}{if TRUE then all values will be log-transformed before plotting}
  \item{labels}{vector specifying names for each array to be plotted}
  \item{label}{if TRUE any controls which have been defined in a genes
    object will be plotted}
  \item{ma.ylim}{specifies a height for the y axis (M values) of the MA-plots. Default is 2}
  \item{sampleSize}{Defines the number of genes to plot.  Default is NULL, which plots everything}
  \item{\dots}{other graphical parameters to be passed}
}

\details{
These are graphical tools for showing the differences that exist between two arrays and they can be used to highlight problems such as bias between arrays as well as identifying genes which are differentially expressed. We constructed these plots as follows. For each bead type, we take the averaged intensity calculated from two different arrays. On an XY plot these values are simply plotted against each other and on an MA plot we plot the mean of the averaged intensities against the corresponding ratio of the average intensities.

In the lower-left section of the plot we see XY plots for all pairwise comparisons between the arrays and in the upper right we have pairwise MA plots. Going down the first column we observe XY plots of array 1 against array 2 and array 1 against array 3 etc. Similarly, in the upper-right corner we can observe pairwise MA plots.

}
\author{Mark Dunning}

\examples{

#data(BSData)

#vec = c (1,2,3)

#plotMAXY(BSData, vec)




}
\keyword{hplot}
