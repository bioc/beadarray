\name{findAllOutliers}
\alias{findAllOutliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find outliers on array}
\description{
Function to find all beads which are outliers for their particular bead type on a given array.

}
\usage{
findAllOutliers(BLData, array, log=FALSE, n=3, ignoreList=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{an BeadLevelList object containing bead level data}
  \item{array}{numeric value specifying which array we want to find outliers on}
  \item{log}{if TRUE the log2 intensities will be used to calculate outliers}
  \item{n}{Specify a cut-off for outliers as n median absolute deviations from the mean. The default value is 3}
\item{ignoreList}{list of ProbeIDs to be omitted from the averaging procedure. These could be Illumina internal controls which are replicated many thousands of times on arrays}
}
\details{

We find the outliers for each bead type on the array in turn using the \code{\link{findBeadStatus}} function and keep a list of the outliers found. By default, outliers for a particular bead type are determined using a 3 MAD cut-off from the mean.

}

\value{
 numeric vector giving a row index from the BeadLevelList for all the beads that are outliers for their bead type.
}

\author{Mark Dunning and Mike Smith}

\seealso{\code{\link{findBeadStatus}}}

\examples{

data(BLData)

length(findAllOutliers(BLData, 1))

length(findAllOutliers(BLData, 1, log=TRUE))

length(findAllOutliers(BLData, 1, n=5))

}

\keyword{methods}