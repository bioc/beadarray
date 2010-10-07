\name{calculateOutlierStats}
\alias{calculateOutlierStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outlier distribution stats
}
\description{
Function that determines the outlier beads on an array and how they are distributed among the segments

}
\usage{
calculateOutlierStats(BLData, array = array, transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, useLocs = TRUE, nSegments = 9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
a \code{\link{beadLevelData}} object
}
  \item{array}{
the number of the array of interest
}
  \item{transFun}{
how the section data is to be transformed prior to calculating outliers
}
  \item{outlierFun}{
a function for calculating outliers
}
  \item{useLocs}{
use locs and sdf information (if available) to determine section layout
}
  \item{nSegments}{
manually set how many segments the section is divided into
}
}
\details{

A section of an expression BeadChip (e.g. the Humanv3 or HumanHT-12) is made up of 9 physcially-separate segments. A useful QA check is to see how the outliers are distributed among these segments. Outliers are beads that have outlying intensities according to some rule that the user can specify. The default (as used by Illumina) is to exclude beads that are more than 3 median absolute deviations from the median.

Note that sections from Sentrix Array Matrix do not have segments, so the results may not be informative

}


\value{
vector with the percentage of beads found in each segment that were determined to be outliers

}


\author{
Mark Dunning
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

data(BLData)

##Artificical example, there are no segments on this type of BeadArray

calculateOutlierStats(BLData, array=1, nSegments=10, useLocs=FALSE)

calculateOutlierStats(BLData, array=2, nSegments=10, useLocs=FALSE)

calculateOutlierStats(BLData, array=3, nSegments=10, useLocs=FALSE)


}
