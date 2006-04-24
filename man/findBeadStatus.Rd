\name{findBeadStatus}
\alias{findBeadStatus}
\alias{getProbeIndices}
\alias{getProbeIndicesC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find Outliers }
\description{
 Function finds all beads which are outliers for a given bead type
}
\usage{
findBeadStatus(BLData, probes, array, log=FALSE, n=3, outputValid = FALSE, intProbeID = NULL, ignoreList=NULL, probeIndex=NULL, startSearch=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object containing bead level data}
  \item{probes}{numeric vector for the ProbeIDs of the bead type we want to find outliers for}
  \item{array}{numeric value for the array we want to find outliers on}
  \item{log}{if TRUE the log2 intensities of each bead will be used to find outliers. Otherwise un-logged data is used}
  \item{n}{numeric value defining a cut-off of 'n' median absolute deviations from the mean to determine outliers. The default value is 3}
  \item{outputValid}{if TRUE the IDs of beads which are not outliers will be output}
  \item{intProbeID}{BLData\$ProbeID coerced to vector of integers.  Never change this, for internal use only}
\item{ignoreList}{list of ProbeIDs to be omitted from the averaging procedure. These could be Illumina internal controls which are replicated many thousands of times on arrays}
\item{probeIndex}{parameter for internal use only}
\item{startSearch}{integer specify where to start searching for a particular ProbeID}

}
\details{

The intensities of each bead with ProbeID 'probe' on the specified array are found and if the 'log' parameter is set to true we do a log2 transformation
on the set of values.

The mean and median absolute deviation of all the bead intensities are then calculated. Outliers are beads which have intensity more than 
'n' median absolute deviations from the mean.

The default setting for the function is the method used by Illumina of un-logged intensities and 'n' = 3.

Any beads which have intensity NA are also counted as outliers

The function returns only the outliers for a bead type unless the
outputValid parameter is specified.

}

\value{
  
  List of of row indices relating to the BeadLevelList dividing the beads of this bead type into two categories.

  \item{valid}{valid beads}
  \item{outliers}{beads which are calculated as outliers}
}

\author{Mark Dunning}

\seealso{ \code{\link{findAllOutliers}} }

\examples{

data(BLData)

findBeadStatus(BLData, 2, 1)

findBeadStatus(BLData, 2, 1, log=TRUE)


##Methods may not always agree

findBeadStatus(BLData, 6140, 1)

findBeadStatus(BLData, 6140, 1, log=TRUE)

##Increasing value of n selects less outliers

findBeadStatus(BLData, 6140, 1, n=4)

findBeadStatus(BLData, 6140, 1, log=TRUE, n=4)

}

\keyword{methods}