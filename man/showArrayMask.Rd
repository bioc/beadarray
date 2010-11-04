\name{showArrayMask}
\alias{showArrayMask}

\title{Show Array Mask}
\description{
Function to display beads masked by BASH. The masked beads are assumed to have a weight of 0 in the specified weights column.
}

\usage{

showArrayMask(BLData,array = 0, override = FALSE, wtsName = "wts", transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod, horizontal=TRUE)

}

\arguments{
  \item{BLData}{A BeadLevelList object.}
  \item{override}{Logical. Plotting a large mask can cause slowdown problems. By default, if more than 200 000 beads are masked, the current mask will not be plotted. You can force the mask to be plotted by setting this argument to \code{TRUE}, however beware as this may cause slower systems to freeze.}
  \item{wtsName}{name under which the bead weights are stored}
  \item{array}{numeric index of the array to plot}
	
  \item{transFun}{function to transform intensities prior to calculating outliers}
  \item{outlierFun}{function to remove outliers}
  \item{horizontal}{if TRUE the resulting image is plotting with the longest edge along the x axis}	

}

\details{

\code{showArrayMask} plots the beads on an array that have been assigned a weight of 0 by BASH in red, and beads determined to be outliers in black.


}

\value{
  None returned
}


\author{Jonathan Cairns and Mark Dunning}

\examples{
data(BLData)

bsh = BASH(BLData, array=1)

BLData = setWeights(BLData, wts=bsh$wts, array=1)

showArrayMask(BLData,1)


}

\keyword{misc}
