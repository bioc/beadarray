\name{imageplot}
\alias{imageplot}
\title{imageplot for beadLevelData object}
\description{
Generates an image plot for data from a \code{beadLevelData} object.
}


\usage{
imageplot(BLData, array = 1, transFun = logGreenChannelTransform, squareSize = NULL, useLocs = TRUE, horizontal = TRUE, low = NULL, high = NULL, ncolors = 100, zlim = NULL, legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{
\code{beadLevelData}
}
  \item{array}{
integer specifying what section to plot
}
  \item{transFun}{
Function that defines how values from the BLData object are to be transformed prior to plotting. 
}
  \item{squareSize}{
Numeric specifying how many pixels in the original image make up each square in the imageplot. If NULL, the function will guess a suitable value from the data. 
}
  \item{useLocs}{
If \code{TRUE} the function will read the locs file associated with the section in order to include the physical properties of the section in the plot
}
  \item{horizontal}{
If \code{TRUE} the image will be plotted so that the longest edge of the section is on the x axis.

}
  \item{low}{
colour to use for lowest intensity
}
  \item{high}{
colour to use for highest intensity
}
  \item{ncolors}{
The number of colour graduations between high and low
}
  \item{zlim}{
numerical vector of length 2 giving the extreme values of 'z'
    to associate with colours 'low' and 'high'.
}
  \item{legend}{
logical, if \code{TRUE}, \code{zlim} and range of data
    is added to plot.
}
  \item{\dots}{
other arguments to plot
}

}
\details{
Produces a standard imageplot for the specified section. The default, transformation \code{\link{logGreenChannelTransform}}, takes the log2 of the green channel. For two channel data, the red channel or log ratio can be plotted by \code{\link{logRedChannelTransform}} or \code{\link{logRatioTransform}} functions can be used. The user can also specify their own functions. 

The default plotting orientation is such that the longest edge of the section is along the x axis. If \code{horizontal = FALSE}, the longest edge will by on the y axis and should match how the corresponding TIFF image from the BeadScan directory is orientated. 

If \code{locs = TRUE} and locs file were made available to \code{readIllumina}, the segments that the section is comprised of will be visible (For expression BeadChips, each section is made of nine physically separate segments). The squareSize parameter will also be set appropriately.


As a result of both having identical function names this function can conflict with the  \code{\link[limma]{imageplot}} method in 'limma'.  If both packages are loaded, the function from whichever package was loaded last takes precedence.  If the 'beadarray' imageplot() function is masking that from 'limma', one can directly call the 'limma' method using the command "limma::imageplot()".  Alternatively, one can detach the 'beadarray' package using "detach(package:beadarray)".  Similar techniques can be used if 'limma' is masking the 'beadarray' method.
}
\value{
A plot is produced on the current graphical device.
}
\author{
Mike Smith and Mark Dunning
}


\examples{
data(BLData)

##Be default the first array is plotted
imageplot(BLData, horizontal = FALSE)

##Can use the squaresize parameter

imageplot(BLData, horizontal=FALSE, squareSize=10)

imageplot(BLData, array=3, horizontal=FALSE, squareSize=10)

##User can specifiy what colours to represent low and high intensity

imageplot(BLData, array=3, horizontal=FALSE, squareSize=10, low="lightgreen", high="darkgreen")


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplots}

