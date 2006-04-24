\name{plotOnSAM}
\alias{plotOnSAM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Show variation between all 96 arrays on a SAM}
\description{
Function to show how quantities vary across all 96 arrays. eg Mean intensity of a certain control probe
}
\usage{
plotOnSAM(values, mx = max(values, na.rm = TRUE), scale = max(values,
na.rm = TRUE), min = 0, main = NULL, label = TRUE, missing_arrays = NULL, colour=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{vector containing 96 numeric values to plot}
  \item{mx}{maximum value to display on y axis of plot}
  \item{scale}{numeric value giving the amount by which to divide all numbers by (eg for log2 intensities this should be 16) to transform to range 0 - 1}
  \item{min}{numeric value giving the minimum value to display on y axis}
  \item{main}{character string giving a title for the plot}
  \item{label}{boolean defining if the arrays are labeled on the plot}
  \item{missing_arrays}{vector of numeric values specifying the index of any arrays that have been removed
    from the SAM.}
  \item{colour}{if TRUE the hexagons will be plotted in colour}
}
\details{
Two plots are produced side-by-side. The first is a plot of the set of values against the index 1-96 and secondly we plot 8 x 12 hexagonal arrays with array 
number 1 being the hexagon in the top-left corner and array 96 in the bottom-right. The colour of hexagon is directly related to the value in v for the particular 
array number. An array which has a higher value in v will be coloured brighter.

The values in v could be the values of a particular control probe across all arrays which can be obtained by the getMeanIntensities function. 
The findAllOutliers function could also be used repeatedly for all arrays and then plotted using this function.
}
\value{
Plot is produced on current graphical device.
}

\author{Mark Dunning}

\keyword{hplot}