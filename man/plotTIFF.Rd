\name{plotTIFF}
\alias{plotTIFF}
\title{Produce plots of the Illumina tiff images}
\description{
  Produces a plot of an Illumina tiff image, which can be useful for
  observing spatial artifacts on an array and checking the alignment of
  spot centres features in the image.
}
\usage{
    plotTIFF(tiff, xrange = c(0, ncol(tiff)-1), yrange = c(0, nrow(tiff)-1), 
    high = "cyan", low = "black", mid = NULL, ncolours = 100, 
    log = TRUE, values = FALSE, textCol = "black", ...)
}

\arguments{
  \item{tiff}{Intended to the the result of \code{\link{readTIFF}}, but
    in reality can be any matrix}
  \item{xrange}{Range of X coordinates to plot.}
  \item{yrange}{Range of Y coordinates to plot.}
  \item{high}{Colour to plot the brightest pixels in the image.}
  \item{low}{Colour to plot the dimmest pixels.}
  \item{mid}{If specified the colour gradient will go from low to mid to
  high.  If not specified then the gradient simply goes from low to
  high.}
  \item{ncolours}{Specify how many steps there should be in the gradient between the high and low colours}
  \item{log}{If TRUE the pixel values are logged before the colour
    gradient is created.}
  \item{values}{When set to TRUE each pixel in the image has it's value
    displayed over it.  This should only be used when displaying a very
    small number of pixels as the text very quickly covers the entire image.}
  \item{textCol}{If \code{values} is \code{TRUE} this argument specifies the colour of the text.}
  %\item{accountForZero}{Sometimes the Illumina tiff images have a very
  %  small number of pixels with an intensity value of zero.  \code{plotTIFF}
  %  scales the colours from the smallest value to the largest.  The
  %  presence of these zeros affects this scaling, resulting in a much
  %  smaller range of colours covering the majority of the pixel values.
  %  This in turn gives an image where it is hard to discern the spots from
  %  the background.  Setting this argument to TRUE adjusts all pixels with a value of zero
  %  to the second smallest value, which should help restore the range of colours.  NEEDS WORK!
  %}
  \item{...}{Other graphical parameters specified in \code{par}}
}
\details{


  This can be very slow, especially when the Cairo graphics library is
  being used.  When using the Cairo library, if one is plotting a large tiff
  with 10s of millions of pixels, the plotting time increases from
  around 20 seconds to 5 minutes on an Intel Xeon E5420.

  If running on a Linux system it is recommended to use:

  \code{x11(type = "Xlib")}

  before running \code{plotTIFF()}, on order to force the quicker
  plotting mechanism.

  Of course it is debatable whether it is useful to plot all of those
  pixels, given that there are far more than can be displayed on a
  normal screen, and future revisions of the code may address this.
}
\value{
A plot is produced on the current graphical device.
}

\author{Mike Smith}

\keyword{hplot}
