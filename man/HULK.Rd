\name{HULK}
\alias{HULK}
\title{HULK - Bead Array Normalization by NEighbourhood Residuals}
\description{
Normalizes an probe intensities by calculating a weighted average
residual based on the residuals of the surrounding probes.
}
\usage{
HULK(BLData, array = 1, neighbours = NULL, invasions = 20, useLocs = TRUE, weightName = "wts", transFun = logGreenChannelTransform)
}

\arguments{
  \item{BLData}{An object of class \code{\link{beadLevelData}}}
  \item{array}{integer specifying which strip/array to process}
  \item{neighbours}{A Neighbours matrix. Optional - if left NULL, it will be computed.}
  \item{invasions}{Integer - Number of invasions used when identifying
    neighbouring beads.}
  \item{useLocs}{If information from an associated locs is to be used}
  \item{weightName}{Column name where bead weights are to be taken from}
  \item{transFun}{Transformation function}
}

\details{
  HULK is a method of intensity normalization based upon the
  BASH framework.  Firstly For each bead a local neighbourhood of beads is
  determined, using the same process as the other BASH functions.

  For each bead a weighted average residual is calculated.  The average
  residual is calculated as the sum of the residuals for each bead in
  the neighbourhood, divided by 1 plus the number of invasions it took to reach
  that bead.  This calculation is made by a call to
  \code{HULKResids}.

  The average residuals are then subtracted from each bead and the
  resulting corrected intensities object are returned. These corrected intensities can be saved in the original beadLevelData object using \code{insertBeadData}
}

\value{
List where each entry is either a vector of corrected intensities or NULL, depending upon which arrays were specified.
}

\author{Mike Smith}

\seealso{\code{\link{BASH}}}

\examples{
data(BLData)
o <- HULK(BLData, 1)


}
\references{
Lynch AG, Smith ML, Dunning MJ, Cairns JM, Barbosa-Morais NL, Tavare S. 
beadarray, BASH and HULK - tools to increase the value of Illumina BeadArray experiments. 
In A. Gusnato, K.V. Mardia, & C.J. Fallaize (eds), Statistical Tools for Challenges in Bioinformatics. 2009 pp. 33-37. Leeds, Leeds University Press.

}
\keyword{misc}
