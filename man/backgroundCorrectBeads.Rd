\name{backgroundCorrectBeads}
\alias{backgroundCorrectBeads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Background Correct a BeadLevelList object }
\description{
Adapted from 'limma' backgroundCorrect function to correct R intensities of a BLData object using the Rb background values.

}

\usage{
backgroundCorrectBeads(BLData, method = "subtract", offset = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{BeadLevelList object containing bead level data}
  \item{method}{character string specifying correction method.  Possible
          values are '"none"', '"subtract"', '"half"', '"minimum"')}
          
  \item{offset}{numeric value to add to intensities}
}
\details{
If 'method="none"' then the corrected intensities are equal to the
     foreground intensities, i.e., the background intensities are
     treated as zero. If 'method="subtract"' then this function simply
     subtracts the background intensities from the foreground
     intensities which is the usual background correction method. 

     The remaining methods are all designed to produce positive
     corrected intensities. If 'method="half"' then any intensity which
     is less than 0.5 after background subtraction is reset to be equal
     to 0.5. If 'method="minimum"' then any intensity which is zero or
     negative after background subtraction is set equal to half the
     minimum of the positive corrected intensities for that array. 
     

     The 'offset' can be used to add a constant to the intensities
     before log-transforming, so that the log-ratios are shrunk towards
     zero at the lower intensities. This may eliminate or reverse the
     usual 'fanning' of log-ratios at low intensities associated with
     local background subtraction.

}
\value{
An 'BeadLevelList' object in which components 'R' and 'G' are background
     corrected 
}


\author{Mark Dunning using limma function}

\examples{

data(BLData)

#default is to simply subtract Rb from R

BLData.bc = backgroundCorrectBeads(BLData)

min(BLData.bc$R)

#Use 'minimum' method to stop negative values

BLData.min = backgroundCorrectBeads(BLData, method="minimum")

min(BLData.min$R)

}

\keyword{methods}