\name{squeezedVarOutlierMethod}
\alias{squeezedVarOutlierMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifier outliers on an array section
}
\description{

An outlier calling method that shrinks the observed variance for a bead-type towards the predicted variance based on all bead-types on the array-section.  

}
\usage{
squeezedVarOutlierMethod(inten, probeList, wts=1, n = 3, predictNlim=14)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inten}{
a list of intensities
}
  \item{probeList}{
the IDs corresponding to each intensity value
}
  \item{wts}{
Weights associated with beads, indicating those recommended for removal by, 
for example, \code{\link{BASH}}
}
  \item{n}{
number of SDs cutoff used
}
  \item{predictNlim}{
how many beads of a bead-type must be present for that bead-type to contribute
to prediction of variances?
}
}

\details{
This function is called within the \code{\link{summarize}} routine of beadarray to exclude outlying beads from an array-section prior to summary. The intensities are not assumed to be on any particular scale and can result from any user-defined transformation function, however a log-transformation is recommended. 

Bead-types that have \code{predictNlim} numbers are used to locally regress bead-type precision against bead-type mean, as well as the squared residual error of bead-type precision against bead-type mean. These are then used as prior values for the distribution of precision to feed into a standard Bayesian calcuation to obtain an estimate of the posterior variance.

Beads with weight zero do not contribute to the outlier calling.
}


\value{
the positions in the original vector that were determined to be outliers

}
\author{
Andy Lynch
}

\seealso{\code{\link{illuminaOutlierMethod}}}

\examples{

if(require(beadarrayExampleData)){

data(exampleBLData)

oList = squeezedVarOutlierMethod(
        logGreenChannelTransform(exampleBLData, 1), 
        getBeadData(exampleBLData, array=1, what="ProbeID")
    )
 
}

}

