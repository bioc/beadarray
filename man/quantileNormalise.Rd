\name{quantileNormalise}
\alias{quantileNormalise}
\alias{normalizeBetweenArrays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantile Normalisation for BeadLevelList objects }
\description{
Normalizes expression intensities so that the intensities or
     log-ratios have similar distributions across a series of arrays.

Function is based on the 'normalizeBetweenArrays' function in 'limma'
}
\usage{
quantileNormalise(BLData, arrays=1:length(BLData$R[1,]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BLData}{can be either an BeadLevelList or BeadSummaryList object}
  \item{arrays}{choose a subset of arrays to normalise} 
}
\details{
Quantile normalization was proposed by Bolstad et al (2003) for
     Affymetrix-style single-channel arrays and by Yang and Thorne
     (2003) for two-color cDNA arrays. 'method="quantile"' ensures that
     the intensities have the same empirical distribution across arrays
     and across channels. 
}

\value{
Produces an BeadLevelList or BeadSummaryList (depending on input type) containing normalised intensity values.}
\author{Mark Dunning}

\seealso{
\code{\link{qsplineNormalise}}
\code{\link{medianNormalise}}
}

\examples{

data(BSData)

BSData.q = quantileNormalise(BSData)

#All distributions are fixed to be the same

boxplot(BSData.q$R~col(BSData.q$R))

}
\keyword{methods}

