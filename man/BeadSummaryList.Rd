\name{BeadSummaryList-class}
\docType{class}
\alias{BeadSummaryList-class}
\title{Stores bead summary data - class}

\description{
A list based class for storing bead summary data.  \code{BeadSummaryList} objects can be created from \code{BeadLevelList} objects or by reading 
directly from pre-processed BeadStudio output by using the \code{\link{readBeadSummaryData}} function.

}

\section{Slots/List Components}{
Objects should contain the following list components:
\tabular{ll}{
\code{R}: \tab Numeric matrix containing foreground intensites for the
red channel.  Rows correspond to beadtypes and columns to arrays.\cr
\code{Rb}: \tab Numeric matrix containing background intensites for the
red channel. \cr
\code{BeadStDev}: \tab Standard deviations for each bead type per array.\cr
\code{NoBeads}: \tab Number of beads of each type per array.\cr
\code{ProbeID}: \tab ID for each probe on the arrays.\cr
\code{nooutliers}: \tab The number of beads flagged as outliers for each
bead type per array. This is only used for bead summary data created by beadarray and not for pre-processed data read from BeadStudio
}
}

\section{Methods}{
\code{BeadSummaryLists} can be subsetted and combined.  They also return dimensions
so functions such as\code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are also defined.  \code{BeadSummaryList}
uses an adapted version of the \code{\link[methods]{show}} method from the Limma class \code{\link[limma:LargeDataObject]{LargeDataObject}}.  This
means that the \code{BeadSummaryList} will print in a relatively compact way.

Other functions which operate on \code{BeadSummaryList} objects include:
\code{\link{plotXY.beads}}
\code{\link{plotMA.beads}}
\code{\link{getMeanIntensities}}
\code{\link{pairwiseMA}}
}

\author{Mark Dunning and Mike Smith}

\keyword{classes}
\keyword{data}