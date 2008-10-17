\name{outlierPlot}
\alias{outlierPlot}

\title{Locations of outliers on an Illumina chip}
\description{
Diagnostic function that reports how many outliers are found on a specific array on a chip. We take advantage of the segmental structure of the array and break-down the number of outliers into 9 sections.
}
\usage{
outlierPlot(BLData, array = array, log = FALSE, plot = FALSE)
}

\arguments{
  \item{BLData}{A BeadLevelList object containing the bead-level data for an Illumina experiment}
  \item{array}{The number of the array of interest}
  \item{log}{if TRUE calculate outliers on the log2 scale.  If FALSE calculate outliers on the original scale}
  \item{plot}{if TRUE a diagnostic plot will be produced, otherwise only the numbers of outliers will be returned.}
}
\details{
The number of outliers are computed for the whole array using the Illumina default method that specifies a cut-off of 3 MADs from the median on either the log2 or original scale. These outliers are then split into 9 different sections on the array (the separation between these sections can usually be seen in the plots).
}

\value{

}
\author{Mark Dunning}
