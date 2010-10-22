\name{beadLevelData-class}
\docType{class}
\alias{beadLevelData-class}
\alias{beadLevelData}
\alias{show,beadLevelData-method}
\alias{dim,beadLevelData-method}
\alias{numBeads,beadLevelData-method}
\alias{[[,beadLevelData,ANY,missing-method}
\alias{boxplot,beadLevelData-method}
\alias{sectionNames,beadLevelData-method}
\alias{combine,beadLevelData,beadLevelData-method}

\title{Class "beadLevelData" }

\description{
    A class for storing red and green channel foreground and background intensities from an Illumina experiment.
}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("beadLevelData")}, but are usually created by \code{\link{readIllumina}}.
}

\section{Slots/List Components}{
  Objects of this class contain the following slots
  \tabular{ll}{
  \code{beadData}:\tab A list of arrays, indexed by array name.  Each item in this list is itself a list, containing enviroments holding the data for that array.\cr
  \code{sectionData}:\tab a list containing information. Each item in the list is a data frame containing one row for each section\cr
  \code{experimentData}:\tab a list containing the annotation of the platform, link to the sdf file and type of data (slide or Sentrix Array Matrix) \cr
  \code{history}: \tab Character vector storing the operations performed on this object.
  }
}

\section{Methods}{
  \describe{
    \item{show(beadLevelData)}{Printing method for BeadLevelList}
    \item{\code{dim(object)}}{The dimension of the \code{BeadLevelList} object}
    \item{\code{sectionNames(object,arrays=NULL)}}{Returns the strip/array names from a}
    \item{\code{numBeads(object,arrays=NULL)}}{Returns the number of beads on selected \code{arrays}}
  }
}

\section{Accessing data from the class}{
  \describe{
    \item{\code{\link{getBeadData}}}{retrieve data}
    \item{\code{\link{insertBeadData}}}{Input or modify existing data}
  }
}

\examples{

data(BLData)

sectionNames(BLData)

head(BLData[[1]])

getBeadData(BLData, array=1, what="Grn")[1:10]

}
\author{Mark Dunning, Mike Smith}
\seealso{\code{\link{readIllumina}}}
\keyword{classes}
