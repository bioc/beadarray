\name{readNonXYData}
\alias{readNonXYData}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read BeadArray information into R}
\description{
Read csv files which contain bead-level data without x and y positions for each bead. An BeadLevelList can be created with limited functionality available 
or used to create an BeadSummaryList (bead summary data).


}
\usage{
readNonXYData(files, max_length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{List of filenames to be read}
  \item{max_length}{numeric value specifying the length of the longest file to be read}	

}
\details{
The function creates an BeadLevelList containing foreground and ProbeID for each bead on each array


}
\value{
BLData object with list components R and ProbeID
}

\author{Mark Dunning}

\keyword{IO}
