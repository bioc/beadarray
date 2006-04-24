\name{readBeadSummaryTargets}
\alias{readBeadSummaryTargets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create targets object}

\description{

Read targets file for a BeadArray experiment into a dataframe.

}
\usage{
readBeadSummaryTargets(file = "beadSummaryTargets.txt", path = NULL,  header=T,sep="" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{character string giving the name of the targets file.}
  \item{path}{character string giving the directory containing the file.
          Can be omitted if the file is in the current working
          directory.}
  \item{header}{if TRUE then header names will be taken from the text file}
  \item{sep}{field separator character}
  
}
\details{
A targets file defines the files containing bead summary information for a BeadArray experiment. Each file may describe either a single array, or a 96-well SAM.
The function looks for a file called beadTargets.txt in the current working directory. It should contain a column heading FileName with each of the files to be read listed below.

}
\value{
A dataframe containing the columns found in the text file
}

\author{Mark Dunning, Mike Smith}

\examples{

#beadTargets = readBeadSummaryTargets()

#can also specify path

#beadTargets = readBeadTargets(file="otherfile.txt", path="otherdirectory")

}
\keyword{IO}
