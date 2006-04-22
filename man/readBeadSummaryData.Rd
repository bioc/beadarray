\name{readBeadSummaryData}
\alias{readBeadSummaryData}
\alias{readCols}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read bead summary data information into R}
\description{
Function to read pre-processed bead summary data into R. Data can be the output from BeadStudio with either a separate file for each array in the experiment, or each file containing all arrays on a SAM.

}
\usage{

readBeadSummaryData(targets=NULL, header=T, sep=",",path=NULL,
                                   columns = list(ProbeID = "TargetID",
                                     AvgSig = "AVG_Signal", Nobeads = "Avg_NBEADS",
                                     Detection="Detection", BeadStDev="BEAD_STDEV"),
                                   other.columns = NULL, skip=7)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targets}{either a list of files created by
    readBeadSummaryTargets or  NULL to read all files in a directory.}
  \item{sep}{character used to seperate each column in the file. Usually either "," or "\t".}
  \item{columns}{list defining the names for each column to be read from the input files}
  \item{other.columns}{list defining any additional columns to be read from the input.} 
  \item{header}{if TRUE then function will read column headers from the input files.}
  \item{path}{optional argument to specify the location of the data to be read. By default, files are read from the current working directory.}
  \item{skip}{number of lines to ignore when reading the input files.}

}
\details{

  This function can be used to read the standard output from the
  Illumina BeadStudio application with default column names for Average
  Signal, Bead Stdev, Number of Beads and Detection and in the
  Column-oriented (Comma delimited) file format. ie the number of rows
  in each file is the number of bead types in the experiment and columns
  are repeated such that each row describes the information for each
  array in the experiment. The files are assumed to contain header
  information created by Illumina and this is skipped when the files are
  read. The skip parameter can be adjusted accordingly for files with a
  non-standard amount of header information.

  The names specified in the columns parameter are used to read
  the file. It is assumed that the average bead intensities for each
  array appear in columns with headings of the form \begin{verbatim}AVG_Signal_ARRAY1,
  AVG_Signal_ARRAY2....AVG_Signal_ARRAYN\end{verbatim} for the N arrays found in the
  file. All other columns headings are found in the same way.

  A separate BeadSummaryList object (described below) is created for
  each file listed in the targets object.
  
}
\value{
An BeadSummaryList object containing the components

R data frame containing the averaged foreground intensities for each
bead type for each array after outlier removal
BeadStDev data frame containing the standard deviation for each bead type for each array
NoBeads data frame containing the number of beads used to calculate the average for each bead type
Detection data frame containing the average detection score for each bead type for each array
other	data frame containing matrices corresponding to other.columns if given
}

\author{Mark Dunning and Mike Smith}

\seealso{
  \code{\link[beadarray:BeadSummaryList]{BeadSummaryList}}
  \code{\link{readBeadSummaryTargets}}
}
\examples{
#To read all bead summary files in current directory simply use
#
#BSData = readBeadSummaryData()
#
#if the files are contained in subdirectory New Folder
#BSData = readBeadSummaryData(path="New Folder")
#
#if all header information has been removed...
#BSData = readBeadSummaryData(skip=0)
#
#
#targets=readBeadSummaryTargets("targets.txt")
#
#BSData=readBeadSummaryData(targets)
#
#
#to also read the MAX_Signal column for each array 
#
#BSData = readBeadSummaryData(other.columns="MAX_Signal")



}
\keyword{IO}