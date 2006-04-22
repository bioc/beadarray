\name{readBeadImages}
\alias{readBeadImages}
\alias{calculateBackground}
\alias{calculateForegroundIntensities}
\alias{sharpen}
\alias{averagePixels}
\alias{read.pgmfile}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read BeadArray information into R}
\description{
Uses image and csv files to load information about each array used in the experiment. 

}
\usage{
readBeadImages(targets, path=NULL, sharpen=TRUE, backgroundSize=17, columns=list(ProbeID="Code", x="x", y="y"),numrow=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targets}{targets vector which defines the images and csv files used to read each array}
  \item{path}{character string specifying the location of files to be read by the function}
  \item{sharpen}{if TRUE then the Illumina sharpening mask will be used prior to foreground averages being calculated}
  \item{backgroundSize}{numeric value which defines the size of the n x n box we will use to calculate local background values}
  \item{columns}{vector specifying the column headings to be read from the csv files}
  \item{numrow}{used to specify the number of rows (ie beads) each input file will contain. If specified in advance this can speed the 
reading of the data.}	

}
\details{
The targets objects can be read in previously by using a 'read.table' command. This defines the location of the csv and pgm files that are required. The column headings
Image1 and xyInfo must be found in the targets object. csv files are used to specify the location and identity of each bead on the array and must contain columns for 
the $x$ and $y$ postion of each bead as well as a ProbeID. The names of these columns can be specified by setting the 'columns' parameter

The foreground and background intensities of each bead are then calculated and stored separately. For foreground calculations the sharpening mask used by Illumina 
is used prior to an averaging over the 9 pixels in a 3 x 3 square closest to the bead centre (we can choose to calculate the values obtained without this mask). For 
the background intensity we look at the pixels in a 17 x 17 square and
calculate the mean of the 5 lowest values inside this square.  If a bead
is particularly close to the edge of the image then it isn't possible to
take a 17 x 17 square for sampling.  In this case the background value is
taken as the mean of all the previous background values and the user is
informed that this has been done.

The function creates an BeadLevelList containing foreground, background intensities for each bead on each array

Function requires the 'pixmap' library to be loaded

}
\value{
BeadLevelList object with list components R, Rb, G, Gb
}

\author{Mark Dunning, Mike Smith}

\examples{

#library(pixmap)

#targets = read.table("targets.txt", header=T)

#targets

#May take a while to run

#BLData.s = readBeadImages(targets)

#Create foreground intensities without using sharpening. Should take less time
#BLData.ns = readBeadImages(targets, sharpen=FALSE)

#Use a 9 x 9 are around each bead to calculate the background

#BLData.ns.9 = readBeadImages(targets, sharpen=FALSE, backgroundSize=9)

}
\keyword{IO}
