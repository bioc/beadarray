\name{summarize}
\alias{summarize}
\title{
    Create a summarized object
}
\description{
Function to summarize the data in a \code{beadLevelData} object into a form more ameanable for downstream analysis (with the same number of observations for each bead type).
}


\usage{
    summarize(BLData, channelList = list(greenChannel), 
    probeIDs=NULL, useSampleFac = FALSE, sampleFac= NULL,
    weightNames = "wts", removeUnMappedProbes = TRUE)
}

\arguments{
    \item{BLData}{
        An object of class beadLevelData
    }
    \item{channelList}{
	List of objects of class illuminaChannel that defines the summarisation to be performed. The default is to do a log2 transformation and report the mean and standard deviation of each bead type after outlier removal.
   }

   \item{probeIDs}{
	Vector of ArrayAddressIDs to be included in the summarized object. The default is to summarize all probes.
	}
  \item{useSampleFac}{
	if TRUE sections belonging to the same biological sample will be combined. The default is to summarize each section separately.
}

  \item{sampleFac}{
	optional character vector giving which a sample identifer for each section

}


 	
   \item{weightNames}{
	name of column in the \code{beadLevelData} to take extract weights 
}


\item{removeUnMappedProbes}{
if TRUE and annotation information is stored in the \code{beadLevelData} object, any ArrayAddressIDs that cannot be mapped to ILMN IDs will be removed.
}


}


\details{
From beadarray version 2.0 onwards, users are allowed more flexibility in how to create summarized data from bead-level data. The \code{illuminaChannel} is a means of allowing this flexibility by definining how summarization will be performed on each array section in the bead-level data object. The three keys steps applied to each section are; 1) use a transform function to get the quantities to be summarized (one value per bead). The most common use-case would be to extract the Green channel intensities and possibly perform a log2 transformation. 2) remove any outliers from this list of values 3) split the values according to ArrayAddressIDs and apply the definied \code{exprFun} and \code{varFun} to the quantities belonging to each ArrayAddress.

Some Illumina chips have multiple sections for the same biological sample; for example the HumanWG-6 chip or the Infinium genotyping chips. For such cases it may be more convenient to produce a summarized object where each column in the output is a different biological sample. This is especially important for genotpying chips where different SNPs are interrogated on the different sections, making a section-based summary problematic.

If the \code{useSampleFac} argument is set to TRUE, beadarray will try and combine sections belonging to the same sample. If the location of the sdf file for the chip is sucessfully stored in the \code{experimentData} slot of the \code{beadLevelData} object, the sdf will be interrogated to determine how samples were allocated to the chip. Otherwise the user can specifiy a sample factor that is the same length as the number of sections. If the sample factor is not supplied, or cannot be determined, then beadarray will summarize each section separately. 

During the course of the summary, ArrayAddressIDs present in the \code{beadLevelData} object will be converted to Illumina IDs (prefix ILMN) if the annotation of the object was set by \code{readIllumina} or \code{setAnnotation}. The rownames of the resulting ExpressionSetIllumina will be set to these new IDs, and the \code{featureData} slot will contain the original and new IDs. Any control probes present in the \code{beadLevelData} object will retain their original ArrayAddressID and the Status vector in \code{featureData} will report if each probe is a control or regular probe. Some ArrayAddressIDs present in the \code{beadLevelData} object may be neither regular probes will ILMN IDs, or control probes. These are internal controls used by Illumina and can be stopped from appearing in the summarized object by choosing the \code{removeUnmappedProbes = TRUE} option. 

The user can specify a vector of ArrayAddressIDs to be summarized using the \code{probeIDs} argument. Otherwise, a unique set of IDs is derived from all the array sections in the \code{beadLevelData} object.


}


\value{
    Returns an object of class ExpressionSetIllumina
}

\author{
    Mark Dunning
}

\examples{

if(require(beadarrayExampleData)){

	data(exampleBLData)

	bsd = summarize(exampleBLData)


	myMean = function(x) mean(x,na.rm=TRUE)
	mySd = function(x) sd(x,na.rm=TRUE)

	greenChannel = new("illuminaChannel", 
	        logGreenChannelTransform, 
	        illuminaOutlierMethod, 
	        myMean, mySd, "G"
	   )

	bsd = summarize(exampleBLData, channelList = list(greenChannel))

	bsd


}

}
	


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}

