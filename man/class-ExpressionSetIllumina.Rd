\name{ExpressionSetIllumina-class}
\docType{class}
\alias{ExpressionSetIllumina-class}
\alias{[,ExpressionSetIllumina-method}
\alias{dim,ExpressionSetIllumina-method}
\alias{exprs,ExpressionSetIllumina-method}
\alias{se.exprs,ExpressionSetIllumina-method}

\alias{dim,ExpressionSetIllumina-method}
\alias{exprs<-}
\alias{exprs<-,ExpressionSetIllumina,matrix-method}

\alias{se.exprs<-}
\alias{se.exprs<-,ExpressionSetIllumina,matrix-method}

\alias{qcData}
\alias{qcData,ExpressionSetIllumina-method}

\alias{show,ExpressionSetIllumina-method}
\alias{nObservations}
\alias{nObservations,ExpressionSetIllumina-method}
\alias{nObservations<-}
\alias{nObservations<-,ExpressionSetIllumina,matrix-method}

\alias{Detection}
\alias{Detection,ExpressionSetIllumina-method}
\alias{Detection<-}
\alias{Detection<-,ExpressionSetIllumina,matrix-method}


\alias{combine,ExpressionSetIllumina,ExpressionSetIllumina-method}

\title{Class "ExpressionSetIllumina"}

\description{
  Container for high-throughput assays and experimental
  metadata. \code{ExpressionSetIllumina} class is derived from
  \code{\link[Biobase:class.eSet]{eSet}}, and requires matrices \code{exprs},
  \code{se.exprs}, \code{nObservations}, \code{Detection} as assay data members. The slots \code{featureData}, \code{phenoData} are accessed in the usual manner using \code{fData} and \code{pData} functions.
	
  For \code{ExpressionSetIllumina} objects created from bead-level data (using the \code{summarize} function), a \code{QC} slot is used to contain any quality control data that was present in the \code{beadLevelData} object. This is a change from previous versions of beadarray, where the intensities of the control probes themselves were stored in this slot. From version 2.0.0 onwards, control probes are stored in the \code{assayData} slot with the regular probes and the \code{featureData} slot has a reference for which rows correspond to controls.

  The \code{ExpressionSetIllumina} class is able to accomodate different channels when created from bead-level data. The \code{channelNames} function may be used to find out what channels are present in the object. The \code{channel} function can be used to select a particular channel, returning an \code{ExpressionSetIllumina} object. 
  

}


\author{Mark Dunning}

\keyword{classes}
