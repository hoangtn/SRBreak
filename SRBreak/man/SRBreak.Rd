\name{SRBreak}
\alias{SRBreak}
\title{
Identify copy-number variable regions
}
\description{
Identify copy-number varialbe regions
}
\usage{
SRBreak(readDepthWindow = 500, chr = NULL, st = NULL, en = NULL, dirBamFile = NULL, genes = NULL, geneNames = NULL, rdQualityMapping = 10, correctGC = TRUE, byGCcontent = 1, useRSamtoolsToCount = FALSE, byMAPPABILITYcontent = 1, mappabilityFile = NULL, detectAllRegion = TRUE, quantileThresholdOfSD = 0.85, lowerCNThreshold = -0.3, upperCNThreshold = 0.3, countThreshold = 5, sigMaTemp = windows/3, NTimes = 50, testType = c("SD", "Count", "positiveCount", "negativeCount"), pemMappingQuality = 0, epsilonPairedOpen = NULL, thresholdOfIntersectionBetweenRDandPEM = 0.7, splitreadMappingQuality = 0, epsilonSplitReadOpen = 160, sdSplitRead = 0.5, usingPairedEnds = TRUE)
}

\arguments{
  \item{readDepthWindow}{
    window sizes for the read-depth based approach
}
\item{chr}{
  name of the chromosome 
}
\item{st}{
  start position of the region
}
\item{en}{
  end position of the region
}
\item{dirBamFile}{
  the directory of BAM files
}
\item{genes}{
  coordinates of genes

}
\item{geneNames}{
    names of genes
}
\item{rdQualityMapping}{
  minimum mappaing quality of the read-depth based approach

}
\item{correctGC}{
  correct GC content or not (default: TRUE)

}
\item{byGCcontent}{
  the percentage of GC content will be corrected

}
\item{useRSamtoolsToCount}{
  use RSamtools to count reads (default: FALSE)
}
\item{byMAPPABILITYcontent}{
  the percentage of mappability content will be corrected
}
\item{mappabilityFile}{
  name of a mappability file
}
\item{detectAllRegion}{
  detect all copy-number variable regions (default: TRUE)
}
\item{quantileThresholdOfSD}{
  a threshold to call CNVRs using standard-deviation information
}
\item{lowerCNThreshold}{
  a threshold to call deletion events
}
\item{upperCNThreshold}{
  a threshold to call duplication events
}
\item{countThreshold}{
  a number to call CNVR

}
\item{sigMaTemp}{
  sigma value for the read-depth based kernel
}
\item{NTimes}{
  the number of times to resample to identify boundaries of a CNVR
}
\item{testType}{
  a test type to call CNVR

}
\item{pemMappingQuality}{
  minimum mapping quality to read paired-end reads
}
  \item{epsilonPairedOpen}{

}
  \item{thresholdOfIntersectionBetweenRDandPEM}{

}
  \item{splitreadMappingQuality}{

}
  \item{epsilonSplitReadOpen}{

}
  \item{sdSplitRead}{

}
  \item{usingPairedEnds}{

}
}
\details{

}
\value{

}
\references{

}
\author{
}
\note{

}


\seealso{

}
