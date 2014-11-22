\name{getSplitPositionForGroup}
\alias{getSplitPositionForGroup}
\title{
  Obtain split-read positions for different groups

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSplitPositionForGroup(dirBamFile, listFile = NULL, windows = 500, chr = NULL, medLeft = NULL, medRight = NULL, qualityThreshold = 0, epsilonOpen = NULL, usingPairedEnds = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dirBamFile}{
%%     ~~Describe \code{dirBamFile} here~~
}
  \item{listFile}{
%%     ~~Describe \code{listFile} here~~
}
  \item{windows}{
%%     ~~Describe \code{windows} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{medLeft}{
%%     ~~Describe \code{medLeft} here~~
}
  \item{medRight}{
%%     ~~Describe \code{medRight} here~~
}
  \item{qualityThreshold}{
%%     ~~Describe \code{qualityThreshold} here~~
}
  \item{epsilonOpen}{
%%     ~~Describe \code{epsilonOpen} here~~
}
  \item{usingPairedEnds}{
%%     ~~Describe \code{usingPairedEnds} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dirBamFile, listFile = NULL, windows = 500, chr = NULL, 
    medLeft = NULL, medRight = NULL, qualityThreshold = 0, epsilonOpen = NULL, 
    usingPairedEnds = TRUE) 
{
    if (is.na(dirBamFile)) 
        dirBamFile <- "./"
    if (substr(dirBamFile, length(dirBamFile), 1) != "/") 
        dirBamFile <- paste(dirBamFile, "/", sep = "")
    if (is.null(listFile)) 
        stop("listFile is not NULL")
    if (is.null(windows)) 
        stop("Please input windows used in the read-depth based approach")
    if (is.null(chr)) 
        stop("Please input chr")
    if (is.null(medLeft)) 
        stop("Please input medLeft obtained from the read-depth based approach")
    if (is.null(medRight)) 
        stop("Please input medRight obtained from the read-depth based approach")
    if (is.null(epsilonOpen)) 
        epsilonOpen <- 2 * windows
    chr <- gsub("chr", "", chr)
    what <- c("pos", "cigar", "mapq")
    which <- RangesList(`2` = IRanges(c(medLeft - epsilonOpen, 
        medRight - epsilonOpen), c(medLeft + epsilonOpen, medRight + 
        epsilonOpen)))
    if ((medLeft + epsilonOpen) >= (medRight - epsilonOpen)) {
        which <- RangesList(`2` = IRanges(c(medLeft - epsilonOpen, 
            medLeft + epsilonOpen + 1), c(medLeft + epsilonOpen, 
            medRight + epsilonOpen)))
    }
    names(which) <- as.character(as.name(chr))
    param <- ScanBamParam(what = what, which = which)
    tempLeft <- tempRight <- NULL
    for (ii in 1:length(listFile)) {
        bam <- scanBam(paste(dirBamFile, listFile[ii], sep = ""), 
            param = param)
        allPos <- rbind(do.call(cbind, bam[[1]]), do.call(cbind, 
            bam[[2]]))
        allPos <- allPos[as.numeric(allPos[, 2]) >= qualityThreshold, 
            ]
        if (!is.matrix(allPos)) 
            allPos <- data.frame(allPos[1], allPos[3], stringsAsFactors = FALSE)
        else allPos <- allPos[, -2]
        allPos <- allPos[grep("S", allPos[, 2]), ]
        if (!is.matrix(allPos)) 
            allPos <- data.frame(allPos[1], allPos[2], stringsAsFactors = FALSE)
        allPos <- allPos[grep("D|I", allPos[, 2], invert = TRUE), 
            ]
        leftPos1 <- rightPos1 <- NULL
        if (!is.matrix(allPos)) 
            allPos <- data.frame(allPos[1], allPos[2], stringsAsFactors = FALSE)
        if (!is.null(allPos)) {
            if (dim(allPos)[1] > 0) {
                leftList <- strsplit(allPos[, 2], "M")
                leftListTemp <- sapply(leftList, function(x) length(x) == 
                  2)
                leftList <- leftList[leftListTemp]
                leftPosTemp <- do.call(rbind, leftList)
                leftPos1 <- data.frame(allPos[, 1][leftListTemp], 
                  leftPosTemp)
                if (dim(leftPos1)[1] > 0) 
                  leftPos1 <- leftPos1[grep("S|D|I", leftPos1[, 
                    2], invert = TRUE), ]
                if (dim(leftPos1)[1] > 0) {
                  leftPos1 <- as.numeric(as.character(leftPos1[, 
                    1])) + as.numeric(as.character(leftPos1[, 
                    2]))
                }
                else leftPos1 <- NULL
                rightList <- strsplit(allPos[, 2], "S")
                rightListTemp <- sapply(rightList, function(x) length(x) == 
                  2)
                rightList <- rightList[rightListTemp]
                rightPosTemp <- do.call(rbind, rightList)
                rightPos1 <- data.frame(allPos[, 1][rightListTemp], 
                  rightPosTemp)
                if (dim(rightPos1)[1] > 0) 
                  rightPos1 <- rightPos1[grep("S|M|D|I", rightPos1[, 
                    2], invert = TRUE), ]
                if (dim(rightPos1)[1] > 0) 
                  rightPos1 <- rightPos1[grep("S|D|I", rightPos1[, 
                    3], invert = TRUE), ]
                if (dim(rightPos1)[1] > 0) {
                  rightPos1 <- as.numeric(as.character(rightPos1[, 
                    1]))
                }
                else rightPos1 <- NULL
            }
        }
        message("tempLeft")
        print(length(tempLeft))
        message("leftPos1")
        print(length(leftPos1))
        tempLeft <- c(tempLeft, leftPos1)
        message("tempRight")
        print(length(tempRight))
        message("rightPos1")
        print(length(rightPos1))
        tempRight <- c(tempRight, rightPos1)
    }
    return(list(tempLeft = tempLeft, tempRight = tempRight))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
