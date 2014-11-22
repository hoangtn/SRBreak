\name{getPairedEndPositionForGroup}
\alias{getPairedEndPositionForGroup}
\title{
  Read paired-end reads into R
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPairedEndPositionForGroup(dirBamFile, listFile = NULL, windows = 500, chr = NULL, medLeft = NULL, medRight = NULL, qualityThreshold = 0, epsilonPairedOpen = NULL, thresholdOfIntersection = 0.9, typeSV = "DEL")
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
  \item{epsilonPairedOpen}{
%%     ~~Describe \code{epsilonPairedOpen} here~~
}
  \item{thresholdOfIntersection}{
%%     ~~Describe \code{thresholdOfIntersection} here~~
}
  \item{typeSV}{
%%     ~~Describe \code{typeSV} here~~
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
    medLeft = NULL, medRight = NULL, qualityThreshold = 0, epsilonPairedOpen = NULL, 
    thresholdOfIntersection = 0.9, typeSV = "DEL") 
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
    if (is.null(epsilonPairedOpen)) 
        epsilonPairedOpen <- 5 * windows
    chr <- gsub("chr", "", chr)
    what <- c("pos", "mpos", "mapq")
    which <- IRanges::RangesList(`2` = IRanges(medLeft - epsilonPairedOpen, 
        medLeft + epsilonPairedOpen))
    names(which) <- as.character(as.name(chr))
    param <- Rsamtools::ScanBamParam(what = what, which = which, 
        flag = scanBamFlag(isPaired = TRUE, isFirstMateRead = TRUE))
    tempLeft <- medLeft
    tempRight <- medRight
    allPos <- lapply(listFile, function(XFile) {
        bam <- Rsamtools::scanBam(paste(dirBamFile, XFile, sep = ""), 
            param = param)
        bamT <- data.frame(bam[[1]]$pos, bam[[1]]$mpos, bam[[1]]$mapq)
        bamT <- bamT[complete.cases(bamT), ]
        bamT <- bamT[bamT[, 3] >= qualityThreshold, ]
        return(bamT[, c(1, 2)])
    })
    allPos <- do.call(rbind, allPos)
    if (typeSV == "DEL") 
        allPos <- allPos[allPos[, 2] > allPos[, 1], ]
    else allPos <- allPos[allPos[, 2] < allPos[, 1], ]
    allPos <- apply(allPos, 1, function(x) {
        if (x[2] < x[1]) {
            xT <- x[1]
            x[1] <- x[2]
            x[2] <- xT
        }
        return(x)
    })
    allPos <- t(allPos)
    tempTest <- IRanges::IRanges(medLeft, medRight)
    tempTestPosition <- apply(allPos, 1, function(x) {
        returnValue <- FALSE
        t10 <- IRanges::IRanges(x[1], x[2])
        t1 <- IRanges::reduce(intersect(t10, tempTest))
        if (length(t1) > 0) {
            testT1 <- width(t1)/width(tempTest)
            testT2 <- width(t1)/width(t10)
            returnValue <- (testT1 >= thresholdOfIntersection) & 
                (testT2 >= thresholdOfIntersection)
        }
        return(returnValue)
    })
    allPos <- allPos[tempTestPosition, ]
    if (!is.matrix(allPos)) 
        allPos <- matrix(allPos, nrow = 1)
    if (!is.null(dim(allPos)) & dim(allPos)[1] > 0) {
        tempLeftVector <- abs(allPos[, 1] - medLeft)
        if (typeSV == "DEL") {
            tempLeft <- round(quantile(allPos[, 1], 0.975), 0)
            tempRight <- round(quantile(allPos[, 2], 0.025), 
                0)
        }
        else {
            tempLeft <- round(quantile(allPos[, 1], 0.005), 0)
            tempRight <- round(quantile(allPos[, 2], 0.975), 
                0)
        }
    }
    return(c(tempLeft, tempRight))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
