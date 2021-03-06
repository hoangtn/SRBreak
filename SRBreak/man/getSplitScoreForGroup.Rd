\name{getSplitScoreForGroup}
\alias{getSplitScoreForGroup}
\title{
  Obtain kernel-based scores for split-read positions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSplitScoreForGroup(dirBamFile, listFile = NULL, typeSV = c("DUP", "DEL"), windows = 500, chr = NULL, medLeft = NULL, medRight = NULL, qualityThreshold = 0, epsilonOpen = NULL, sdSplitRead = 0.5, usingPairedEnds = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dirBamFile}{
%%     ~~Describe \code{dirBamFile} here~~
}
  \item{listFile}{
%%     ~~Describe \code{listFile} here~~
}
  \item{typeSV}{
%%     ~~Describe \code{typeSV} here~~
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
  \item{sdSplitRead}{
%%     ~~Describe \code{sdSplitRead} here~~
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
function (dirBamFile, listFile = NULL, typeSV = c("DUP", "DEL"), 
    windows = 500, chr = NULL, medLeft = NULL, medRight = NULL, 
    qualityThreshold = 0, epsilonOpen = NULL, sdSplitRead = 0.5, 
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
    typeSV <- match.arg(typeSV)
    outSplitPosition <- getSplitPositionForGroup(dirBamFile = dirBamFile, 
        listFile = listFile, chr = chr, medLeft = medLeft, medRight = medRight, 
        windows = windows, qualityThreshold = qualityThreshold, 
        epsilonOpen = epsilonOpen, usingPairedEnds = usingPairedEnds)
    leftScore <- rightScore <- NULL
    leftPos <- outSplitPosition$tempLeft
    rightPos <- outSplitPosition$tempRight
    if (typeSV == "DUP") {
        leftPos <- rightPos <- c(leftPos, rightPos)
    }
    leftPos <- leftPos[abs(leftPos - medLeft) <= epsilonOpen]
    rightPos <- rightPos[abs(rightPos - medRight) <= epsilonOpen]
    leftPos <- table(leftPos)
    print(leftPos)
    rightPos <- table(rightPos)
    print(rightPos)
    for (kk in 1:length(leftPos)) {
        leftScore[kk] <- sum(leftPos * dnorm(x = as.numeric(names(leftPos)), 
            mean = as.numeric(names(leftPos))[kk], sd = sdSplitRead))
    }
    names(leftScore) <- names(leftPos)
    for (kk in 1:length(rightPos)) {
        rightScore[kk] <- sum(rightPos * dnorm(x = as.numeric(names(rightPos)), 
            mean = as.numeric(names(rightPos))[kk], sd = sdSplitRead))
    }
    names(rightScore) <- names(rightPos)
    leftBreakFromSR <- as.numeric(names(sort(leftScore, decreasing = TRUE)[1]))
    rightBreakFromSR <- as.numeric(names(sort(rightScore, decreasing = TRUE)[1]))
    if (length(leftBreakFromSR) == 0) 
        leftBreakFromSR <- medLeft
    if (length(rightBreakFromSR) == 0) 
        rightBreakFromSR <- medRight
    return(list(splitBreak = c(leftBreakFromSR, rightBreakFromSR), 
        leftPos = sort(leftPos), rightPos = sort(rightPos), leftScore = sort(leftScore), 
        rightScore = sort(rightScore)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
