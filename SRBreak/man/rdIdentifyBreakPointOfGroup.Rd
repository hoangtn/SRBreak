\name{rdIdentifyBreakPointOfGroup}
\alias{rdIdentifyBreakPointOfGroup}
\title{
  Identify breakpoints for different groups using read-depth information
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rdIdentifyBreakPointOfGroup(dataMatrix, classM, upperCNThreshold = 0.4, lowerCNThreshold = -0.4, windows = 500, sigMaTemp = windows/3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrix}{
%%     ~~Describe \code{dataMatrix} here~~
}
  \item{classM}{
%%     ~~Describe \code{classM} here~~
}
  \item{upperCNThreshold}{
%%     ~~Describe \code{upperCNThreshold} here~~
}
  \item{lowerCNThreshold}{
%%     ~~Describe \code{lowerCNThreshold} here~~
}
  \item{windows}{
%%     ~~Describe \code{windows} here~~
}
  \item{sigMaTemp}{
%%     ~~Describe \code{sigMaTemp} here~~
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
function (dataMatrix, classM, upperCNThreshold = 0.4, lowerCNThreshold = -0.4, 
    windows = 500, sigMaTemp = windows/3) 
{
    outData <- NULL
    for (ii in 1:length(table(classM))) {
        subER1 <- cbind(classM, dataMatrix[pmatch(names(classM), 
            rownames(dataMatrix)), ])
        tempData <- subER1[subER1[, 1] == ii, ]
        if (is.null(dim(tempData))) {
            tempData <- tempData[-1]
            tempData <- matrix(tempData, nrow = 1)
            rownames(tempData) <- names(classM[classM == ii])
        }
        else tempData <- tempData[, -1]
        scorePos <- NULL
        tempScore <- NULL
        nSample <- dim(tempData)[1]
        tempTakeScore <- matrix(0, ncol = 2, nrow = nSample)
        for (k1 in 1:nSample) {
            tempData1 <- tempData[sample(1:nSample, nSample, 
                replace = TRUE), ]
            if (!is.matrix(tempData1)) 
                tempData1 <- matrix(tempData1, ncol = ncol(tempData))
            for (kk in 1:length(tempPos)) {
                tempScore[kk] <- sum(abs(tempData1[, kk + 1] - 
                  tempData1[, kk]))
            }
            names(tempScore) <- tempPos
            for (kk in 1:length(tempPos)) {
                scorePos[kk] <- sum(tempScore * dnorm(x = tempPos, 
                  mean = tempPos[kk], sd = sigMaTemp))
            }
            names(scorePos) <- tempPos
            tempOutScoreAA <- sort(scorePos, decreasing = TRUE)[1:2]
            tempOutScoreAA <- sort(as.numeric(names(tempOutScoreAA)))
            tempTakeScore[k1, ] <- tempOutScoreAA
        }
        aa <- rownames(tempData)
        breakSout <- apply(tempTakeScore, 2, median)
        tempTakeScoreL <- sort(tempTakeScore[, 1])
        tempTakeScoreR <- sort(tempTakeScore[, 2])
        breakSout[1] <- tempTakeScoreL[floor((length(tempTakeScore) + 
            1)/2)]
        breakSout[2] <- tempTakeScoreR[floor((length(tempTakeScore) + 
            1)/2)]
        leftPos <- as.numeric(rownames(mSD3[mSD3[, 1] == breakSout[1], 
            ]))
        rightPos <- as.numeric(rownames(mSD3[mSD3[, 2] == breakSout[2], 
            ]))
        groupSubMatrix <- subRegionMatrix[, leftPos:rightPos]
        if (is.matrix(groupSubMatrix)) {
            groupSubMatrix <- groupSubMatrix[pmatch(aa, rownames(groupSubMatrix)), 
                ]
        }
        else {
            groupSubMatrix <- groupSubMatrix[pmatch(aa, names(groupSubMatrix))]
        }
        if (is.matrix(groupSubMatrix)) {
            scoreOfGroup <- apply(groupSubMatrix, 1, mean)
        }
        else {
            scoreOfGroup <- mean(groupSubMatrix)
        }
        CNStatus <- ifelse(scoreOfGroup > upperCNThreshold, "DUP", 
            ifelse(scoreOfGroup < lowerCNThreshold, "DEL", "NORMAL"))
        tempOut <- data.frame(aa, rep(breakSout[1], length(aa)), 
            rep(breakSout[2], length(aa)), scoreOfGroup, CNStatus)
        outData <- rbind(outData, tempOut)
    }
    colnames(outData) <- c("Name", "Start", "End", "Score", "Status")
    return(outData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
