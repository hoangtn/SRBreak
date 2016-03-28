SRBreak <- function(readDepthWindow = 500,
                    chr = NULL, st  = NULL, en = NULL, chrNameForSplitRead = NULL,
                    dirBamFile = NULL, genes = NULL, geneNames = NULL, dirCoordinate = NULL,
                    rdQualityMapping = 0, correctGC = TRUE, byGCcontent= 1, useRSamtoolsToCount = FALSE,
                    byMAPPABILITYcontent = 1, mappabilityFile = NULL, useMixtureModel2ClusterGroup = FALSE,

                    detectAllRegion = FALSE, quantileThresholdOfSD = 0.85, lowerCNThreshold = -0.5,
                    upperCNThreshold = 0.5, countThreshold = NULL, 
                    sigMaTemp = readDepthWindow/3, NTimes = 50, testType = c("Count", "SD", "positiveCount", "negativeCount"),
                    
                    pemMappingQuality = 0,  epsilonPairedOpen = NULL, thresholdOfIntersectionBetweenRDandPEM = 0.8,

                    splitreadMappingQuality = 0, epsilonSplitReadOpen = 2*readDepthWindow,
                    sdSplitRead = 0.5, usingPairedEnds = TRUE,
                    singleSample = FALSE,
                    nIncreaseSampleSize = NULL , minLengthSV = 5000,
                    inputRawReadCountMatrix = NULL,
                    epsilonCovDET = -1,
                    NTimesThreshold = 20,
                    NtransferToOtherPackage = 20000,
			referenceGenome = "BSgenome.Hsapiens.UCSC.hg19",
                    reference_fasta = NULL,
                    printOut = FALSE, gcSmallThreshold = 0.001, nCore = 1,
                    segmentalDuplicationFile = NULL, adjustStartPosition = FALSE
                    ){


#####################RD approach###########################################################
###########################################################################################
############################################################################################

    ###Make a pseudo-gene
    if (is.null(genes))
        genes <- c((st + en)/2, (st + en)/2 + 1000)
    if (is.null(geneNames))
        geneNames <- "gene"
    

    if (is.null(dirCoordinate)){
        dirCoordinate <- "TempAll"
        warning("No input for dirCoordinate")
    }
    if (is.null(dirBamFile)){
        dirBamFile <- "."
        warning("No input for dirBamFile")
    }

    dir.create(dirCoordinate)
    objectCNVrd2 <- new("CNVrd2", windows = readDepthWindow, chr = chr,
                        st = st, en = en, dirBamFile = dirBamFile,
                        dirCoordinate = dirCoordinate,
                        genes = genes, geneNames = geneNames)

        #############################Obtain regions with no GC bases
      chr <- as.character(chr)

    if(is.null(reference_fasta)){

          st1 <- seq(st, en, by = windows)
                        
                        if (max(st1) < en)
                          st1 <- c(st1, en)
                        st1[-c(1, length(st1))] <- st1[-c(1, length(st1))] - 1
                        
                        gr1 <- GRanges(chr, IRanges(start = st1[-length(st1)], end = st1[-1]))
                        system.time(b1 <- getSeq(Hsapiens, gr1)) 
                        
                        gcContentInSegment <- apply(letterFrequency(b1, letters = c("G", "C")), 1, sum)
                        
                        gcContentInSegment <- ifelse(is.na(gcContentInSegment), 0, gcContentInSegment)
                        
                      } else{
                          names(referenceGenome) <- toupper(names(referenceGenome))

                          message("names(referenceGenome) ")
                          print(names(referenceGenome))
                          
                          positionChr <- grep(toupper(chr), names(referenceGenome))[1]

                          message("position: ")
                          print(positionChr)
                          tempG <- referenceGenome[positionChr]
                          
                          st1 <- seq(st, en, by = windows)
                          
                          if (max(st1) < en)
                            st1 <- c(st1, en)
                          st1[-c(1, length(st1))] <- st1[-c(1, length(st1))] - 1
                          
                          gr1 <- GRanges(chr, IRanges(start = st1[-length(st1)], end = st1[-1]))
                          system.time(b1 <- getSeq(tempG, gr1)) 
                          
                          gcContentInSegment <- apply(letterFrequency(b1, letters = c("G", "C")), 1, sum)
                          
                          gcContentInSegment <- ifelse(is.na(gcContentInSegment), 0, gcContentInSegment)
                          
                          
                          
    }
    gcContentInSegment <- gcContentInSegment/windows
    bZero <- which(gcContentInSegment < gcSmallThreshold )

    gr0 <- IRanges(start = st1[-length(st1)], end = st1[-1])
    gr2 <- gr0[bZero, ]
    grZeroOut <- reduce(gr2)

    newST1 <- start(grZeroOut)[1]
    newST2 <- end(grZeroOut)[1]
    rm(gr0)
    rm(gr2)
    rm(bZero)
    rm(gcContentInSegment)

    if (adjustStartPosition){
    if ((newST1 <= objectCNVrd2@st) & (newST2 >= objectCNVrd2@st)) {
        objectCNVrd2@st <- newST2
        message("Start position is changed to: ", newST2)
    }
    }


    #############Use CNVrd2 to count read
    if (!is.null(inputRawReadCountMatrix)){
        rawcntMatrix0 <- inputRawReadCountMatrix
    } else {
    rawcntMatrix0 <- CNVrd2::countReadInWindow(Object = objectCNVrd2,
                                       rawReadCount = TRUE, qualityThreshold = rdQualityMapping,
                                       correctGC = correctGC, byGCcontent= byGCcontent,
                                       useRSamtoolsToCount = useRSamtoolsToCount, 
					referenceGenome = referenceGenome,
					reference_fasta = reference_fasta, nCore = nCore)
}
    ############Correct mappability bias

    message("Correcting mappability bias")
    rawcntMatrix01 <- correctMappability(readCountMatrix = rawcntMatrix0,
                                     chr = objectCNVrd2@chr, start = objectCNVrd2@st,
                                     end = objectCNVrd2@en,
                                     mappabilityFile = mappabilityFile, 
                                     byMAPPABILITYcontent = byMAPPABILITYcontent  )


    #######For single sample, this step is aimed to make a pseu-do matrix of multiple samples (not good)
   
    if (!is.null(nIncreaseSampleSize)){
        tempNameSample <- rownames(rawcntMatrix01)


        newMatrixSample <- rawcntMatrix01

        for (ij1 in 2:nIncreaseSampleSize){
            tempRowSample <- rawcntMatrix01[1, ] + rnorm(length(rawcntMatrix01[1, ]), 0, 0.001)
            newMatrixSample <- rbind(newMatrixSample, tempRowSample)

            file.copy(from = paste(dirBamFile, tempNameSample, sep = ""), to = paste(dirBamFile, tempNameSample, ".", ij1,  "bam", sep = ""))
            file.copy(from = paste(dirBamFile, tempNameSample, ".bai", sep = ""), to = paste(dirBamFile, tempNameSample, ".", ij1, "bam.bai", sep = ""))
        }

        
        rownames(newMatrixSample) <- paste(tempNameSample, ".", 1:nIncreaseSampleSize, "bam", sep = "")
        rownames(newMatrixSample)[1] <- tempNameSample

        
        rawcntMatrix01 <- newMatrixSample

    }






    ###############################################################
    ###Transform data to use CNVrd2
    rawcntMatrix <- t(apply(rawcntMatrix01, 1, function(x){
        temp <- x
        if (median(x) > 0)
            temp <- x/median(x)
        temp <- temp - 1
        return(temp)}))

    if (is.null(countThreshold))
        countThreshold <- max(1, as.integer(0.075*dim(rawcntMatrix)[1]))
##Segment read counts into different regions
    resultSegment <- CNVrd2::segmentSamples(Object = objectCNVrd2,
                                    stdCntMatrix = rawcntMatrix)
    ############################################################
    ###Identify CNVRs
    polymorphicRegion <- CNVrd2::identifyPolymorphicRegion(Object = objectCNVrd2,
                   segmentObject = resultSegment,
                   plotPolymorphicRegion = FALSE)


    #####Identify CNVRs based on criteria
    testOut <- detectBreakPointFromRD(polymorphicObject = polymorphicRegion,
                                      genes = objectCNVrd2@genes,
                                      windows = objectCNVrd2@windows,
                                      detectAllRegion = detectAllRegion,
                                      quantileThreshold =  quantileThresholdOfSD,
                                      upperCNThreshold = upperCNThreshold,

                                      lowerCNThreshold = lowerCNThreshold,
                                      countThreshold = countThreshold,
                                      testType = testType,
                                      sigMaTemp = sigMaTemp, NTimes = NTimes,
                                      useMixtureModel2ClusterGroup = useMixtureModel2ClusterGroup,
                                      minLengthSV = minLengthSV,
                                      epsilonCovDET = epsilonCovDET,
                                      NTimesThreshold = NTimesThreshold,
                                      NtransferToOtherPackage = NtransferToOtherPackage,
                                      printOut = printOut,
                                      singleSample = singleSample)

############################################################################################
############################################################################################
###################Paired-end mapping and split-read sections
###############################################

    reportOutData <- NULL
    for (kk in 1:length(testOut)){

        testOutPairedEnd <- testOut[[kk]]
        ###########Use paired-end approach
        if (usingPairedEnds == TRUE){
            message("You're using paired-end information. It takes longer time to analyse")
            testOutPairedEnd <- detectBreakpointFromPairedEnds(resultFromRD = testOut[[kk]],
                                   dirBamFile = dirBamFile,
                                   windows = objectCNVrd2@windows, chr = objectCNVrd2@chr,
                                   qualityThreshold = pemMappingQuality,
                                   epsilonPairedOpen = epsilonPairedOpen,
                                   thresholdOfIntersection = thresholdOfIntersectionBetweenRDandPEM,
                                                               printOut = printOut)
        }
        

        ##########Use split-read approach
        finalOutAll <- getSplitReadBreakpoint(resultFromRD = testOutPairedEnd,
                                   dirBamFile = dirBamFile,
                                   windows = objectCNVrd2@windows, chr = objectCNVrd2@chr,
                                   qualityThreshold = splitreadMappingQuality,
                                   epsilonOpen = epsilonSplitReadOpen,
                                   sdSplitRead = sdSplitRead, usingPairedEnds = usingPairedEnds,
					chrNameForSplitRead = chrNameForSplitRead)

        #############Combine all outputs
        reportOutData <- rbind(reportOutData, finalOutAll)
        }

    ###########Produce outputs
    breakpointFinal <- IRanges::unique(reportOutData[, c(5, 7, 8)])
    breakpointFinal <- breakpointFinal[breakpointFinal[, 1] != "NORMAL", ]


    if (!is.matrix(breakpointFinal))
        breakpointFinal <- matrix(breakpointFinal, nrow = 1)
    
    
    breakpointDataFrame <- matrix("N", ncol = dim(rawcntMatrix)[1] + 3,
                              nrow = dim(breakpointFinal)[1])


    
    colnames(breakpointDataFrame) <- c("chr", "start", "end", rownames(rawcntMatrix))


    
    breakpointDataFrame[, c(2, 3)] <- breakpointFinal[, c(2, 3)]

    breakpointDataFrame[, 1] <- rep(objectCNVrd2@chr, dim(breakpointFinal)[1])

    shortReport <- reportOutData[, c(1, 5, 7, 8)]

    shortReport <- shortReport[shortReport[, 2] != "NORMAL", ]


    if (!is.matrix(shortReport))
        shortReport <- matrix(shortReport, nrow = 1)

    if (!is.null(shortReport) & dim(shortReport)[1] > 0){
        for (jj in 1:dim(breakpointDataFrame)[1]){
        tempShortData <- shortReport[(shortReport[, 3] == breakpointDataFrame[jj, 2])
                               & (shortReport[, 4] == breakpointDataFrame[jj, 3]),]

        if (is.data.frame(tempShortData) | is.matrix(tempShortData)){

            indexColumn <- pmatch(tempShortData[, 1], colnames(breakpointDataFrame))
            breakpointDataFrame[jj, indexColumn] <- tempShortData[, 2]
            } else {
                indexColumn <- pmatch(tempShortData[1], colnames(breakpointDataFrame))
                breakpointDataFrame[jj, indexColumn] <- tempShortData[2]

                }
        }}


    if (dim(breakpointDataFrame)[1] > 1)
        breakpointDataFrame <- breakpointDataFrame[order(breakpointDataFrame[, 2]),]

    breakpointDataFrame <- breakpointDataFrame[(as.numeric(breakpointDataFrame[, 3]) - as.numeric(breakpointDataFrame[, 2])) >= minLengthSV,]


    ##################Remove files#################################################

    ###############################################################################
    if (!is.null(nIncreaseSampleSize)){
    
        for (ij1 in 2:nIncreaseSampleSize){
     
            file.remove(paste(dirBamFile, tempNameSample, ".", ij1,  "bam", sep = ""))
            file.remove(paste(dirBamFile, tempNameSample, ".", ij1, "bam.bai", sep = ""))
        }

        breakpointDataFrame <- breakpointDataFrame[, c(1, 2, 3, 4)]
    }

############################################################################################
    ################################################################################
    ###Add segmental Duplication information
    outputCor <- IRanges(as.integer(breakpointDataFrame[, 2]), as.integer(breakpointDataFrame[, 3]))
    indexOut <- outputCor %outside% grZeroOut
    breakpointDataFrame <- breakpointDataFrame[indexOut, ]

    if (!is.null(segmentalDuplicationFile)){
         dupFile <- read.table(segmentalDuplicationFile)
         sDupA <- IRanges(as.integer(dupFile[, 2]), as.integer(dupFile[, 3]))
         sDupA <- reduce(sDupA)
         print(head(sDupA))
             outputCor <- IRanges(as.integer(breakpointDataFrame[, 2]), as.integer(breakpointDataFrame[, 3]))

    indexOut <- outputCor %outside% sDupA
    breakpointDataFrame <- breakpointDataFrame[indexOut, ]

    }



    
    return(list(svResult = breakpointDataFrame, objectSRBreak = objectCNVrd2,
           polymorphicRegionObject = polymorphicRegion, resultSegment = resultSegment))


}
                    

getSplitReadBreakpoint <- function(resultFromRD = NULL,
                                   dirBamFile = NULL,
                                   windows = 500, chr = NULL,
                                   qualityThreshold = 0,
                                   epsilonOpen = NULL,
                                   sdSplitRead = 0.5,
                                   usingPairedEnds = TRUE,
					chrNameForSplitRead = NULL){
    if (is.na(dirBamFile))
        dirBamFile <- "./"
    if (substr(dirBamFile, length(dirBamFile), 1) != "/")
        dirBamFile <- paste(dirBamFile, "/", sep = "")
    if (is.null(windows))
        stop("Please input windows used in the read-depth based approach")
    if (is.null(chr))
        stop("Please input chr")
    if (is.null(epsilonOpen))
        epsilonOpen <- 2*windows

    ############################################################
    normalGroup <- resultFromRD[resultFromRD[, dim(resultFromRD)[2] ] == "NORMAL", ]
    resultFromRD<- resultFromRD[resultFromRD[, dim(resultFromRD)[2] ] != "NORMAL", ]

    
    tempGroup <- rep(1, dim(resultFromRD)[1])
    nAA <- 2

    if (length(tempGroup) > 1){
    for (kk in 2:(dim(resultFromRD)[1])){
        if ((resultFromRD[kk, 2] != resultFromRD[kk - 1, 2]) |
            (resultFromRD[kk, 3] != resultFromRD[kk - 1, 3]) |
            (resultFromRD[kk, 5] != resultFromRD[kk - 1, 5])){
            tempGroup[kk:(length(tempGroup))] <- nAA
            nAA <- nAA + 1
        }
}}
    names(tempGroup) <- resultFromRD[, 1]
    groupFromRD <- as.numeric(names(table(tempGroup)))

###Group From read-depth approach
    groupTableFromRD <- cbind(resultFromRD, tempGroup)
    Nnormal <- dim(normalGroup)[1]
    finalResultSplitAndRead <- data.frame(normalGroup,
                                      rep("NORMAL", Nnormal),
                                      rep("NORMAL", Nnormal),
                                      rep("NORMAL", Nnormal))
    colnames(finalResultSplitAndRead) <- paste("V", 1:dim(finalResultSplitAndRead)[2], sep = "")
    tempDataFrameOut <- NULL

    for (gG in groupFromRD){
        subGroupFromRD <- names(tempGroup[tempGroup == gG])
        SubgroupTableFromRD <- groupTableFromRD[groupTableFromRD[, dim(groupTableFromRD)[2]] == gG, ]
        tempSplitOut <- getSplitScoreForGroup(dirBamFile = dirBamFile,
                                              listFile = SubgroupTableFromRD[, 1],
                                              chr = chr, windows = windows,
                                              typeSV = toupper(as.character(SubgroupTableFromRD[, 5])[1]),
                                              medLeft =  SubgroupTableFromRD[, 2][1],
                                              medRight = SubgroupTableFromRD[, 3][1],
                                              sdSplitRead = 0.5,
                                              epsilonOpen = epsilonOpen,
                                              usingPairedEnds = usingPairedEnds,
						chrNameForSplitRead = chrNameForSplitRead)
        
        tempDataFrameOut <- rbind(tempDataFrameOut, data.frame(SubgroupTableFromRD,
                                                               rep(tempSplitOut$splitBreak[1],
                                                                   dim(SubgroupTableFromRD)[1]),
                                   rep(tempSplitOut$splitBreak[2], dim(SubgroupTableFromRD)[1])))
        }

    if (!is.null(tempDataFrameOut)){
    colnames(tempDataFrameOut) <- paste("V", 1:dim(finalResultSplitAndRead)[2], sep = "")
    ppOut <- rbind(apply(finalResultSplitAndRead, 2, as.character),
                   apply(tempDataFrameOut, 2, as.character))
} else ppOut <- apply(finalResultSplitAndRead, 2, as.character)


    ###################################
    return(ppOut)
       
}
######################################################################
getSplitScoreForGroup <- function(dirBamFile, listFile = NULL,
                                  typeSV = c("DUP", "DEL"),
                                  windows = 500, chr = NULL,
                                  medLeft = NULL, medRight = NULL,
                                  qualityThreshold = 0,
                                  epsilonOpen = NULL,
                                  sdSplitRead = 0.5,
                                  usingPairedEnds = TRUE,
				chrNameForSplitRead = NULL){
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
        epsilonOpen <- 2*windows

    typeSV <- match.arg(typeSV)

    #############Obtain split positions
    outSplitPosition <- getSplitPositionForGroup(dirBamFile = dirBamFile,
                                             listFile = listFile,
                                             chr = chr,
                                             medLeft =  medLeft,
                                             medRight = medRight,
                                                 windows = windows,
                                                 qualityThreshold = qualityThreshold,
                                                 epsilonOpen = epsilonOpen,
                                                 usingPairedEnds = usingPairedEnds,
						chrNameForSplitRead = chrNameForSplitRead)
    leftScore <- rightScore <- NULL
    leftPos <- outSplitPosition$tempLeft
    rightPos <- outSplitPosition$tempRight
    if (typeSV == "DUP"){
        leftPos <- rightPos <- c(leftPos, rightPos)
    }

    leftPos <- leftPos[abs(leftPos - medLeft) <= epsilonOpen]

    rightPos <- rightPos[abs(rightPos - medRight) <= epsilonOpen]

    leftPos <- table(leftPos)
    
    rightPos <- table(rightPos)

    
    ###########Calculate scores
    for (kk in 1:length(leftPos)){
        leftScore[kk] <- sum(leftPos*dnorm(x = as.numeric(names(leftPos)),
                mean = as.numeric(names(leftPos))[kk], sd = sdSplitRead))
        }
    names(leftScore) <- names(leftPos)

   for (kk in 1:length(rightPos)){
        rightScore[kk] <- sum(rightPos*dnorm(x = as.numeric(names(rightPos)),
                mean = as.numeric(names(rightPos))[kk], sd = sdSplitRead))
        }
    names(rightScore) <- names(rightPos)

    #######################

    leftBreakFromSR <- as.numeric(names(sort(leftScore, decreasing = TRUE)[1]))
    rightBreakFromSR <- as.numeric(names(sort(rightScore, decreasing = TRUE)[1]))

    if (length(leftBreakFromSR) == 0)
        leftBreakFromSR <- medLeft
    if (length(rightBreakFromSR) == 0)
        rightBreakFromSR <- medRight

    return(list(splitBreak = c(leftBreakFromSR, rightBreakFromSR),
                leftPos = sort(leftPos), rightPos = sort(rightPos),
           leftScore = sort(leftScore), rightScore = sort(rightScore)))
}


###############################################################################
###############################################################################

getSplitPositionForGroup <- function(dirBamFile, listFile = NULL, windows = 500, chr = NULL, medLeft = NULL,
                                     medRight = NULL, qualityThreshold = 0, epsilonOpen = NULL,
                                     usingPairedEnds = TRUE, chrNameForSplitRead = NULL){
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
        epsilonOpen <- 2*windows
    
    chr <- gsub("chr", "", chr)

   if (!is.null(chrNameForSplitRead))
     chr <- chrNameForSplitRead
	

    what <- c("pos", "cigar", "mapq")
    

    which <- RangesList('2' = IRanges(c(medLeft - epsilonOpen, medRight - epsilonOpen),
                            c(medLeft + epsilonOpen, medRight + epsilonOpen)))
    if ((medLeft + epsilonOpen) >= (medRight - epsilonOpen)){
        which <- RangesList('2' = IRanges(c(medLeft - epsilonOpen, medLeft + epsilonOpen + 1),
                                c(medLeft + epsilonOpen, medRight + epsilonOpen)))
    }
    names(which) <- as.character(as.name(chr))
    param <- ScanBamParam( what = what, which = which)
########Left/Right split positions#################################
    tempLeft <- tempRight <- NULL
#################Reading BAM files###############################
    for (ii in 1:length(listFile)){
        bam <- scanBam(paste(dirBamFile, listFile[ii], sep = ""),  param=param)
        allPos <- rbind(do.call(cbind, bam[[1]]), do.call(cbind, bam[[2]]))
        allPos <- allPos[as.numeric(allPos[, 2]) >= qualityThreshold,]

        if (!is.matrix(allPos))
            allPos <- data.frame(allPos[1], allPos[3], stringsAsFactors = FALSE)

        else allPos <- allPos[, -2]
        
        ##Search CIGAR strings including 'S'
        allPos <- allPos[grep("S", allPos[, 2]),]
        if (!is.matrix(allPos))
            allPos <- data.frame(allPos[1], allPos[2], stringsAsFactors = FALSE)

        allPos <- allPos[grep("D|I", allPos[, 2], invert  = TRUE),]
        leftPos1 <- rightPos1 <- NULL

        if (!is.matrix(allPos))
            allPos <- data.frame(allPos[1], allPos[2], stringsAsFactors = FALSE)

        
        if (!is.null(allPos)){
        if (dim(allPos)[1] > 0){
                 
############################################################
############Left positions##################################
##Search -M-S CIGAR strings
            leftList <- strsplit(allPos[, 2], "M")
            leftListTemp <- sapply(leftList, function(x) length(x) == 2)
            leftList <- leftList[leftListTemp]
            leftPosTemp <- do.call(rbind, leftList)
            leftPos1 <- data.frame(allPos[, 1][leftListTemp], leftPosTemp)
    
##Remove positons having 'S' in the first part
            if (dim(leftPos1)[1] > 0)
                leftPos1 <- leftPos1[grep("S|D|I", leftPos1[, 2], invert = TRUE),]
    
##Obtain left positions
            if (dim(leftPos1)[1] > 0){
                leftPos1 <- as.numeric(as.character(leftPos1[, 1])) + as.numeric(as.character(leftPos1[, 2])) - 1
                } else leftPos1 <- NULL
############################################################
############Right positions##################################
##Search -M-S CIGAR strings
            rightList <- strsplit(allPos[, 2], "S")
            rightListTemp <- sapply(rightList, function(x) length(x) == 2)
            rightList <- rightList[rightListTemp]
            rightPosTemp <- do.call(rbind, rightList)
            rightPos1 <- data.frame(allPos[, 1][rightListTemp], rightPosTemp)
                
##Remove positons having 'S' in the first part
            if (dim(rightPos1)[1] > 0)
                rightPos1 <- rightPos1[grep("S|M|D|I", rightPos1[, 2], invert = TRUE),]
            if (dim(rightPos1)[1] > 0)
                rightPos1 <- rightPos1[grep("S|D|I", rightPos1[, 3], invert = TRUE),]
               
##Obtain right positions
            if (dim(rightPos1)[1] > 0) {
                rightPos1 <- as.numeric(as.character(rightPos1[, 1])) - 1
                
                } else rightPos1 <- NULL
            }}


        tempLeft <- c(tempLeft, leftPos1)
        

        

        tempRight <- c(tempRight, rightPos1)

       
}
    return(list(tempLeft = tempLeft, tempRight = tempRight))
}


####################################################################
#####################NEW FUNCTION###################################

detectBreakpointFromPairedEnds <- function(resultFromRD = NULL,
                                   dirBamFile = NULL,
                                   windows = 500, chr = NULL,
                                   qualityThreshold = 0,
                                   epsilonPairedOpen = NULL,
                                           thresholdOfIntersection = 0.9,
                                           insertSize = 500, readLength = 100,
                                           printOut = FALSE, chrNameForSplitRead = NULL){

    if (printOut)
        message("Using paired-end information")
    
    if (is.na(dirBamFile))
        dirBamFile <- "./"
    if (substr(dirBamFile, length(dirBamFile), 1) != "/")
        dirBamFile <- paste(dirBamFile, "/", sep = "")
    if (is.null(windows))
        stop("Please input windows used in the read-depth based approach")
    if (is.null(chr))
        stop("Please input chr")
    if (is.null(epsilonPairedOpen))
        epsilonPairedOpen <- 2*windows

    ############################################################
    normalGroup <- resultFromRD[resultFromRD[, dim(resultFromRD)[2] ] == "NORMAL", ]
    resultFromRD<- resultFromRD[resultFromRD[, dim(resultFromRD)[2] ] != "NORMAL", ]
    tempGroup <- rep(1, dim(resultFromRD)[1])
    nAA <- 2

    if (length(tempGroup) > 1){
    for (kk in 2:(dim(resultFromRD)[1])){
        if ((resultFromRD[kk, 2] != resultFromRD[kk - 1, 2]) |
            (resultFromRD[kk, 3] != resultFromRD[kk - 1, 3])){
            tempGroup[kk:(length(tempGroup))] <- nAA
            nAA <- nAA + 1
        }
}}
    names(tempGroup) <- resultFromRD[, 1]
    groupFromRD <- as.numeric(names(table(tempGroup)))

###Group From read-depth approach
    groupTableFromRD <- cbind(resultFromRD, tempGroup)
    Nnormal <- dim(normalGroup)[1]


    tempDataFrameOut <- NULL

    for (gG in groupFromRD){

        subGroupFromRD <- names(tempGroup[tempGroup == gG])
        SubgroupTableFromRD <- groupTableFromRD[groupTableFromRD[, dim(groupTableFromRD)[2]] == gG, ]
        
        tempPairedEndOut <- getPairedEndPositionForGroup(dirBamFile = dirBamFile,
                                              listFile = SubgroupTableFromRD[, 1],
                                              chr = chr, windows = windows,
                                              medLeft =  SubgroupTableFromRD[, 2][1],
                                              medRight = SubgroupTableFromRD[, 3][1],
                                              thresholdOfIntersection = thresholdOfIntersection,
                                              typeSV = toupper(as.character(SubgroupTableFromRD[, 5])[1]),
						chrNameForSplitRead = chrNameForSplitRead)


        SubgroupTableFromRD$Start <- tempPairedEndOut[1]
        SubgroupTableFromRD$End <- tempPairedEndOut[2]
        tempDataFrameOut <- rbind(tempDataFrameOut, SubgroupTableFromRD)
        
        }

    #############################Merge Groups

    if (!is.null(tempDataFrameOut)){
        
    mergeTemp <- NULL
    mergeRemove <- tempDataFrameOut

    while (dim(mergeRemove)[1] >= 1){
        xTemp  <- mergeRemove[1, ]
        tempT1 <- mergeRemove[as.character(mergeRemove[, 5]) == as.character(xTemp[, 5]),]
        tempT1 <- tempT1[abs(tempT1[, 2] - xTemp[, 2]) <= (insertSize - 2*readLength),]
        tempT1 <- tempT1[abs(tempT1[, 3] - xTemp[, 3]) <= (insertSize - 2*readLength),]


        
        if (tempT1[1, 5] == "DEL"){
            tempT1[, 2] <- round(quantile(tempT1[, 2], 0.9))
            tempT1[, 3] <- round(quantile(tempT1[, 3], 0.1))
        } else {
            tempT1[, 2] <- round(quantile(tempT1[, 2], 0.1))
            tempT1[, 3] <- round(quantile(tempT1[, 3], 0.9))
        }
            
        mergeTemp <- rbind(mergeTemp, tempT1)

        mergeRemove <- mergeRemove[-pmatch(rownames(tempT1), rownames(mergeRemove)),]

    }
    #########################################

    tempDataFrameOut <- mergeTemp
    tempDataFrameOut <- tempDataFrameOut[, -dim(tempDataFrameOut)[2]]
    colnames(tempDataFrameOut) <- paste("V", 1:dim(tempDataFrameOut)[2], sep = "")
}    
    colnames(normalGroup) <- paste("V", 1:dim(normalGroup)[2], sep = "")

    ppOut <- rbind(tempDataFrameOut, normalGroup)
    

    

    rownames(ppOut) <- as.character(ppOut[, 1])

    ###################################
    return(ppOut)
       
}
######################################################################

###############################################################################
###########################NEW FUNCTION########################################
###############################################################################

getPairedEndPositionForGroup <- function(dirBamFile, listFile = NULL, windows = 500, chr = NULL,
                                         medLeft = NULL, medRight = NULL, qualityThreshold = 0,
                                         epsilonPairedOpen = NULL, thresholdOfIntersection = 0.9,
                                         typeSV = "DEL", chrNameForSplitRead = NULL){
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
        epsilonPairedOpen <- 5*windows
    
    chr <- gsub("chr", "", chr)

    what <- c("pos", "mpos", "mapq")
    if (!is.null(chrNameForSplitRead))
	chr <- chrNameForSplitRead    

    which <- IRanges::RangesList('2' = IRanges(medLeft - epsilonPairedOpen, medLeft + epsilonPairedOpen))
    names(which) <- as.character(as.name(chr))
    param <- Rsamtools::ScanBamParam( what = what, which = which,
                          flag=scanBamFlag(isPaired=TRUE, isFirstMateRead = TRUE))
########Left/Right split positions#################################
    tempLeft <- medLeft
    tempRight <- medRight
#################Reading BAM files###############################
    allPos <- lapply(listFile, function(XFile){
        bam <- Rsamtools::scanBam(paste(dirBamFile, XFile, sep = ""),  param=param)
        bamT <- data.frame(bam[[1]]$pos, bam[[1]]$mpos, bam[[1]]$mapq)
        
        bamT <- bamT[complete.cases(bamT), ]
        bamT <- bamT[bamT[, 3] >= qualityThreshold,]
        return(bamT[, c(1, 2)])})

    allPos <- do.call(rbind, allPos)

    if (typeSV == "DEL")
        allPos <- allPos[allPos[, 2] > allPos[, 1],]
    else
        allPos <- allPos[allPos[, 2] < allPos[, 1],]
    
    allPos <- apply(allPos,1, function(x) {
        if (x[2] < x[1]){
            xT <- x[1]
            x[1] <- x[2]
            x[2] <- xT}
        return(x)})

    allPos <- t(allPos)
    tempTest <- IRanges::IRanges(medLeft, medRight)
    tempTestPosition <- apply(allPos, 1, function(x){
        returnValue <- FALSE
        
        t10 <- IRanges::IRanges(x[1], x[2])
        t1 <- IRanges::reduce(IRanges::intersect(t10, tempTest))
        if (length(t1) > 0){
            testT1 <- width(t1)/width(tempTest)
            testT2 <- width(t1)/width(t10)

            returnValue <- (testT1 >= thresholdOfIntersection) & (testT2 >= thresholdOfIntersection)
        }
        return(returnValue)})
          
    
    allPos <- allPos[tempTestPosition, ]
    if (!is.matrix(allPos))
         allPos <- matrix(allPos, nrow = 1)
    ######################################Test left and right##############
        if (!is.null(dim(allPos)) & dim(allPos)[1] > 0){
                    tempLeftVector <- abs(allPos[, 1] - medLeft)

            if (typeSV == "DEL"){
                tempLeft <- round(quantile(allPos[, 1], 0.975), 0)
                tempRight <- round(quantile(allPos[, 2], 0.025), 0)
                } else {
                tempLeft <- round(quantile(allPos[, 1], 0.005), 0)
                tempRight <- round(quantile(allPos[, 2], 0.975), 0)
            }
}

            

    return(c(tempLeft, tempRight))
}



#############################################
############################################
###############NEW FUNCTION#################
############################################################################
############################NEW FUNCTION###################################        
rdIdentifyBreakPointOfGroup <- function(dataMatrix, classM,
                                        upperCNThreshold = 0.4,
                                        lowerCNThreshold = -0.4,
                                        windows = 500,
                                        sigMaTemp = windows/3, useMixtureModel2ClusterGroup = FALSE){
    outData <- NULL

    for (ii in 1:length(table(classM))){
        #Combine classes and the matrix
        subER1 <- cbind(classM, dataMatrix[pmatch(names(classM), rownames(dataMatrix)),])
        #Retain only a sub-matrix in class being considered
        tempData <- subER1[subER1[, 1] ==ii, ]
        #tempData <- matrix(tempData, ncol = dim(subER1)[2])
        #tempData <- matrix(tempData, ncol = length(tempData))
        if (is.null(dim(tempData))){
            tempData <- tempData[-1]
            tempData <- matrix(tempData, nrow = 1)
            rownames(tempData) <- names(classM[classM == ii])
            } else
                tempData <- tempData[ , -1]

        scorePos <- NULL
        tempScore <- NULL
        ###Table random sample
        nSample <- dim(tempData)[1]
        tempTakeScore <- matrix(0, ncol = 2, nrow = nSample)

        for (k1 in 1:nSample){
                    
            tempData1 <- tempData[sample(1:nSample, nSample, replace=TRUE),]
              if (!is.matrix(tempData1))
                tempData1 <- matrix(tempData1, ncol = ncol(tempData))
            for (kk in 1:length(tempPos)){
                tempScore[kk] <- sum(abs(tempData1[, kk + 1] - tempData1[, kk]))
                }
            names(tempScore) <- tempPos
            for (kk in 1:length(tempPos)){
                scorePos[kk] <- sum(tempScore*dnorm(x = tempPos, mean = tempPos[kk], sd = sigMaTemp))

                }
            names(scorePos) <- tempPos
            tempOutScoreAA <- sort(scorePos, decreasing = TRUE)[1:2]
            tempOutScoreAA <- sort(as.numeric(names(tempOutScoreAA)))
            tempTakeScore[k1, ] <- tempOutScoreAA
            }

        #aa <- substr(rownames(tempData), 1, 7)
        aa <- rownames(tempData)
        breakSout <- apply(tempTakeScore, 2, median)
        tempTakeScoreL <- sort(tempTakeScore[, 1])
        tempTakeScoreR <- sort(tempTakeScore[, 2])
        breakSout[1] <- tempTakeScoreL[floor((length(tempTakeScore) + 1)/2)]
        breakSout[2] <- tempTakeScoreR[floor((length(tempTakeScore) + 1)/2)]
        ###Obtain positions of two breakpoints
        leftPos <- as.numeric(rownames(mSD3[mSD3[, 1] == breakSout[1],]))
        rightPos <- as.numeric(rownames(mSD3[mSD3[, 2] == breakSout[2],]))
        groupSubMatrix <- subRegionMatrix[, leftPos:rightPos]

####################Take segmentation results
        if (is.matrix(groupSubMatrix)){
            groupSubMatrix <- groupSubMatrix[pmatch(aa, rownames(groupSubMatrix)),]
            } else{
    ###If leftPos == rightPos
                groupSubMatrix <- groupSubMatrix[pmatch(aa, names(groupSubMatrix))]}
        if (is.matrix(groupSubMatrix)){
            scoreOfGroup <- apply(groupSubMatrix, 1, mean)
            } else{
                scoreOfGroup <- mean(groupSubMatrix)
                }
        
        CNStatus <- ifelse(scoreOfGroup > upperCNThreshold, "DUP",
                   ifelse(scoreOfGroup < lowerCNThreshold, "DEL", "NORMAL"))
        
##Add into data.frame
        tempOut <- data.frame(aa, rep(breakSout[1], length(aa)),
                                     rep(breakSout[2], length(aa)),
                      scoreOfGroup, CNStatus)
        outData <- rbind(outData, tempOut)
        }
    
colnames(outData) <- c("Name", "Start", "End", "Score", "Status")

    ##################This option often results in errors
    if (useMixtureModel2ClusterGroup){
        scoreOfGroup <- outData$Score
        names(scoreOfGroup) <- outData$Name
        objectCluster <- new("clusteringCNVs", x = scoreOfGroup, k = 3)
        groupCNVofScore <- CNVrd2::groupCNVs(Object = objectCluster, autoDetermineGroup = TRUE)
        tempGroup <- groupCNVofScore$allGroups$Classification
        checkingGroup <- abs(sapply(split(scoreOfGroup, tempGroup), median))
        normalGroup <- as.numeric(names(which(checkingGroup == min(checkingGroup))))
        tempGroupOut <- ifelse(tempGroup == normalGroup, "NORMAL",
                       ifelse(tempGroup > normalGroup, "DUP", "DEL"))
        outData$Status <- tempGroupOut
              
    }
    
    return(outData)
}


####################################################################
#############################NEW FUNCTION###########################
####################################################################

detectBreakPointFromRD <- function(polymorphicObject,
                                   windows = 500, 
                                   genes, 
                                   quantileThreshold = 0.85,
                                   countThreshold = 5,
                                   sigMaTemp = windows/3,
                                   upperCNThreshold = 0.4,
                                   lowerCNThreshold = -0.4,
                                   detectAllRegion = FALSE,
                                   NTimes = 50,
             testType = c("Count", "positiveCount", "negativeCount", "SD"),
                                   useMixtureModel2ClusterGroup = FALSE, minLengthSV = 5000,
                                   epsilonCovDET = 0,
                                   NTimesThreshold = 20, NtransferToOtherPackage = 10,
                                   printOut = FALSE, singleSample = FALSE){

    testType <- match.arg(testType)

   
    listBreakPointOut <- NULL

    polymorphicRegion <- polymorphicObject
    
    subRegionMatrix <- polymorphicRegion$subRegionMatrix
    subRegion <- polymorphicRegion$subRegion

    if (testType == "Count")
        outSignal <- apply(subRegionMatrix, 2, function(x)
                           length(x[(x >= upperCNThreshold) | (x <= lowerCNThreshold)]))
    

    if (testType == "SD")
        outSignal <- apply(subRegionMatrix, 2, sd)
    if (testType == "positiveCount")
        outSignal <- apply(subRegionMatrix, 2, function(x) length(x[x >= upperCNThreshold]))
    if (testType == "negativeCount")
        outSignal <- apply(subRegionMatrix, 2, function(x) length(x[x <= lowerCNThreshold]))
###########################################################################
###########################################################################
    mSD <- data.frame(polymorphicRegion$subRegion, outSignal)
    mSD1 <- mSD[mSD[, 3] >= countThreshold, ]

    if (testType == "SD"){
        sdThreshold = quantile(mSD[, 3], quantileThreshold)
        mSD1 <- mSD[mSD[, 3] >= sdThreshold, ]
        }
###########################################################################
##Reduce to polymorphic regions############################################

    if (dim(mSD1)[1] > 0){
        
    mSD2 <- IRanges::reduce(IRanges(mSD1[, 1], mSD1[, 2]))

    geneMatrix <- matrix(genes, ncol = 2, byrow = TRUE)


    if (detectAllRegion){
        geneMatrix <- data.frame(start(mSD2), end(mSD2))

        if (printOut){
        message("gene Matrix: ")
        print(geneMatrix)
        message("=============")
    }
        
        geneMatrix <- geneMatrix[geneMatrix[, 2] - geneMatrix[, 1] >= minLengthSV,]
    }
    

    if (printOut)
        print(geneMatrix)

    
#######################################################################
####Scan for all genes################################################    
    for (kG in 1:dim(geneMatrix)[1]){


        gene <- as.numeric(geneMatrix[kG, ])
        tempGene <- IRanges::intersect(IRanges(gene[1], gene[2]), mSD2)

        if (printOut){
        message("Analysing the region: ")
        print(tempGene)
    }

        

        if (length(tempGene) > 0){

        
        tempGene <- tempGene[width(tempGene) == max(width(tempGene)),]
        tempGene <- mSD2[subjectHits(IRanges::findOverlaps(tempGene, mSD2)),]


###Find the idex of of the common region
        mSD3 <- mSD1[(mSD1[, 1] >= start(tempGene)) & (mSD1[, 2] <= end(tempGene)), ]
###Sub-region matrix
        subRegionMatrix <- polymorphicRegion$subRegionMatrix
###########################################################################    
###Got a submatrix only including the common region########################
        positionToPick <- as.numeric(rownames(mSD3))
        subSD2 <- subRegionMatrix[, positionToPick]

####################Function to classfify the data into different groups###
###########################################################################
        if (printOut)
            message("Running the clustering process")
        if (is.null(dim(subSD2))){

            if (printOut)
                message("Using Mclust to cluster for one-dimension data")

            if (sd(subSD2) < 0.01)
                classM <- rep(1, length(subSD2))
            else
                classM <-  mclust::Mclust(subSD2)$classification
             
            
        } else {

            if (singleSample)
                classM <- rep(1, length(subSD2[, 1]))
            else {

                if (max(dim(subSD2)) < NtransferToOtherPackage){
                    if (printOut)
                    message("Using Mclust to cluster for multi-dimension data")
                    
                    classM <-  mclust::Mclust(subSD2, modelNames="EII")$classification
                } else {
                    if (printOut)
                    message("Using HDclassif to cluster")
                    classM <- HDclassif::hddc(subSD2)$class
                }

            }
            
                
        }

        if (printOut){
        message("Number of Class M: ")
        print(table(classM))
    }

        names(classM) <- rownames(subRegionMatrix)
        classM1 <- cbind(classM, subSD2)
####################################################################
####################################################################
        tempPos <- unique(c(mSD3[, 1], mSD3[, 2]))
        tempPosForSub <- as.numeric(rownames(mSD3))
#########Add one left value and one right value
        tempPosForSub <- c(tempPosForSub[1] - 1, tempPosForSub, tempPosForSub[length(tempPosForSub)] + 1)

        if (tempPosForSub[1] == 0)
            tempPosForSub[1] <- 1
        if (tempPosForSub[length(tempPosForSub) - 1] == dim(subRegionMatrix)[2])
            tempPosForSub[length(tempPosForSub)] <- tempPosForSub[length(tempPosForSub) - 1]

###################################################################    
############subER: a matrix of common regions######################
        subER <- subRegionMatrix[, tempPosForSub]

##Change the first and last column to zero
        subER[, 1] <- rep(0, dim(subER)[1])
        subER[, dim(subER)[2]] <- rep(0, dim(subER)[1])
#############################################
###Finding breakpoints
        dataMatrix <- subER


            outData <- NULL


    for (ii in as.numeric(names(table(classM)))){
        #Combine classes and the matrix
        subER1 <- cbind(classM, dataMatrix[pmatch(names(classM), rownames(dataMatrix)),])

        #Retain only a sub-matrix in class being considered
        tempData <- subER1[subER1[, 1] ==ii, ]
        
        #tempData <- matrix(tempData, ncol = dim(subER1)[2])
        #tempData <- matrix(tempData, ncol = length(tempData))
        if (is.null(dim(tempData))){
            tempData <- tempData[-1]
            tempData <- matrix(tempData, nrow = 1)
            rownames(tempData) <- names(classM[classM == ii])
            } else
                tempData <- tempData[ , -1]

        scorePos <- NULL
        tempScore <- NULL
        ###Table random sample
        nSample <- dim(tempData)[1]
        tempTakeScore <- matrix(0, ncol = 2, nrow = nSample)

        if (is.null(NTimes))
            NTimes <- nSample
        if (NTimesThreshold > nSample)
            NTimesThreshold <- nSample

        if (nSample > 50)
            NTimesThreshold <- 5
        if (printOut){
            message(paste("\nRunning the resampling process: ", NTimesThreshold, " times\n", sep = ""))
            message(paste("Running the resampling process with nSample =  ", nSample, sep = ""))
        }


        for (k1 in 1:NTimesThreshold){

            if (printOut)
                message("Resample with k1 = ", k1)
            if (k1 == 1)
                tempData1 <- tempData
            else
                tempData1 <- tempData[sample(1:nSample, NTimes*nSample, replace=TRUE),]
              if (!is.matrix(tempData1))
                tempData1 <- matrix(tempData1, ncol = ncol(tempData))
            for (kk in 1:length(tempPos)){
                tempScore[kk] <- sum(abs(tempData1[, kk + 1] - tempData1[, kk]))
                }
            names(tempScore) <- tempPos
            for (kk in 1:length(tempPos)){
                scorePos[kk] <- sum(tempScore*dnorm(x = tempPos, mean = tempPos[kk], sd = sigMaTemp))

                }
            names(scorePos) <- tempPos
            tempOutScoreAA <- sort(scorePos, decreasing = TRUE)[1:2]
            tempOutScoreAA <- sort(as.numeric(names(tempOutScoreAA)))

                tempTakeScore[k1, ] <- tempOutScoreAA
            }

        
        if (is.matrix(tempTakeScore) | is.data.frame(tempTakeScore))
            tempTakeScore <- tempTakeScore[tempTakeScore[, 1] < tempTakeScore[, 2], ]

        
        
        #aa <- substr(rownames(tempData), 1, 7)
        aa <- rownames(tempData)

        breakSout <- tempTakeScore

     
#        message("breakSout: ", breakSout)
        if (!is.null(dim(tempTakeScore))){
        tempTakeScoreL <- sort(tempTakeScore[, 1])
        tempTakeScoreR <- sort(tempTakeScore[, 2])
        breakSout[1] <- tempTakeScoreL[floor((length(tempTakeScore) + 1)/2)]
        breakSout[2] <- tempTakeScoreR[floor((length(tempTakeScore) + 1)/2)]
    }
        ###Obtain positions of two breakpoints
        leftPos <- as.numeric(rownames(mSD3[mSD3[, 1] == breakSout[1],]))
        rightPos <- as.numeric(rownames(mSD3[mSD3[, 2] == breakSout[2],]))
        groupSubMatrix <- subRegionMatrix[, leftPos:rightPos]

####################Take segmentation results
        if (is.matrix(groupSubMatrix)){
            groupSubMatrix <- groupSubMatrix[pmatch(aa, rownames(groupSubMatrix)),]
            } else{
    ###If leftPos == rightPos
                groupSubMatrix <- groupSubMatrix[pmatch(aa, names(groupSubMatrix))]}
        if (is.matrix(groupSubMatrix)){
            scoreOfGroup <- apply(groupSubMatrix, 1, mean)
            } else{
                scoreOfGroup <- mean(groupSubMatrix)
                }
        CNStatus <- ifelse(scoreOfGroup >= upperCNThreshold, "DUP",
                   ifelse(scoreOfGroup <= lowerCNThreshold, "DEL", "NORMAL"))
        
##Add into data.frame
        tempOut <- data.frame(aa, rep(breakSout[1], length(aa)),
                                     rep(breakSout[2], length(aa)),
                      scoreOfGroup, CNStatus)
        outData <- rbind(outData, tempOut)

    }
    } else outData <- data.frame(rownames(subRegionMatrix),
                                 rep(gene[1], dim(subRegionMatrix)[1]),
                                 rep(gene[2], dim(subRegionMatrix)[1]),
                                 rep(0, dim(subRegionMatrix)[1]),
                                 rep("NORMAL", dim(subRegionMatrix)[1]))

        
    
colnames(outData) <- c("Name", "Start", "End", "Score", "Status")

        testStatus <- sort(names(table(outData$Status)))
        
       # if ((length(testStatus) > 1) | (testStatus[1] != "NORMAL"))
            listBreakPointOut[[kG]] <- outData

    }

    if (!is.null(listBreakPointOut))
    listBreakPointOut <- listBreakPointOut[!unlist(lapply(listBreakPointOut, is.null))]

} else listBreakPointOut <- NULL
    
    
    return(listBreakPointOut)
}
        
        
############################################################################
############################################################################
#######################NEW FUNCTION#########################################        
rdIdentifyBreakPointOfGroup <- function(dataMatrix, classM,
                                        upperCNThreshold = 0.4,
                                        lowerCNThreshold = -0.5,
                                        windows = 500,
                                        sigMaTemp = windows/3,
                                        NTimesThreshold = 20){
    outData <- NULL

    message("###Running rdIdentifyBreakPointOfGroup###")
    for (ii in 1:length(table(classM))){
        #Combine classes and the matrix
        subER1 <- cbind(classM, dataMatrix[pmatch(names(classM), rownames(dataMatrix)),])
        #Retain only a sub-matrix in class being considered
        tempData <- subER1[subER1[, 1] ==ii, ]
        #tempData <- matrix(tempData, ncol = dim(subER1)[2])
        #tempData <- matrix(tempData, ncol = length(tempData))
        if (is.null(dim(tempData))){
            tempData <- tempData[-1]
            tempData <- matrix(tempData, nrow = 1)
            rownames(tempData) <- names(classM[classM == ii])
            } else
                tempData <- tempData[ , -1]

        scorePos <- NULL
        tempScore <- NULL
        ###Table random sample
        nSample <- dim(tempData)[1]
        tempTakeScore <- matrix(0, ncol = 2, nrow = nSample)

        for (k1 in 1:nSample){
                    
            tempData1 <- tempData[sample(1:nSample, nSample, replace=TRUE),]
              if (!is.matrix(tempData1))
                tempData1 <- matrix(tempData1, ncol = ncol(tempData))
            for (kk in 1:length(tempPos)){
                tempScore[kk] <- sum(abs(tempData1[, kk + 1] - tempData1[, kk]))
                }
            names(tempScore) <- tempPos
            for (kk in 1:length(tempPos)){
                scorePos[kk] <- sum(tempScore*dnorm(x = tempPos, mean = tempPos[kk], sd = sigMaTemp))

                }
            names(scorePos) <- tempPos
            tempOutScoreAA <- sort(scorePos, decreasing = TRUE)[1:2]
            tempOutScoreAA <- sort(as.numeric(names(tempOutScoreAA)))
            tempTakeScore[k1, ] <- tempOutScoreAA
            }

        #aa <- substr(rownames(tempData), 1, 7)
        aa <- rownames(tempData)
        breakSout <- apply(tempTakeScore, 2, median)
        tempTakeScoreL <- sort(tempTakeScore[, 1])
        tempTakeScoreR <- sort(tempTakeScore[, 2])
        breakSout[1] <- tempTakeScoreL[floor((length(tempTakeScore) + 1)/2)]
        breakSout[2] <- tempTakeScoreR[floor((length(tempTakeScore) + 1)/2)]
        ###Obtain positions of two breakpoints
        leftPos <- as.numeric(rownames(mSD3[mSD3[, 1] == breakSout[1],]))
        rightPos <- as.numeric(rownames(mSD3[mSD3[, 2] == breakSout[2],]))
        groupSubMatrix <- subRegionMatrix[, leftPos:rightPos]

####################Take segmentation results
        if (is.matrix(groupSubMatrix)){
            groupSubMatrix <- groupSubMatrix[pmatch(aa, rownames(groupSubMatrix)),]
            } else{
    ###If leftPos == rightPos
                groupSubMatrix <- groupSubMatrix[pmatch(aa, names(groupSubMatrix))]}
        if (is.matrix(groupSubMatrix)){
            scoreOfGroup <- apply(groupSubMatrix, 1, mean)
            } else{
                scoreOfGroup <- mean(groupSubMatrix)
                }
        CNStatus <- ifelse(scoreOfGroup > upperCNThreshold, "DUP",
                   ifelse(scoreOfGroup < lowerCNThreshold, "DEL", "NORMAL"))
        
##Add into data.frame
        tempOut <- data.frame(aa, rep(breakSout[1], length(aa)),
                                     rep(breakSout[2], length(aa)),
                      scoreOfGroup, CNStatus)
        outData <- rbind(outData, tempOut)
        }
    
colnames(outData) <- c("Name", "Start", "End", "Score", "Status")
    return(outData)
}

correctMappability <- function(readCountMatrix, chr = NULL,
                               start = NULL, end = NULL, byMAPPABILITYcontent = 1,
                               mappabilityFile = NULL){

    ###If not use mappability file
    afterCorrectMappability <- readCountMatrix

    if (!is.null(mappabilityFile)){
        message("\nThere is a mappability file\n")

        mapData <- read.table(mappabilityFile)[, 3]
    mappabilityn <- mapData
    mappabilityList <- list()

    if (max(mappabilityn) <= 1)
        mappabilityn <- mappabilityn*100

        
     windowbyMAPPABILITYcontent <- seq(0, 100, by = byMAPPABILITYcontent)
        if (windowbyMAPPABILITYcontent[length(windowbyMAPPABILITYcontent)] <= 100)
            windowbyMAPPABILITYcontent[length(windowbyMAPPABILITYcontent)] <- 101
        
        dataframeToCorrect <- data.frame(windowbyMAPPABILITYcontent[-length(windowbyMAPPABILITYcontent)],
                                 windowbyMAPPABILITYcontent[-1])

    mappabilityList <- apply(dataframeToCorrect, 1, function(x){
        tempMap <- which((mappabilityn >= x[1]) & (mappabilityn < x[2]))
        return(tempMap)

        })
        

        mappabilityList <- unique(mappabilityList)

        mappabilityList <- mappabilityList[lapply(mappabilityList,length)>0]


        
    lengthMAPPABILITY <- length(mappabilityList)


#####################function################################
    correctMAPPABILITYforRow <- function(xRow){
        medianAll <- median(xRow)
        for (jj in 1:lengthMAPPABILITY){
            x = mappabilityList[[jj]]

            if (!is.null(x)){
                x1 = xRow[x]
                medianRegion <- median(x1)
                if (medianRegion != 0){
                    xRow[x] <- x1*medianAll/medianRegion
                    }
                else
                    xRow[x] <- xRow[x]
                }}
        return(xRow)
    }


    afterCorrectMappability <- t(apply(readCountMatrix, 1, correctMAPPABILITYforRow))
    }


    return(afterCorrectMappability)
      
}
