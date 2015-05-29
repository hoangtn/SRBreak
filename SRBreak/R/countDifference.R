##This function is used to count differences between populations
countDifference <- function(data, popName, regions = NULL,
                         upperThreshold = 0.2,
                         lowerThreshold = -0.2){

    ##data is a data.frame/matrix of segmentation results
    ##popName is a data.frame including two columns: sample names and their corresponding populations

    
    popName <- popName[pmatch(rownames(data), popName[, 1]), ]

    ###Calculate number of populations
    nPops <- length(table(popName[, 2]))
 

    Vst <- apply(data, 2, function(x){

        x <- ifelse((x > lowerThreshold) & (x < upperThreshold), 0, 1)
                 
          bTemp <- data.frame(data = x,
                              pop = as.character(popName[, 2])) ##Combine a column and population information
          
          nameValueVST <- c()
          vResults <- c()
          indexVST = 1 ##Set index for the first pair

          for (iST in 1:(nPops-1)){ #iST for the first pop
              for (jST in (iST+1):nPops){ #jST for the second pop
                  bTempIJ <- bTemp[(bTemp[, 2] == popNames[iST]) |
                             (bTemp[, 2] == popNames[jST]), ] #Extract the two pops
                  

                  bVSTtemp <- sapply(split(bTempIJ, bTempIJ$pop),
                                     function(x) dim(x[x[, 1] > 0, ])[1])

                  vResults[indexVST] <- abs(bVSTtemp[1] - bVSTtemp[2])

                  nameValueVST[indexVST] <- paste(popNames[iST], "-", popNames[jST], sep = "")
                  indexVST <- indexVST + 1
              }}
          names(vResults) <- nameValueVST

          return(vResults)
          })

    if (!is.null(regions)){
        Vst <- data.frame(regions, Vst)
        colnames(Vst) <- c("Start", "End", "Vst")
    }
    
    return(Vst = Vst)
}
