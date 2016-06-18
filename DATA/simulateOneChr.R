library("BSgenome.Hsapiens.UCSC.hg19")
###Extract the chromosome 21
chr <- "chr21"
tempChr <- unmasked(Hsapiens[[chr]])

##Read structural-variant (SV) file
inputList <- read.table("ListSVofChr21ForsimulateOneChrFile.txt", header = FALSE)

##Extract start and end positions of SV regions
startP <- inputList[, 1]
endP <- inputList[, 2]

##Convert to start and end segments to extract from Chr 21
endPosition <- c(startP, length(tempChr))
startPosition <- c(1, as.numeric(endP))
outChr <- NULL
for (ii in 1:length(startPosition)){

    outChr <- c(outChr, tempChr[startPosition[ii]:endPosition[ii]])

}

finalChr <- do.call(c, outChr)
ff <- file(paste("Chr21.simulateDeletionsFrom1000Genomes.fa", sep = ""), "w")
write('>chr21_Del', ff)
writeLines(paste(finalChr), ff)
close(ff)

