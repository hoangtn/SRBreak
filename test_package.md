## This is a faster way to check results of SRBreak

```{}
git clone https://github.com/hoangtn/CNVrd2.git 

git clone https://github.com/hoangtn/SRBreak.git

git clone https://github.com/hoangtn/DataFile.git


cd SRBreak/
cat SRBreak/R/SRBreak.R > ../allSourceFile.R
cd ../CNVrd2/
cat CNVrd2/R/*R >> ../allSourceFile.R 
cd ..

```

Then, go inside R and install three packages: VariantAnnotation, BSgenome.Hsapiens.UCSC.hg19, mclust

```{}
##Install packages

source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

install.packages("mclust")

##Load SRBreak source

source("allSourceFile.R")
windows = 500
dirBamFile = "DataFile/BamSimulatedData/"
st = 101100001
en = 102100000
chr = "chr1"
library("VariantAnnotation") 
library("BSgenome.Hsapiens.UCSC.hg19")
library("mclust")

system.time(outputSRBreak <- SRBreak(readDepthWindow = windows,##read-depth window size
                                     
                                     chr = chr, ##Chromosome name
                                                                                                            
                                     st = st, ##Start position
                                     
                                     en = en, ##End position
                                     
                                     dirBamFile = dirBamFile, ##Bam files' directory
                                     
                                     detectAllRegion = TRUE, ##Set this = TRUE in order to obtain all CNV regions
                                     
                                     rdQualityMapping = 0, ##Mapping quality (used in read-depth method)
                                                             
                                     testType = "Count", ##Test type
                                     
                                     correctGC = TRUE, ##Correct GC content
                                     
                                     upperCNThreshold = 0.25, ##Larger than this threshold is duplication
                                     
                                     lowerCNThreshold = -0.25, ##Smaller than this threshold is deletion
                                     
                                     countThreshold = 1, ##Number of duplications/deletions 
                                     
                                     minLengthSV = 1000, ##Minimum length of a duplication/deletion event

                                    #inputRawReadCountMatrix = readCountMatrix,
                                     usingPairedEnds = FALSE, ##Not use paired-end information,
                                                            
                                     thresholdOfIntersectionBetweenRDandPEM = 0.9 #Not useful if only single-end reads used
                                     
                                    ))

rawoutputSV <- outputSRBreak$svResult

###Write out
write.table(rawoutputSV, "SVresults.txt", quote = FALSE)


##Next steps are aimed to compare simulated breakpoints and predicted breakpoints
##We can see all simulated breakpoints by extracting sample names

sampleNames <- dir(dirBamFile, "bam$")

sampleOfSV <- sampleNames[grep("Dup|Del", sampleNames)]
sampleOfSV

listSV <- unique(t(sapply(strsplit(sampleOfSV, ".", fixed = TRUE), function(x) x[4:5])))
    
    
simulatedSV  <- IRanges(as.integer(listSV[, 1]), as.integer(listSV[, 2]))
    
sort(simulatedSV)

outputSV <- IRanges(as.integer(rawoutputSV[, 2]), as.integer(rawoutputSV[, 3]))

###Extract predicted breakpoints around simulated breakpoints
###Check whether the simulated breakpoints are in these results 
outputFinal <- outputSV[(start(outputSV) > (101545220 - 5000)) & (end(outputSV) < (101630000 + 5000)) ]

outputFinal

```
