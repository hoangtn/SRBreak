## This is a faster way to install/check results of SRBreak

Users need to install *samtools* if required.


```{}
git clone https://github.com/hoangtn/CNVrd2.git

git clone https://github.com/hoangtn/SRBreak.git


cd SRBreak/
cat SRBreak/R/SRBreak.R > ../allSourceFile.R
cd ../CNVrd2/
cat CNVrd2/R/*R >> ../allSourceFile.R
#sed -i .bak "s/CNVRd2:://g" ../allSourceFile.R
#rm ../allSourceFile.R.bak
cd ..

```

Then, go inside R and install three packages: Rsamtools, VariantAnnotation, BSgenome.Hsapiens.UCSC.hg19, mclust

```{}
##Install packages
source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
biocLite("VariantAnnotation")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("DNAcopy")
install.packages("mclust")
```
## Test for whole chromosome chr21

Users should download whole-chromosome files from:

https://www.dropbox.com/sh/wejg4r37kdkjokc/AADglfiSzny1bJ4w1ROAVERoa?dl=0

After that, unzip the files into one directory (e.g., Chr21BamFile)


```{}
unzip Chr21BamFile.zip -d Chr21BamFile
```
### In R

```{}

source("allSourceFile.R")

library("VariantAnnotation")
library("BSgenome.Hsapiens.UCSC.hg19")
library("mclust")
library("Rsamtools")

#link file: https://www.dropbox.com/sh/wejg4r37kdkjokc/AADglfiSzny1bJ4w1ROAVERoa?dl=0
#
ePos <- 48129895 ##End position
windows = 500

dirBamFile = "./Chr21BamFile/"

segmentalDuplicationFile = paste0(dirBamFile, "GRCh37GenomicSuperDup.tab.onlyChr21")
##Segmental duplication information can be obtained from here: http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
##and should be formatted as the file GRCh37GenomicSuperDup.tab.onlyChr21 for each chromosome

st = 9500001
en = ePos
st = 9500001 ##Start position; if we don't know, we can set st = 1 and SRBreak can adjust this value inside to remove unknown regions
chr = "chr21" ##Note: these bam files are only from chr21

system.time(outputSRBreak <- SRBreak(readDepthWindow = windows,##read-depth window size

                                     chr = chr, ##Chromosome name

                                     st = st, ##Start position

                                     en = en, ##End position

                                     dirBamFile = dirBamFile, ##Bam files' directory

                                     detectAllRegion = TRUE, ##Set this = TRUE in order to obtain all CNV regions

                                     rdQualityMapping = 0, ##Mapping quality (used in read-depth method)

                                     testType = "Count", ##Test type

                                     correctGC =  TRUE, #   FALSE, ##Correct GC content

                                     upperCNThreshold = 0.25, ##Larger than this threshold is duplication

                                     lowerCNThreshold = -0.25, ##Smaller than this threshold is deletion

                                     countThreshold = 2, ##Number of duplications/deletions: should be >= 2

                                     minLengthSV = 1000, ##Minimum length of a duplication/deletion event

#                                     mappabilityFile = "wgEncodeCrgMapabilityAlign100mer.bigWig.Window.500.start.1.end.48129895.txt",
##Mappability information is not a compulsory option, it should be turned off
                                     usingPairedEnds = FALSE, ##Not use paired-end information,

                                     thresholdOfIntersectionBetweenRDandPEM = 0.9 #Not useful if only single-end reads used
                                     , segmentalDuplicationFile = segmentalDuplicationFile,
                                     adjustStartPosition = TRUE ##Adjust positions
                                    ))

rawoutputSV <- outputSRBreak$svResult

dim(rawoutputSV)


finalResult <- rawoutputSV
dim(finalResult)

```

This file [./test_chr21_SRBreak.ipynb](test_chr21_SRBreak.ipynb) shows results for the test data. The test results include: whole chr21, a specific region, and a specific region for two samples.
