# SRBreak.

SRBreak is a read-depth and split-read package written in R for identifying copy-number variants
in next-generation sequencing datasets.

### Note: SBReak was designed to work for multiple samples. It can work for >= 2 samples, but we suggest that users should use >= 5 samples as in the work tested in our paper.

If users only have one sample from their real data, one possible solution is that users can download some samples from the 1000 Genomes project (known SVs) as reference samples.

The package needs chromosome, start and end positions to run. This information should be in the BAM files of users. 

SRBreak will try to adjust the start position (to remove 'NNN' regions) if the option *adjustStartPosition* is set to TRUE. If just analysing for one specific region, this option should be set to FALSE.  

### SRBreak only analyses one chromosome. Therefore, if users have multiple chromosomes then he/she should analyse them separately. 

### Segmental duplication information: SRBreak needs this information to refine final results.

Segmental duplication information can be obtained from here: http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
and should be formatted the same as the file [GRCh37GenomicSuperDup.tab.onlyChr21](./GRCh37GenomicSuperDup.tab.onlyChr21) for each chromosome.


Please follow this example [test_package.md](./test_package.md) to learn more about parameters used inside the package.


## Installation

Users can run the package by using one of these two ways:

### 1. From source codes: [test_package.md](./test_package.md)


### 2. Install as follows

SRBreak relies on the CNVrd2, therefore please install the package CNVrd2 first:

Note: firstly, you need to install jags (http://mcmc-jags.sourceforge.net/), samtools (http://samtools.sourceforge.net/) and then:

*git clone https://github.com/hoangtn/CNVrd2.git*

*cd CNVrd2; R CMD INSTALL CNVrd2*

*git clone https://github.com/hoangtn/SRBreak.git*

*cd SRBreak; R CMD INSTALL SRBreak*




## Test Data

Download the test data from: https://github.com/hoangtn/DataFile

### Test directly from source code:  [test_package.md](./test_package.md)

### If you already installed the package, an example can be seen here: [TestSRBreak.ipynb](./TestSRBreak.ipynb)

## Go to [DATA](./DATA) to reproduce simulated data

## Docker

If you want to use a dockerized version, please download the image hoangtn/srbreak:

*sudo docker pull hoangtn/srbreak*

*sudo docker run -it hoangtn/srbreak /bin/bash*

*https://registry.hub.docker.com/u/hoangtn/srbreak/*

