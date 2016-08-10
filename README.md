# SRBreak.

SRBreak is a read-depth and split-read package written in R for identifying copy-number variants 
in next-generation sequencing datasets.

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

## Go to [DATA](./DATA) to reproduce simulated data

## Docker

If you want to use a dockerized version, please download the image hoangtn/srbreak:

*sudo docker pull hoangtn/srbreak*

*sudo docker run -it hoangtn/srbreak /bin/bash*

*https://registry.hub.docker.com/u/hoangtn/srbreak/*

