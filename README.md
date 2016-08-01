# SRBreak.

SRBreak is a read-depth and split-read package written in R for identifying copy-number variants 
in next-generation sequencing datasets.

## Installation

SRBreak relies on the CNVrd2, therefore please install the package CNVrd2 first:

*git clone https://github.com/hoangtn/CNVrd2.git*

*R CMD INSTALL CNVrd2*  

*git clone https://github.com/hoangtn/SRBreak.git*

*R CMD INSTALL SRBreak*

## Test Run

[TestSRBreak.ipynb](./TestSRBreak.ipynb) explains how to run the package.

## Test Data

Download the test data from: https://github.com/hoangtn/DataFile

## Go to [DATA](./DATA) to reproduce simulated data

## Docker

If you want to use a dockerized version, please download the image hoangtn/srbreak:

*sudo docker pull hoangtn/srbreak*

*sudo docker run -it hoangtn/srbreak /bin/bash*

*https://registry.hub.docker.com/u/hoangtn/srbreak/*

