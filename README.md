# Judith A. Sclafani - Portfolio

The follow collection contains samples of R scripts written for scientific projects. Examples have been selected to showcase my statistical and data processing capabilities.

## Table of Contents

## Building taxonomic abundance matrices
### **SKILLS:** dataframe and data file manipulation, data culling for quality and statistical robustness, function writing
This script is a part of a larger project to evaluate whether changes in abundance distributions through geologic time can be explained by neutral ecological theory. To do this, I wrote a program to transform a download from the Paleobiology Database into the format required by a previously developed maximum likehood estimation program written in PARI/GP. 

I needed to automate the analysis of thousands of fossil collections throughout geologic time. To do this, I developed a workflow, executed in bash, to cull the database download, split it by collection into multiple abundance matrices, save those matrices as .csv files, combine the abundance data with the additional information needed to run the PARI/GP program, excute the PARI/GP code, and combine the results into a new dataframe for additional analysis in R. 

Included here is the R script to perform the first part of this workflow, from culling the database download to writing the abundance matrices [PBDBtoAbundMatrices.r](/PBDBtoAbundMatrices.r)


## Evolutionary, ecological, and morphological similiarity 
### **SKILLS:** statistical analysis, ggplot


