# Judith A. Sclafani - Portfolio

The follow collection contains samples of R scripts written for scientific projects. Examples have been selected to showcase my statistical and data processing capabilities.

## Table of Contents

## Building taxonomic abundance matrices
### **SKILLS:** dataframe and data file manipulation, data culling for quality and statistical robustness, function writing
This script is a part of a larger project to evaluate whether changes in abundance distributions through geologic time can be explained by neutral ecological theory. I needed to automate the analysis of thousands of fossil collections downloaded from the Paleobiology Database using a maximum likelihood estimation program written in PARI/GP. To do this, I developed a workflow, executed in bash, to first, use R to cull the database download, split it by collection into multiple abundance matrices, and save those matrices as .csv files. Then, combine the abundance data with the additional information needed to run the PARI/GP program, excute the PARI/GP code, and combine the results into a new dataframe for additional analysis in R. 

Included here is the R script to perform the first part of this workflow, from culling the database download to writing the abundance matrices [PBDBtoAbundMatrices.r](/PBDBtoAbundMatrices.r)


## Evolutionary, ecological, and morphological similiarity 
### **SKILLS:** statistical analysis, ggplot
This script contains code that I have regularly used in my research to quantify species evolution, extinction, and environmmental preference. Specifically, the goal was to evaluate whether the most morphologically and ecologically similiar taxa are the most evolutionarily related. To do this, I calculated the pairwise branch length distance between tips on a phylogenetic tree. I also calculated the similarity (as pairwise distances) between metrics that quantify morphology and ecology. This script was written to work with the brachiopod data I collected during my dissertation and then further developed during my postdoctoral research. 

Included here is an R script where I compiled pairwise distance calculations and the procedure to plot results 
