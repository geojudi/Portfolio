# Judith A. Sclafani - Portfolio

The following collection contains samples of R scripts written for scientific projects. Examples have been selected to showcase my statistical and data processing capabilities.

## Building taxonomic abundance matrices
### **SKILLS:** dataframe and data file manipulation, function writing, commenting for publication
This script is a part of a larger project to evaluate whether changes in abundance distributions through geologic time can be explained by neutral ecological theory. I needed to automate the analysis of thousands of fossil collections downloaded from the Paleobiology Database using a maximum likelihood estimation program written in PARI/GP. To do this, I developed a workflow, executed in bash, to first, use R to cull the database download, split it by collection into multiple abundance matrices, and save those matrices as .csv files. Then, combine the abundance data with the additional information needed to run the PARI/GP program, excute the PARI/GP code, and combine the results into a new dataframe for additional analysis in R. 

Included here is the R script to perform the first part of this workflow, from culling the database download to writing the abundance matrices [PBDBtoAbundMatrices.r](/PBDBtoAbundMatrices.r)

## Comparing extinction and origination data to random expectation
### **SKILLS:** ordination, resampling
This script is a part of my dissertation research in which I quantified morphological change at a mass extinction event. My goal with this part of the project was to ordinate morphological character data to determine how variability within a phylogenetic tree changed across the extinction. I then bootstrapped the centroid of ordination space (morphospace) to determine 1) if extinct taxa in pre-extinction morphospace were randomly distributed (they were) and 2) if newly originated taxa in post-extinction morphospace were randomly distributed (they weren't). 

Included here is the script that I used to perform these analyses [RandomExtinctionTest.r](RandomExtinctionTest.r)

## Body size evolution
### **SKILLS:** statistical model fitting, ggplot
This script contains code that I recently developed to evaluate brachiopod body size evolution over geologic time. Specifically, the goal was to determine, for each order, which evolutionary model best explains change in size. To do this, I transformed the raw data into an object of a class required by the model-fitting package. Then, I automated fitting and comparing multiple models for all orders. I have also included an example of code to plot the body size data using ggplot.

Included here is an R script that contains this procedure [BodySizeTimeSeriesAnalysis.r](BodySizeTimeSeriesAnalysis.r)
