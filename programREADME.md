This folder contains all of the files needed to run superautomaticPARI.sh, a program that estimates Hubble's theta for abundance data downloaded from the PaleobioDB. 

last updated: Sept 9, 2021
created: Aug 29, 2014
Judi Sclafani 

Contents of this readme: 
I. description of files in this folder 
II. protocol for running superautomaticPARI.sh



#########################################
I. Files in this folder:
1.  PaleobioDBtoCSV.r - an R script that culls a large PaleobioDB download file to remove taxa that are not-suspension feeders and collections that contain non-numeric abundance data, transforms the culled data into an abundance matrix, and then, splits the abundance matrix by grouping collections by reference number, 10 million year time bin, and depositional environment. The outputs of this are script are: 
	1. culledPBDB.csv - the culled version of the input matrix
	2. Abundances.csv - an abundance matrix containing all downloaded data
	3. CollectionInfo.csv - a matrix with associated categorical data for each collection
	4. AbundMatrices - a folder containing all of the split abundance matrices
		file names for these matrices are formatted: "ReferenceNumber_10myAgeBin_environment"
		NA at the end of some file names is a meaningless artifact of correcting for spacing and punctuation in depositional environment data
		
	UPDATE (Sept 2021): this code was created under the old PBDB download column naming system. Column names are now slightly different and the code has been updated to reflect this.
	***NOTE BEFORE RUNNING ANALYSIS*** This is set up to run all taxonomic classes and orders designated as suspension feeders. To analyze a different trophic level, the function RemoveNonSuspensionFeeders that is in the R script will need to be adjusted. 
	
As written, this program will run any PBDB download file with abundance data, time bins, and environment fields. For the program to run correctly, the name of the file should be  pbdbdata-occs.csv and the column names should be unaltered from those output by the Paleobiology Database download screen. Data should be downloaded without comments (or comment columns should be removed in Excel because punctuation in the comments can offset the column alignment in the .csv file. 
	
	
2. superautomaticPARI.sh - the bash code to execute the R script, create and setup folders  and input files in the format needed for the Etienne (2007) theta estimation program, and run the theta estimation in PARI/gp


3. TemplateFilesAlphaDiv - folder containing files needed to run the program automaticPARI.sh, which creates and sets up the folders and input files needed to run Etienne's theta estimation program for a single abundance matrix. This program is executed by superautomaticPARI.sh so that it can be iterated over a multiple abundance matrices. Files in this folder are: 
	1. automaticPARI.sh - the bash script for the program
	2. OriginalFiles - a folder containing Etienne's PARI/gp scripts
	3. PARIabundance.r - R script to transform an abundance matrix into PARI format
	
	
	
##########################################
Steps for running superautomaticPARI.sh: 
1. copy all of the files from TemplateFilesSuperAuto into a new folder that will hold all of the inputs and outputs of the program
2. place within the folder a PaleobioDB download file, making sure the file name is pbdbdata-occs.csv 
3. in a terminal window, set the working directory to the folder you created in step 1
4. in the command line type: sh superautomaticPARI.sh 
5. Patiently await results. It will take many hours for a large file to run.

*** If you are not using a PBDB data download and have your own abundance matrices, you can ignore the parts of the code that manipulate the PBDB download file. Instead the procedure should be as follows: 

1. copy all of the files from TemplateFilesSuperAuto into a new folder that will hold all of the inputs and outputs of the program
2. Make a folder within the folder created in step 1 and name it AbundMatrices
3. within 'AbundMatrices' make a folder for each collection for which you want to calculate a theta value. Place the abundance matrix in the folder.
4. open superautomaticPARI.sh and delete lines of code for handling PBDB data, noted by the BEGIN and END comments
5. in a terminal window, set the working directory to the folder you created in step 1
6. in the command line type: sh superautomaticPARI.sh 
7. Patiently await results. It will take many hours for a large file to run.
	
