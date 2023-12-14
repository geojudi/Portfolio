# Convert PALEOBIOLOGY DATABASE DOWNLOAD TO ABUNDANCE MATRICES

# This script will generate abundance matrices for all collections in a Paleobiology Databse (PBDB) download file. The download file must contain the abundance, time bins, and environment fields as collections are split by collection number, time, and environment. For the program to run correctly, the name of the file should be  pbdbdata-occs.csv and the column names should be unaltered from those output by the PBDB. Data should be downloaded without comments or comment columns should be removed in Excel because punctuation in the comments can offset the column alignment in the .csv file. 

# input: a .csv file downloaded from the PBDB. This file should be named pbdbdata-occs.csv

# output: multiple .csv files, each of which contains an abundance matrix. The number of output files is determined by number of collections in downloaded file and the number of time bins and depositional environments in each unique collection

#### NOTE BEFORE RUNNING ANALYSIS: This is set up to run all taxonomic classes and orders designated as suspension feeders. To analyze a different trophic level, the function RemoveNonSuspensionFeeders will need to be adjusted. 

###############################################################################################################################
# last updated: Sept 9, 2021

# UPDATE (Sept 2021): this code was created under the old PBDB download column naming system. Column names are now slightly different, and the code has been updated to reflect this.

# created: Aug 29, 2014 by Judith Sclafani
###############################################################################################################################


###############################################################################################################################
################################################### NECESSARY FUNCTIONS #######################################################
###############################################################################################################################

# Remove taxa that are not classified as suspension feeders from the download

RemoveNonSuspensionFeeders <- function(pbdbData) {
	ClassesAllowed <- c("Scaphopoda", "Anthozoa", "Echinoidea", "Bivalvia", "Hyolitha", "Hyolithomorpha", "Orthothecimorpha", "Rostroconchia", "Lingulata", "Paterinata", "Chileata", "Kutorginata", "Obolellata", "Rhynchonellata", "Strophomenata", "Archaeocyatha", "Irregulares", "Regulares", "Calcarea", "Demospongea", "Heteractinida", "Hexactinellida", "Stromatoperoidea", "Gymnolaemata", "Stenolaemata", "Blastoidea")
    
	OrdersAllowed <- c("Thoracica", "Vermetidae", "Calyptraeidae", "Euomphalina", "Murchisoniina", "Coenothecalia", "Gorgonacea", "Helioporacea", "Cystiphyllida", "Heterocorallia","Stauriida", "Auloporida", "Favositida", "Halysitida", "Heliolitida", "Lichenariida", "Sarcinulida","Tetradiida", "Actiniaria")

	ClassesPresent <- as.vector(pbdbData$class)
	OrdersPresent <- as.vector(pbdbData$order)

	ClassIndices <- which(ClassesPresent %in% ClassesAllowed)
	OrderIndices <- which(OrdersPresent %in% OrdersAllowed)

	Indices <- sort(union(ClassIndices, OrderIndices))

	AbundanceMatrix <- pbdbData[Indices, ]
	
	AbundanceMatrix
}

# Generate a matrix of abundance values from a PBDB download. Formats matrix with genera as columns and collections as rows

GenerateAbundanceMatrix <- function(pbdbData){
	UniqueCollections <- as.vector(unique(pbdbData$collection_no))
	UniqueGenera <- as.vector(unique(pbdbData$genus))
	Abundances <- matrix(data=0, nrow=length(UniqueCollections), ncol=length(UniqueGenera), dimnames=list(UniqueCollections, UniqueGenera))
	
	for (i in 1:nrow(pbdbData)){ 
		CollectionNumber <- which(UniqueCollections==pbdbData$collection_no[i])
		GenusNumber <- which(UniqueGenera==pbdbData$genus[i]) 
		
		AbundanceValue <- as.numeric(as.vector(pbdbData$abund_value[i]))
		Abundances[CollectionNumber,GenusNumber] <- Abundances[CollectionNumber,GenusNumber] + AbundanceValue
	}
		Abundances
}
#################################################################################################################################


## CREATE MATRICES OF ABUNDANCE DATA AND COLLECTION INFO ##

# set data file to read in data
InputFile <- 'pbdbdata-occs.csv'
pbdbData <- read.csv(file=InputFile, as.is=TRUE)

# remove any data that does not have an associated time bin or environment
pbdbData <- pbdbData[which(as.vector(pbdbData$time_bin)!=""), ]
pbdbData <- pbdbData[which(as.vector(pbdbData$environment)!=""), ]

# remove taxa that are not in the list of suspension feeders
pbdbData <- RemoveNonSuspensionFeeders(pbdbData)


## remove collections with problematic abundance values

# first identify values that contain non-numeric characters
unallowableValues <- unique(pbdbData$abund_value)[grep('[^0-9]', unique(pbdbData$abund_value))]

# then determine which collection contains those values
collectionstoRemove <- pbdbData$collection_no[which(pbdbData$abund_value %in% unallowableValues)]

# remove an entire collection if it contains any non-numeric abundance values
indices <- which(pbdbData$collection_no %in% collectionstoRemove)
pbdbData <- pbdbData[-indices, ]


## transform data download to abundance matrix
AbundanceMatrix <- GenerateAbundanceMatrix(pbdbData)
collectionSize <- rowSums(AbundanceMatrix)

# cull abundance matrix to include only collections of sizes larger than 5
# this determination is based on the summary statistics and to improve statistical robustness
indices <- which(collectionSize > 5)
AbundanceMatrix <- AbundanceMatrix[indices, ] 

# cull abundance matrix to remove collections that have fewer than 2 genera
numGenera <- c()
for(i in 1:nrow(AbundanceMatrix)){
    numGenera[i] <- length(which(AbundanceMatrix[i, ]!=0))
}

indices <- which(numGenera >1)
AbundanceMatrix <- AbundanceMatrix[indices, ]

# remove any genera that, after the culling, have an abundance of zero 
AbundanceMatrix <- AbundanceMatrix[ ,which(colSums(AbundanceMatrix)!=0)]

## make collection info data table

# make a list of unique collections
collections <- unique(pbdbData$collection_no)

# set and fill desired information fields
age <- c()
paleolat <- c()
paleolong <- c()
environment <- c()
country <- c()
state <- c()
plate <- c()
lithification <- c()
lithology <- c()
referenceno <- c()

for(i in 1:length(collections)){ 
	age[i] <- as.vector(unique(pbdbData$Xtime_bin[which(pbdbData$collection_no==collections[i])]))
	paleolat[i] <- as.vector(unique(pbdbData$paleolat[which(pbdbData$collection_no==collections[i])]))
	paleolong[i] <- as.vector(unique(pbdbData$paleolng[which(pbdbData$collection_no==collections[i])]))
	environment[i] <- as.vector(unique(pbdbData$environment[which(pbdbData$collection_no==collections[i])]))
	country[i] <- as.vector(unique(pbdbData$cc[which(pbdbData$collection_no==collections[i])]))
	state[i] <- as.vector(unique(pbdbData$state[which(pbdbData$collection_no==collections[i])]))
	plate[i] <- as.vector(unique(pbdbData$geoplate[which(pbdbData$collection_no==collections[i])]))
	lithification[i] <- as.vector(unique(pbdbData$lithification[which(pbdbData$collection_no==collections[i])]))
	lithology[i] <- as.vector(unique(pbdbData$lithology1[which(pbdbData$collection_no==collections[i])]))
	referenceno[i] <- as.vector(unique(pbdbData$reference_no[which(pbdbData$collection_no==collections[i])]))
}

# combine fields generated above, collection number, and collection size into data frame
CollectionInfo <- cbind(collections, age, paleolat, paleolong, environment, country, state, plate, lithification, lithology, referenceno, collectionSize)
colnames(CollectionInfo) <- c("collectionNo", "age", "paleolat", "paleolong", "environment", "country", "state", "plate", "lithification", "lithology", "ReferenceNo", "collectionSize")

# transform the CollectionInfo matrix to a dataframe so that dollar notation can be used later
CollectionInfo <- as.data.frame(CollectionInfo)

# this is to make sure the collection info and abundance matrices contain the same collections
CollectionInfo <- CollectionInfo[which(CollectionInfo$collectionNo %in% rownames(AbundanceMatrix)), ]
CollectionInfo <- droplevels(CollectionInfo)

## write culled paleobiodb data, abundance matrix, collection info matrix

write.csv(pbdbData, file='culledPBDB.csv')
write.csv(AbundanceMatrix, file="Abundances.csv")
write.csv(CollectionInfo, file='CollectionInfo.csv')


### MAKE .CSV FILES OF ABUNDANCE DATA
## these files serve as the data input for running automaticPARI.sh in bash
## Data is organized as collections and are split by reference number first
## reference numbers with collections from more than one age or depositional environment are then split accordingly

# set working directory to a previously created folder to contain the .csv files
setwd('AbundMatrices')

# split the collection data by reference number
splitRefNo <- split(CollectionInfo, CollectionInfo$ReferenceNo)

## this is to identify the names of environments that have characters that aren't allowed in file names
# generally these characters are punctuation
# these problems are corrected within the for loop
probEnvir <- as.vector(unique(CollectionInfo$environment)[grep('[[:punct:]]', unique(CollectionInfo$environment))])


## this loop makes the .csv files

# loop through each reference number 
for (i in 1:length(splitRefNo)){
    RefNo <- names(splitRefNo)[i]
    collecIndex <- which(CollectionInfo$ReferenceNo==RefNo)
		
    # make matrix, include desired collections, remove genera that have an abundance of zero
    matrix <- AbundanceMatrix[collecIndex,which(colSums(AbundanceMatrix[collecIndex, , drop=FALSE])!=0), drop=FALSE]

    # determine number of time bins in reference number
    age <- as.vector(unique(splitRefNo[[i]]$age))

    # for each time bin within a reference number
    for(j in 1:length(age)){

        # identify collections within reference number that are from the time bin
        collecAge <- splitRefNo[[i]]$collectionNo[which(splitRefNo[[i]]$age==age[j])]
        ageIndex <- which(rownames(matrix) %in% collecAge)
        ageData <- matrix[ageIndex, which(colSums(matrix[ageIndex, , drop=FALSE])!=0), drop=FALSE]

        # determine number of depositional environments samples within age bin are from
        ageEnvir <- as.vector(unique(splitRefNo[[i]]$environment[ageIndex]))
			
        # for each depositional environment within a time bin
        for(k in 1:length(ageEnvir)){
					
            # identify collections within reference number and time bin that are from the depositional environment
            envirIndex <- which(splitRefNo[[i]]$environment[ageIndex]==ageEnvir[k])
            envirData <- ageData[envirIndex, which(colSums(ageData[envirIndex, ,drop=FALSE])!=0), drop=FALSE]
					
            # for environment names that contain punctuation
            if(ageEnvir[k] %in% probEnvir){
						
                # this corrects for punctuation in environment names by using only a portion of the full name
                correctedEnvir <- strsplit(ageEnvir[k], '[[:punct:]]')[[1]][1]
							
                # remove spaces from the environment and age names to make file name
                envirFileName <- paste(strsplit(correctedEnvir, ' ')[[1]][1], strsplit(correctedEnvir, ' ')[[1]][2], sep='')
                ageFileName <- paste(strsplit(age[j], ' ')[[1]][1], strsplit(age[j], ' ')[[1]][2], sep='')
                FileName <- paste(RefNo, ageFileName, envirFileName, ".csv", sep='_')
                write.csv(envirData, file=FileName)
						
            # for environment names that do no contain punctuation
            } else {
                # remove spaces from the environment and age names to make file name
                envirFileName <- paste(strsplit(ageEnvir[k], ' ')[[1]][1], strsplit(ageEnvir[k], ' ')[[1]][2], sep='')
                ageFileName <- paste(strsplit(age[j], ' ')[[1]][1], strsplit(age[j], ' ')[[1]][2], sep='')
                FileName <- paste(RefNo, ageFileName, envirFileName, ".csv", sep='_')
                write.csv(envirData, file=FileName)
            }
        }
    }
}	

