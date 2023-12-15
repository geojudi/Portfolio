#############################################################################################
##########     Compare extinction and origination to random expectations       ##############
#############################################################################################
# This code ordinates categorical morphological traits for all taxa in a phylogenetic analysis, while keeping track of their duration in geologic time. It then uses a bootstrap to determine whether morphological change at a mass extinction event differs from a random expectation. 

# Input: two data tables
#### 1. a morphology character matrix - containing trait designations for taxa in a phylogenetic tree
MorphologyFile <- 'https://raw.githubusercontent.com/geojudi/Portfolio/master/data/Strophid_M_STAGE.txt'
#### 2. a presence absence matrix where rows are taxa, columns are geologic stages, and cells are filled with 1 (present) or 0 (absent)
StagesFile <- 'https://raw.githubusercontent.com/geojudi/Portfolio/master/data/Strophid_D_STAGE.txt'

# Output: results of the ordination and bootstrap analyses
#############################################################################################


#############################################################################################
##########################      Required Packages    ########################################
#############################################################################################
# for PCO calculation
library(ecodist)

#############  read in the morphology and stage data and format it for analysis  ###############
m <- read.table(MorphologyFile)

# make as a matrix AND REMOVES FIRST COLLUMN IE THE TAXA NAMES
M <- as.matrix(m[,-c(1)])

# make rownames vector and adds back into matrix
rownames(M)<-(m[,1]) 			

# read in the stratigraphic range data (stage level)
pa <- read.table(StagesFile,header=T)

# make as a matrix AND REMOVES FIRST COLLUMN IE THE TAXA NAMES
PA <- as.matrix(pa[,-c(1)])

# make rownames vector and adds back into matrix
rownames(PA) <- (pa[,1])
################################################################################################


###################    Ordinate morphological character data   ################################
# make a Bray-Curtis distance matrix	
DistM <- bcdist(M)

# ordinate using Principal Coordinates Analysis (PCO)
ResultsPCO <- pco(DistM)	

# separate out coordinates to a new object
MS <- ResultsPCO$vectors
################################################################################################


#############    Bootstrap morphospace values across a mass extinction   #######################
### separate out data from the geologic stages of interest
# before the extinction. all Hirnantian taxa are present in the Katian
MS_Hirn <- MS[which(PA[,which(colnames(PA)=='Katian')]==1),1:2]

# after the extinction
MS_Rhudd <- MS[which(PA[,which(colnames(PA)=='Rhuddanian')]==1),1:2]

# calculate the centroid in ordination space for all taxa present both before and after the extinction
Hirn_Centroid <- c(mean(MS_Hirn[,1]), mean(MS_Hirn[,2]))
Rhudd_Centroid <- c(mean(MS_Rhudd[,1]), mean(MS_Rhudd[,2]))

# determine which taxa survived, went extinct at, and originated after the event
Survivors <- rownames(MS_Hirn)[which(rownames(MS_Hirn) %in% rownames(MS_Rhudd))]
Extinct <- rownames(MS_Hirn)[-which(rownames(MS_Hirn) %in% rownames(MS_Rhudd))]
Originators <- rownames(MS_Rhudd)[-which(rownames(MS_Rhudd) %in% rownames(MS_Hirn))]

# calculate the centroid of each group
CentroidExt <- c(mean(MS_Hirn[which(rownames(MS_Hirn) %in% Extinct), 1]), mean(MS_Hirn[which(rownames(MS_Hirn) %in% Extinct), 2]))
CentroidOrig <- c(mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originators), 1]), mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originators), 2]))
CentroidSurv <- c(mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Survivors), 1]), mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Survivors), 2]))

# determine the number of taxa that go extinct and originate
NumExtinct <- length(Extinct)
NumOrig <- length(Originators)

# use the number of taxa that go extinct to resample the pre-extinction data with replacement and calculate a centroid
BootstrapExtinctCentroid <- matrix(data=0, nrow=1000, ncol=2)
for(i in 1:1000){
	RandomExtinctionIndex <- sample(seq(1,nrow(MS_Hirn)), replace=T, NumExtinct)
	ExtinctRand <- rownames(MS_Hirn)[RandomExtinctionIndex]
	CentroidExtinct <- c(mean(MS_Hirn[which(rownames(MS_Hirn) %in% ExtinctRand), 1]), mean(MS_Hirn[which(rownames(MS_Hirn) %in% ExtinctRand), 2]))
	BootstrapExtinctCentroid[i,1] <- CentroidExtinct[1]
	BootstrapExtinctCentroid[i,2] <- CentroidExtinct[2]
}

# use the number of taxa that originate to resample the post-extinction data with replacement and calculate a centroid
BootstrapOrigCentroid <- matrix(data=0, nrow=1000, ncol=2)
for(i in 1:1000){
	RandomOriginationIndex <- sample(seq(1,nrow(MS_Rhudd)), replace=T, NumOrig)
	Originate <- rownames(MS_Rhudd)[RandomOriginationIndex]
	CentroidOriginate <- c(mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originate), 1]), mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originate), 2]))
	BootstrapOrigCentroid[i,1] <- CentroidOriginate[1]
	BootstrapOrigCentroid[i,2] <- CentroidOriginate[2]
}
