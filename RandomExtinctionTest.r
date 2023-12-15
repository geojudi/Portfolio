# to figure out if taxa extinction/recovery at each stage is random in the morphospace
this code requires the following matrices made in 'R Disparity CG script.r'
#     1. MS = contains PCO coordinates used to make the morphospace
#     2. PA = contains presence/absence data of each taxa in each stage
# these objects can be loaded into the R workspace using: 



MS_Kat <- MS[which(PA[,which(colnames(PA)=='Katian')]==1),1:2]
MS_Hirn <- MS[which(PA[,which(colnames(PA)=='Hirnantian')]==1),1:2]
MS_Rhudd <- MS[which(PA[,which(colnames(PA)=='Rhuddanian')]==1),1:2]



# to combine katian and hirnantian data. all hirnantian taxa are in the katian
MS_Hirn <- MS_Kat



Hirn_Centroid <- c(mean(MS_Hirn[,1]), mean(MS_Hirn[,2]))
Rhudd_Centroid <- c(mean(MS_Rhudd[,1]), mean(MS_Rhudd[,2]))


Survivors <- rownames(MS_Hirn)[which(rownames(MS_Hirn) %in% rownames(MS_Rhudd))]
Extinct <- rownames(MS_Hirn)[-which(rownames(MS_Hirn) %in% rownames(MS_Rhudd))]
Originators <- rownames(MS_Rhudd)[-which(rownames(MS_Rhudd) %in% rownames(MS_Hirn))]

NumExtinct <- length(Extinct)
NumOrig <- length(Originators)
CentroidExt <- c(mean(MS_Hirn[which(rownames(MS_Hirn) %in% Extinct), 1]), mean(MS_Hirn[which(rownames(MS_Hirn) %in% Extinct), 2]))
CentroidOrig <- c(mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originators), 1]), mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originators), 2]))
CentroidSurv <- c(mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Survivors), 1]), mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Survivors), 2]))

# bootstrap to see if extinction is random

BootstrapExtinctCentroid <- matrix(data=0, nrow=1000, ncol=2)

for(i in 1:1000){
	RandomExtinctionIndex <- sample(seq(1,nrow(MS_Hirn)), replace=T, NumExtinct)
	ExtinctRand <- rownames(MS_Hirn)[RandomExtinctionIndex]
	CentroidExtinct <- c(mean(MS_Hirn[which(rownames(MS_Hirn) %in% ExtinctRand), 1]), mean(MS_Hirn[which(rownames(MS_Hirn) %in% ExtinctRand), 2]))
	BootstrapExtinctCentroid[i,1] <- CentroidExtinct[1]
	BootstrapExtinctCentroid[i,2] <- CentroidExtinct[2]
}


# bootstrap to see if origination is random

BootstrapOrigCentroid <- matrix(data=0, nrow=1000, ncol=2)

for(i in 1:1000){
	RandomOriginationIndex <- sample(seq(1,nrow(MS_Rhudd)), replace=T, NumOrig)
	Originate <- rownames(MS_Rhudd)[RandomOriginationIndex]
	CentroidOriginate <- c(mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originate), 1]), mean(MS_Rhudd[which(rownames(MS_Rhudd) %in% Originate), 2]))
	BootstrapOrigCentroid[i,1] <- CentroidOriginate[1]
	BootstrapOrigCentroid[i,2] <- CentroidOriginate[2]
}
