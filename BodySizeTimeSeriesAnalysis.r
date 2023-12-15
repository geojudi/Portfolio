
#########################################################################################
####################   FIT EVOLUTIONARY MODELS TO BODY SIZE DATA   ######################
#########################################################################################
# This code calculates the change in mean and variance of body size over time for an order and fits the change in mean to several common evolutionary models.

# Input: a .csv containing body size data for genera within orders
InputFile <- 'AllBrachsCleaned.csv'

# Output: a list containing, for each order, Akaike weights comparing strength of support for each model
#########################################################################################



##########################################################################################
##########################     LOAD REQUIRED PACKAGES     ################################
##########################################################################################
# contains functions to fit time series of evolutionary data to common evolutionary models
library(paleoTS)

# adds the geologic time scale to plots
library(deeptime)

# dependency for deeptime
library(zoo)

# for plotting
library(ggplot2)
##########################################################################################



##############  Read in body size (BS) data    #############################

AllBS <- read.csv(InputFile)
############################################################################



#################  Perform data culling procedures  #########################

# split data by orders
OrderSplit <- split(AllBS,AllBS$order)

# remove any rows without an order name
OrderSplit <- OrderSplit[-which(names(OrderSplit)=="")]

# remove orders with fewer than 20 genera
no_in_ord <- unlist(lapply(OrderSplit,nrow))
OrderSplit <- OrderSplit[-which(names(OrderSplit) %in% names(which(no_in_ord<20)))]

# return to a dataframe format for plotting
culledOrdBS <-do.call(rbind.data.frame,OrderSplit)
############################################################################



####################  Plot body size data  ################################
# plot body size data through time for each order
ggplot()+
	# make grey lines to mark the 5 mass extinction events
	geom_vline(xintercept = 443.8,color='grey')+
 	geom_vline(xintercept = 372.2,color='grey')+
 	geom_vline(xintercept = 251.9,color='grey')+
  	geom_vline(xintercept = 201.3,color='grey')+
  	geom_vline(xintercept = 66,color='grey')+
  	
  	# plot body size with a different panel for each order
  	geom_segment(data=culledOrdBS,aes(x=fad_age,xend=lad_age,y=log(calc_max_vol),yend=log(calc_max_vol), color=order))+ 
  	facet_wrap(facets=vars(order),ncol=4)+
  	labs(y="log(volume)", x="mya")+
  	
  	# put the geologic time scale on the x axis
  	coord_geo(xlim=c(541,0),ylim=c(-5,15),height=unit(.5,'line'),size=1.7)+
 	scale_x_reverse()+
 	
 	# make the plots look cleaner
  	theme(panel.background=element_rect(fill='white',color='black'))+
  	theme(panel.grid=element_blank())+
  	theme(strip.background=element_blank())+
  	theme(panel.spacing=unit(.2,'lines'))+
  	theme(legend.position = "none")
############################################################################



##########  Convert body size data to an object of class paleoTS  ##########

# set up a list to contain data. Each element equates to an order in the dataset 
PaleoTSformat <- vector("list",length(OrderSplit))
names(PaleoTSformat) <- names(OrderSplit)

# create a vector to represent 10 my time bins 
time <-seq(541,1,-10)

# make time vector end at 0 my (present day)
time <- c(time,0)

# convert body size for each order by each 10 my time bin
for(i in 1:length(OrderSplit)){
  # extract a single order
  taxon <- OrderSplit[[i]]
  taxon.name <- names(OrderSplit)[i]
  
  # determine the order's time span by matching first and last occurrence dates to 10 my time bins
  ages <- sort(unique(c(unique(taxon$fad_age), unique(taxon$lad_age))), decreasing=T)

  start <- which(time<=max(ages))
  end <- which(time>=min(ages))
  
  t2 <- time[intersect(start, end)]
  
  # for each time bin, calculate biovolume mean and variance and the number of genera present
  meanVol <- c()
  varVol <- c()
  n <- c()
  
  for(j in 1:length(t2)){
  	sizes <- c()
  	for(k in 1:nrow(taxon)){
  	
  		# if a taxon is present in time bin j, extract its volume
  		
  		if(taxon$fad_age[k]>=t2[j]&taxon$lad_age[k]<=t2[j]){
  			sizes[k] <- log(taxon$calc_max_vol[k])
  		} else {
  			k <- k+1
  		}
  	}
  		
  	# remove NAs	
  	sizes <- sizes[!is.na(sizes)]
  	
	# calculate biovolume mean, variance, and number for genera present in time bin j
    if(length(sizes)==0){
    	meanVol[j] <- NA
    	varVol[j] <- NA
    	n[j] <- NA
    } else {
    	meanVol[j] <- mean(sizes)
   	 	varVol[j] <- var(sizes)
    	n[j] <- length(sizes)
    }
  }	
  
  # give NAs a value close to zero so that the model fitting function can run
  varVol[is.na(varVol)] <- 0.00001
  meanVol[is.na(meanVol)] <- 0.00005
  
  PaleoTSformat[[i]] <- as.paleoTS(mm=meanVol,vv=varVol,nn=n,tt=t2)
}
############################################################################



############  Use paleoTS to fit evolutionary models to data  #############
###### set up lists to hold model results
# AIC results from the model comparison
modelresults <- vector("list",length(PaleoTSformat))

# fit models
for(i in 1:length(PaleoTSformat)){

 	# unbiased random walk
  	urw <- fitSimple(PaleoTSformat[[i]],model="URW",pool=F)
  
  	# generalized random walk
  	grw <- fitSimple(PaleoTSformat[[i]],model="GRW",pool=F)
  
  	# stasis 
  	stasis <- fitSimple(PaleoTSformat[[i]],model="Stasis",pool=F)

	# Ornstein-Uhlenbeck
  	ou <- fitSimple(PaleoTSformat[[i]],model="OU",pool=F)
  
  	# shift in evolutionary mode
  	st.urw <- fitModeShift(PaleoTSformat[[i]],minb=3,pool=F,order='Stasis-RW',rw.model='URW')
  	st.grw <- fitModeShift(PaleoTSformat[[i]],minb=3,pool=F,order='Stasis-RW',rw.model='GRW')
  	urw.st <- fitModeShift(PaleoTSformat[[i]],minb=3,pool=F,order='RW-Stasis',rw.model='URW')
  	grw.st <- fitModeShift(PaleoTSformat[[i]],minb=3,pool=F,order='RW-Stasis',rw.model='GRW')
  
	# calculate Akaike weights for model results  
	modelresults[[i]] <- compareModels(urw,grw,stasis,st.stasis,ou,st.urw,st.grw,urw.st,grw.st)
}
names(modelresults) <- names(PaleoTSformat)
