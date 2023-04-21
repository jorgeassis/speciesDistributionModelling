## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
##
##                                      #####
##                                ####  #####
##                                ####       
##          ####                         
##         ##################             
##           ##################           
##       #######################
##   ##################################   
##  #######################################
##  ######################################
##  ###################################### 
##  ####################################
##  ##################################     
##  ####################                   
##  ###################                    
##  ##################                     
##  #################                      
##  ###############                                     
##      
##  theMarineDataScientist
##
##  github.com/jorgeassis
##  medium.com/themarinedatascientist
##  medium.com/@jorgemfa
##
## -------------------------------------------------------------------------------
##
##  SDM 3.0
##  R Pipelines for Marine Species Distribution Modelling
##
## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

setwd("/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Scripts")
source("NicheModelling/0_Config.R")
source("NicheModelling/Dependencies/Main.Functions.R")

# ----------------------------------------------

intertidal.mask.file <- "/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Data/Spatial/Coastline.tif"
bathymetry.file <- "/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Data/Spatial/BathymetryDepthMin.tif"

# ----------------------------------------------

study.extent <- extent(8,22,35,50)

Ocean.temperature.Max <- list.files(external.rasters.directory,pattern = "Ocean.temperature", full.names = TRUE)
Ocean.temperature.Max <- Ocean.temperature.Max[grep("Pred.Max",Ocean.temperature.Max)]
Ocean.temperature.Max <- stack(Ocean.temperature.Max)
Ocean.temperature.Max <- crop(Ocean.temperature.Max,study.extent)
Ocean.temperature.Max <- mask(Ocean.temperature.Max,shape)
names(Ocean.temperature.Max) <- 2000:2016

Dummy <- subset(Ocean.temperature.Max,1)
Dummy[1:ncell(Dummy)] <- rep(sample(1:100,15),ncell(Dummy))[1:ncell(Dummy)]

D <- subset(Dummy,1)
D[1:ncell(D)] <- rep(1.5,ncell(D))[1:ncell(D)]

# ----------------------------------------------

bathymetry <- raster(bathymetry.file)
bathymetry <- crop(bathymetry,study.extent)
bathymetry[bathymetry < -2] <- NA
bathymetry[bathymetry >= -2] <- 1

intertidal.mask <- raster(intertidal.mask.file)
intertidal.mask <- crop(intertidal.mask,study.extent)
intertidal.mask[is.na(intertidal.mask)] <- 0
shape <- sum(intertidal.mask,bathymetry)
shape[shape <= 0] <- NA
shape[!is.na(shape)] <- 1
shape <- mask(shape,shapefile(region.mask))
plot(shape)

# ----------------------------------------------

data <- read.xls("/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Data/Mechanistic/Experiments.xlsx")

presences <- read.occurrence.records(presence.records.file,exclude.na=FALSE,exclude.duplicated=FALSE,exclude.low.resolution=FALSE,exclude.out.of.range=FALSE,presence.records.mask.file=NULL,extra.columns=c("Year.of.record","Year.of.comparison","Occurrence.Data"))
presences[,1:2] <- relocate.coordinates.na(presences=presences[,1:2],rasters=shape,bathymetry=NULL,maximum.distance=500,use.species.depth=FALSE)
presences <- presences[which(!is.na(presences$Occurrence.Data)),]

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

models <- list()

for(m in 1:4) {
  
  if( m == 1) {
    project.name <- "Croacia_Ante"
    data.1 <- data[data$Location == "Coracia (Ante)",]
  }
  if( m == 2) {
    project.name <- "Croacia_Lilj"
    data.1 <- data[data$Location == "Cracia (Lilj)",]
  }
  if( m == 3) {
    project.name <- "Italy_Trieste"
    data.1 <- data[data$Location == "Italy (Trieste)",]
  }
  if( m == 4) {
    project.name <- "Italy_Venice"
    data.1 <- data[data$Location == "Italy (Venice)",]
  }
  
  # ----------------------------------------------------------------------------------------------------------
  
  Dummy.vect <- rep(0,length(data.1$Grwth.T2.))
  Dummy.vect[c(1,length(Dummy.vect))] <- 1
  data.to.model <- data.frame(Response=data.1$Grwth.T2.,
                              Temperature=as.numeric(as.character(data.1$Experiment)),
                              Dummy=Dummy.vect)
  
  data.to.model[data.to.model$Response <= 0,1] <- 0
  
  # ----------------------------------------------------------------------------------------------------------
  
  brt.learning.complex.span <- c(0.005,0.001,0.0005)
  brt.max.tree.depth <- 1:2
  
  cross.validation.results <- data.frame()
  parameters <- expand.grid(l.rate = brt.learning.complex.span,tree.c = brt.max.tree.depth)
  
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  parameters.results <- foreach( i = 1:nrow(parameters), .combine=rbind , .packages = c("dismo", "rJava","SDMTools","ENMeval")) %dopar% {
    
    l.rate <- parameters$l.rate[i]
    tree.c <- parameters$tree.c[i]
    
    model <- gbm.step( data=data.to.model, 
                       gbm.x = 2:3,
                       gbm.y = 1, 
                       family = "gaussian", 
                       plot.main = FALSE,
                       tree.complexity = tree.c, 
                       learning.rate = l.rate, 
                       bag.fraction = 0.5, 
                       n.folds=10,
                       step.size=50,
                       max.trees=5000,
                       silent=TRUE,
                       var.monotone = c(-1,-1),
                       verbose=FALSE)
    
    mydata <- data.frame(obs = data.to.model[,1] , fit = model$fitted)
    fit <- glm( obs ~ fit, data=mydata , family = gaussian)
    return(data.frame( tree.complexity = tree.c, learning.rate = l.rate, Dev= 1 - Dsquared(model = fit) , AIC = fit$aic ))
    
  }
  
  stopCluster(cl)
  rm(cl)
  
  # --------------------
  
  best.model <- parameters.results[which.max(parameters.results$Dev),]
  
  model.brt <- gbm.step( data=data.to.model, 
                         gbm.x = 2:3,
                         gbm.y = 1, 
                         family = "gaussian", 
                         plot.main = FALSE,
                         tree.complexity = best.model$tree.complexity, 
                         learning.rate = best.model$learning.rate, 
                         bag.fraction = 0.5, 
                         n.folds=10,
                         step.size=50,
                         max.trees=5000,
                         silent=TRUE,
                         var.monotone = c(-1,-1),
                         verbose=FALSE)
  
  models <- c(models,list(model.brt))
  
}

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Ensemble

y.span <- 2000:2016
ensemble <- list()

for( y in y.span ) {
  
  ensemble.i <- list()
  
  for(m in 1:4) {
    
    # rasters <- stack(subset(Ocean.temperature.Max,which(gsub("X","",names(Ocean.temperature.Max)) == y)),Dummy)
    # names(rasters) <- models[[m]]$var.names
    # rasters <- mask(rasters,shape)
    # 
    
    rasters <- subset(Ocean.temperature.Max,which(gsub("X","",names(Ocean.temperature.Max)) == y))
    
    if( y == 2016 ) { rasters <- rasters - 0.32 }
    
    rasters <- stack(rasters,D)
    names(rasters) <- model.brt$var.names
    rasters <- mask(rasters,shape)
    
    # 23.5
    # 
    # model.predicted <- subset(rasters,1)
    # model.predicted[model.predicted >= temp] <- 0
    # model.predicted[ model.predicted != 0] <- 1
    
    model.predicted <- predict( rasters , models[[m]] , n.trees=models[[m]]$n.trees,type="response")
    
    span.vals <- plot(models[[m]],1, return.grid=TRUE)$y
    span.vals.var <- numeric(length(span.vals)-1)
    for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- abs(span.vals[i+1] - span.vals[i]) }
    
    breaks <- which.max(span.vals.var) + 1
    model.predicted.t <- span.vals[breaks]
    
    model.predicted[model.predicted < model.predicted.t] <- 0
    model.predicted[model.predicted != 0] <- 1
    
    ensemble.i <- c(ensemble.i,list(model.predicted))
    
  }
  
  ensemble.i.t <- calc(stack(ensemble.i),median)
  ensemble <- c(ensemble,list(ensemble.i.t))
  
}

# -------------------------

model.predicted <- calc(stack(ensemble),mean)
plot(model.predicted)
model.predicted <- model.predicted - min(getValues(model.predicted),na.rm=T)
model.predicted <- model.predicted / max(getValues(model.predicted),na.rm=T)
plot(model.predicted)

writeRaster(model.predicted,filename=paste0(results.directory,"/Results/Mechanistic model/Suitability_",y.span[1],"_",y.span[length(y.span)],".tif"),format="GTiff",overwrite=T)

# ----------------------------------------
# 1950 | 2100

rasters <- raster("/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Data/Environment 1950/Temperature_Lt_max.tif")
rasters <- raster("/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Data/Environment 2100/Temperature_Lt_max_RCP45.tif")

rasters <- crop(rasters,study.extent)
rasters <- mask(rasters,shape)
rasters <- stack( rasters ,D )
names(rasters) <- models[[1]]$var.names

model.predicted <- calc(stack(sapply(1:4,function(x) { 
  
  span.vals <- plot(models[[x]],1, return.grid=TRUE)$y
  span.vals.var <- numeric(length(span.vals)-1)
  for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- abs(span.vals[i+1] - span.vals[i]) }
  
  breaks <- which.max(span.vals.var) + 0
  model.predicted.t <- span.vals[breaks]
  
  model.predicted <- predict( rasters , models[[x]] , n.trees=models[[x]]$n.trees,type="response") 
  model.predicted[model.predicted < model.predicted.t] <- 0
  model.predicted[model.predicted != 0] <- 1
  plot(model.predicted)
  return(model.predicted)
})),mean)

plot(model.predicted)
model.predicted <- model.predicted - min(getValues(model.predicted),na.rm=T)
model.predicted <- model.predicted / max(getValues(model.predicted),na.rm=T)
plot(model.predicted)

writeRaster(model.predicted,filename=paste0(results.directory,"/Results/Mechanistic model/Suitability_2090_2100.tif"),format="GTiff",overwrite=T)

# -------------------------
# No dispersal

final.ensemble.sum[final.ensemble.sum < length(y.span)] <- 0
final.ensemble.sum[final.ensemble.sum >= length(y.span)] <- 1
plot(final.ensemble.sum)

writeRaster(final.ensemble.sum,filename=paste0(results.directory,"/Results/Mechanistic model/Suitability_2016.tif"),format="GTiff",overwrite=T)

presences.i <- presences[presences$Year.of.record == 2016,]
points(presences.i[presences.i$Occurrence.Data == 1,],col="Black")
points(presences.i[presences.i$Occurrence.Data == 0,],col="Gray")

write.csv(presences.i,file=paste0(results.directory,"/Results/Mechanistic model/Records_2016.csv"))

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

y.span <- 2000:2016
dispersal.span <- 0:(length(y.span)-1)

# Temporal dispersal

final.accuracy <- data.frame()

for( disp in dispersal.span ) {
  
  final.accuracy.i <- data.frame()
  
  for( y in which(y.span %in% c(2016)) ) { #  1:length(y.span)
    
    if( disp >= y ) { next }
    
    index.i <- which(presences$Year.of.record == (y.span)[y])
    
    if( length(index.i) > 0) {
      
      observed <- presences$Occurrence.Data[index.i]
      
      reconstruction.before <- matrix(NA,ncol=length(y.span),nrow=length(index.i))
      
      for(r in 1:ncol(reconstruction.before)) {
        
        predicted.before <- extract( ensemble[[r]] , presences[index.i,1:2] )
        predicted.before[predicted.before != 1] <- 0
        reconstruction.before[,r] <- predicted.before
        
      }
      
      predicted <- numeric()
      
      for(l in 1:nrow(reconstruction.before)) {
        
        r <- rle(reconstruction.before[l,])
        r.l <- r$lengths[r$values == 0]
        
        if( TRUE %in% ( r.l > disp )  ) { predicted[l] <- 0 }
        if( ( ! FALSE %in% ( r.l <= disp ) ) & reconstruction.before[l,ncol(reconstruction.before)] == 1  ) { predicted[l] <- 1 }
        if( reconstruction.before[l,ncol(reconstruction.before)] == 0  ) { predicted[l] <- 0 }
        
        
      } 
      
      accuracy.i <- data.frame(observed=observed,predicted=predicted)
      final.accuracy.i <- rbind(final.accuracy.i,accuracy.i[complete.cases(accuracy.i),])
      
    }
    
  }
  
  if( nrow(accuracy.i) > 0) {
    
    accuracy.i <- accuracy( final.accuracy.i$observed , final.accuracy.i$predicted , threshold =100 )
    accuracy.i <- accuracy.i[accuracy.i$threshold != 0 ,]
    accuracy.i <- accuracy.i[ which.max( accuracy.i$sensitivity + accuracy.i$specificity ),]
    final.accuracy <- rbind(final.accuracy , data.frame( dispersal=disp , accuracy=accuracy.i) )
    
  }
}

# ----------------

final.accuracy
tss.final <- max(final.accuracy$accuracy.sensitivity+final.accuracy$accuracy.specificity-1)
tss.final

# ----------------

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(final.accuracy$dispersal,final.accuracy$accuracy.sensitivity+final.accuracy$accuracy.specificity-1,lty=2,col="#878787",type="l",ylab="",xlab="Temporal dispersal (year)",axes=FALSE)
axis(2,las=2,col="white",col.ticks="Black")
axis(1,las=0,col="white",col.ticks="Black")
axis(1,las=0,col="white",col.ticks="Black",at=2,label="2")
box()
title(ylab="Variation in accuracy (TSS)",mgp=c(4,1,0))
abline(v = 2, lty=1, col="Gray")

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

y.span <- 2000:2016
reconstruction.area.unlimited <- numeric(length(y.span))
reconstruction.area <- numeric(length(y.span))

for(r in 1:ncol(reconstruction.before)) {
  area.predicted <- ensemble[[r]]
  area.predicted[area.predicted != 1] <- 0
  area.predicted <- raster::area(area.predicted) * area.predicted
  reconstruction.area.unlimited[r] <- sum(getValues(area.predicted),na.rm=TRUE) # km2
}

for(r in 1:ncol(reconstruction.before)) {
  area.predicted <- calc( subset(stack(ensemble) , 1:r) , min)
  area.predicted[area.predicted != 1] <- 0
  area.predicted <- raster::area(area.predicted) * area.predicted
  reconstruction.area[r] <- sum(getValues(area.predicted),na.rm=TRUE) # km2
}

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(y.span,reconstruction.area.unlimited/1000,col="black",type="l",ylab="",xlab="Time (year)",axes=FALSE,ylim=c(0,max(reconstruction.area.unlimited)/1000),lty=3)
axis(2,las=2,col="white",col.ticks="Black")
axis(1,las=0,col="white",col.ticks="Black")
box()
title(ylab="Suitable habitats (10^3 km2)",mgp=c(4,1,0))
lines(y.span,reconstruction.area/1000,lty=2)
#Suitability_Time : 12-12

# --------------------------------------------------------------------------------------------------------