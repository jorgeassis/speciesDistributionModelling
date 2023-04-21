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

setwd("/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Scripts/NicheModelling")
#setwd("~/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Scripts/NicheModelling")
source("0_Config.R")
source("Dependencies/Main.Functions.R")

# ------------------

# options( java.parameters = "-Xmx4g" )

# ------------------------------------------------------------------------------------

study.extent <- extent(8,22,35,50)

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

Nitrate.Min <- list.files(external.rasters.directory,pattern = "Nitrate", full.names = TRUE)
Nitrate.Min <- stack(Nitrate.Min)
Nitrate.Min <- crop(Nitrate.Min,study.extent)
Nitrate.Min <- mask(Nitrate.Min,shape)
names(Nitrate.Min) <- 2000:2016

Ocean.temperature.Min <- list.files(external.rasters.directory,pattern = "Ocean.temperature", full.names = TRUE)
Ocean.temperature.Min <- Ocean.temperature.Min[grep("Pred.Min",Ocean.temperature.Min)]
Ocean.temperature.Min <- stack(Ocean.temperature.Min)
Ocean.temperature.Min <- crop(Ocean.temperature.Min,study.extent)
Ocean.temperature.Min <- mask(Ocean.temperature.Min,shape)
names(Ocean.temperature.Min) <- 2000:2016

Ocean.temperature.Max <- list.files(external.rasters.directory,pattern = "Ocean.temperature", full.names = TRUE)
Ocean.temperature.Max <- Ocean.temperature.Max[grep("Pred.Max",Ocean.temperature.Max)]
Ocean.temperature.Max <- stack(Ocean.temperature.Max)
Ocean.temperature.Max <- crop(Ocean.temperature.Max,study.extent)
Ocean.temperature.Max <- mask(Ocean.temperature.Max,shape)
names(Ocean.temperature.Max) <- 2000:2016

Air.temperature.Min <- list.files(external.rasters.directory,pattern = "Air.temperature", full.names = TRUE)
Air.temperature.Min <- Air.temperature.Min[grep("Pred.Min",Air.temperature.Min)]
Air.temperature.Min <- stack(Air.temperature.Min)
Air.temperature.Min <- crop(Air.temperature.Min,study.extent)
Air.temperature.Min <- mask(Air.temperature.Min,shape)
names(Air.temperature.Min) <- 2000:2016

Air.temperature.Max <- list.files(external.rasters.directory,pattern = "Air.temperature", full.names = TRUE)
Air.temperature.Max <- Air.temperature.Max[grep("Pred.Max",Air.temperature.Max)]
Air.temperature.Max <- stack(Air.temperature.Max)
Air.temperature.Max <- crop(Air.temperature.Max,study.extent)
Air.temperature.Max <- mask(Air.temperature.Max,shape)
names(Air.temperature.Max) <- 2000:2016

Dummy <- Ocean.temperature.Max
Dummy[1:ncell(Dummy)] <- rep(1:2,ncell(Dummy))[1:ncell(Dummy)]

D <- subset(Dummy,1)
D[1:ncell(D)] <- rep(1.5,ncell(D))[1:ncell(D)]

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

project.name <- "Temperature"

presences <- read.occurrence.records(presence.records.file,exclude.na=TRUE,exclude.duplicated=FALSE,exclude.low.resolution=FALSE,exclude.out.of.range=TRUE,presence.records.mask.file=NULL,extra.columns=c("Year.of.record","Year.of.comparison","Occurrence.Data"))
absences <- presences[ which( ! is.na(presences$Year.of.record) & presences$Occurrence.Data == 0 ) , - which(colnames(presences) == "Occurrence.Data") ]
presences <- presences[  which( ! is.na(presences$Year.of.record) & presences$Occurrence.Data == 1 ) , - which(colnames(presences) == "Occurrence.Data") ]

# --------------------

absences <- absences[absences$Year.of.record >= 2000,]
presences <- presences[presences$Year.of.record >= 2000,]
absences[,1:2] <- relocate.coordinates.na(presences=absences[,1:2],rasters=shape,bathymetry=NULL,maximum.distance=500,use.species.depth=FALSE)
presences[,1:2] <- relocate.coordinates.na(presences=presences[,1:2],rasters=shape,bathymetry=NULL,maximum.distance=500,use.species.depth=FALSE)
nrow(absences)
nrow(presences)

# ------------------------------------------------------------------------------------

presences <- cbind(presences,data.frame(Nitrate.Min = NA ,
                                        Ocean.temperature.Min = NA ,
                                        Ocean.temperature.Max = NA ,
                                        Air.temperature.Min = NA ,
                                        Air.temperature.Max = NA,
                                        Dummy=NA))

absences <- cbind(absences,data.frame(Nitrate.Min = NA ,
                                      Ocean.temperature.Min = NA ,
                                      Ocean.temperature.Max = NA ,
                                      Air.temperature.Min = NA ,
                                      Air.temperature.Max = NA,
                                      Dummy=NA))

variable.monotonic.response <- data.frame(t(c(1,1,-1,1,-1,+1)))
colnames(variable.monotonic.response) <- colnames(presences)[-c(1:4)] ; variable.monotonic.response

# ------------------------------------------------------------------------------------

vars <- colnames(variable.monotonic.response)
vars

# Subset data

vars <- vars[c(3,6)]
vars

presences <- presences[presences$Year.of.record %in% c(2015,2016),c(1:4,which(colnames(presences) %in% vars))]
absences <- absences[absences$Year.of.record %in% c(2015,2016),c(1:4,which(colnames(absences) %in% vars))]

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

for( i in 1:nrow(presences) ) {
  
  Year.of.record <- presences$Year.of.record[i]
  Year.of.comparison <- presences$Year.of.comparison[i]
  
  if(is.na(Year.of.comparison)) { Year.of.comparison <- Year.of.record }
  
  coords <- presences[ i, 1:2]
  
  for( v in 1:length(vars)) {
    
    raster.var <- subset(get(vars[v]) ,  grep(Year.of.record,names(get(vars[v]))) ) # grep(Year.of.comparison,names(get(vars[v]))):
    
    raster.var <- subset(get(vars[v]) ,  grep(Year.of.record-1,names(get(vars[v]))):grep(Year.of.record,names(get(vars[v]))) ) # grep(Year.of.comparison,names(get(vars[v]))):
    
    if( variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ][v] > 0 ) { raster.var <- calc(raster.var,max) }
    if( variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ][v] < 0 ) { raster.var <- calc(raster.var,min) }
    
    presences[i, which(colnames(presences) == vars[v]) ] <- as.numeric(extract(raster.var,coords))
    
  }
}

# ------------------------

for( i in 1:nrow(absences) ) {
  
  Year.of.record <- absences$Year.of.record[i]
  Year.of.comparison <- absences$Year.of.comparison[i]
  coords <- absences[ i, 1:2]
  
  for( v in 1:length(vars)) {
    
    raster.var <- subset(get(vars[v]) ,  grep(Year.of.record,names(get(vars[v]))) ) # grep(Year.of.comparison,names(get(vars[v]))):
    
    raster.var <- subset(get(vars[v]) ,  grep(Year.of.record-1,names(get(vars[v]))):grep(Year.of.record,names(get(vars[v]))) ) # grep(Year.of.comparison,names(get(vars[v]))):
    
    if( variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ][v] > 0 ) { raster.var <- calc(raster.var,min) }
    if( variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ][v] < 0 ) { raster.var <- calc(raster.var,max) }
    
    # raster.var <- subset(get(vars[v]) , grep(Year.of.comparison,names(get(vars[v]))):grep(Year.of.record,names(get(vars[v]))) ) #
    # 
    # if( variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ][v] > 0 ) { raster.var <- calc(raster.var,min) }
    # if( variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ][v] < 0 ) { raster.var <- calc(raster.var,max) }
    
    absences[i, which(colnames(absences) == vars[v]) ] <- as.numeric(extract(raster.var,coords))
    
  }
}

absences <- absences[absences$Year.of.comparison != absences$Year.of.record , ]

# ------------------------

data.to.model <- rbind( data.frame(Presence=rep(1,nrow(presences)),presences) ,
                        data.frame(Presence=rep(0,nrow(absences)),absences) )

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

predictors <- colnames(variable.monotonic.response)
cross.validation.results <- data.frame()

parameters <- expand.grid(l.rate = brt.learning.complex.span,tree.c = 1:length(vars))

cl <- makeCluster(n.cores, "SOCK")
registerDoSNOW(cl)

parameters.results <- foreach( i = 1:nrow(parameters), .combine=rbind , .packages = c("dismo", "rJava","SDMTools","ENMeval")) %dopar% {
  
  l.rate <- parameters$l.rate[i]
  tree.c <- parameters$tree.c[i]
  
  model <- NULL
  while(is.null(model)) {
    
    model <- gbm.step( data=data.to.model, 
                       gbm.x = which(colnames(data.to.model) %in% vars),
                       gbm.y = 1, 
                       family = "bernoulli", 
                       plot.main = FALSE,
                       tree.complexity = tree.c, 
                       learning.rate = l.rate, 
                       bag.fraction = 0.5, 
                       n.folds=10,
                       step.size=10,
                       max.trees=1000,
                       silent=TRUE,
                       var.monotone = variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ] ,
                       verbose=FALSE)
  }
  
  mydata <- data.frame(obs = data.to.model[,1] , fit = model$fitted)
  
  fit <- glm( obs ~ fit, data=mydata , family = binomial)
  return(data.frame( tree.complexity = tree.c, learning.rate = l.rate, Dev= 1 - model$cv.statistics$deviance.mean , AIC=fit$aic))
  
}

stopCluster(cl)
rm(cl)

# ------------------------

best.model <- parameters.results[which.max(parameters.results$Dev),]
model.brt <- NULL

while(is.null(model.brt)) {
  model.brt <- gbm.step( data=data.to.model, 
                         gbm.x = which(colnames(data.to.model) %in% vars),
                         gbm.y = 1, 
                         family = "bernoulli", 
                         plot.main = FALSE,
                         tree.complexity = best.model$tree.complexity, 
                         learning.rate = best.model$learning.rate, 
                         bag.fraction = 0.5, 
                         n.folds=10,
                         step.size=10,
                         max.trees=1000,
                         silent=TRUE,
                         var.monotone = variable.monotonic.response[ which(colnames(variable.monotonic.response) %in% vars) ] ,
                         verbose=FALSE)
}

# ------------------------

plot(model.brt,1,type="response")
plot(model.brt,2,type="response")

# ------------------------

plot(data.to.model[,1],model.brt$fitted)
1 - model.brt$cv.statistics$deviance.mean

# ------------------------

sum.brt <- summary.model(model.brt,print.data=FALSE)
colnames(sum.brt) <- c("Variable","Contrib.BRT")
sum.brt

# -----------------------

cor(data.to.model$Air.temperature.Max,data.to.model$Ocean.temperature.Max)

# 10 10 
#Temp_PartialPlot_TempMax

model.brt$var.names[1]

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(plot(model.brt,1, return.grid=TRUE,type = "response")[,1],plot(model.brt,1, return.grid=TRUE,type = "response")[,2],lty=1,col="#878787",type="l",ylab="",xlab="Temperature (ºC)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Effect on response",mgp=c(4,1,0)) 
abline(v = 24.89, lty=3, col="#878787")
abline(v = 25.15, lty=3, col="#878787")

model.brt$var.names[2]

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(plot(model.brt,2, return.grid=TRUE,type = "response")[,1],plot(model.brt,2, return.grid=TRUE,type = "response")[,2],lty=1,col="#878787",type="l",ylab="",xlab="Temperature (ºC)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Effect on response",mgp=c(4,1,0)) 
abline(v = 9.65, lty=3, col="#878787")
abline(v = 12.3, lty=3, col="#878787")

model.brt$var.names[1]

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(plot(model.brt,1, return.grid=TRUE,type = "response")[,1],plot(model.brt,1, return.grid=TRUE,type = "response")[,2],lty=1,col="#878787",type="l",ylab="",xlab="Temperature (ºC)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Effect on response",mgp=c(4,1,0)) 
abline(v = 24.4, lty=3, col="#878787")
abline(v = 25.7, lty=3, col="#878787")

# ----------------------------------------

names(rasters) <- model.brt$var.names
model.predicted <- predict( rasters , model.brt , n.trees=model.brt$n.trees,type="response")
plot(model.predicted)

writeRaster(model.predicted,filename=paste0(results.directory,"/Results/Correlative model/Suitability_2000_2016.tif"),format="GTiff",overwrite=T)

# ----------------------------------------
# 1950 | 2100

rasters <- raster("/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Data/Environment 1950/Temperature_Lt_max.tif")
rasters <- raster("/Volumes/Laminaria/Dropbox/Manuscripts/Fucus virsoides and Climate Change/Data/Environment 2100/Temperature_Lt_max_RCP45.tif")

rasters <- crop(rasters,study.extent)
rasters <- mask(rasters,shape)
rasters <- stack( rasters ,D )
names(rasters) <- model.brt$var.names

span.vals <- plot(model.brt,1, return.grid=TRUE,type="response")$y
span.vals.var <- numeric(length(span.vals)-1)
for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- abs(span.vals[i+1] - span.vals[i]) }

breaks <- which.max(span.vals.var) + 24
model.predicted.t <- span.vals[breaks]

model.predicted <- predict( rasters , model.brt, n.trees=model.brt$n.trees,type="response") 
model.predicted[model.predicted < model.predicted.t] <- 0
model.predicted[model.predicted != 0] <- 1

plot(model.predicted)
model.predicted <- model.predicted - min(getValues(model.predicted),na.rm=T)
model.predicted <- model.predicted / max(getValues(model.predicted),na.rm=T)
plot(model.predicted)

writeRaster(model.predicted,filename=paste0(results.directory,"/Results/Correlative model/Suitability_2090_2100.tif"),format="GTiff",overwrite=T)

# ------------------------------------------------------------------------------------
# Ensemble

y.span <- 2000:2016
ensemble <- list()

for( y in y.span ) {
  
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
  
  model.predicted <- predict( rasters , model.brt , n.trees=model.brt$n.trees,type="response")
  
  span.vals <- plot(model.brt,1, return.grid=TRUE,type="response")$y
  span.vals.var <- numeric(length(span.vals)-1)
  for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- abs(span.vals[i+1] - span.vals[i]) }
  
  breaks <- which.max(span.vals.var) + 8
  model.predicted.t <- span.vals[breaks]
  
  model.predicted[model.predicted < model.predicted.t] <- 0
  model.predicted[model.predicted != 0] <- 1
  
  ensemble <- c(ensemble,list(model.predicted))
  
}

# -------------------------

model.predicted <- calc(stack(ensemble),mean)
plot(model.predicted)
model.predicted <- model.predicted - min(getValues(model.predicted),na.rm=T)
model.predicted <- model.predicted / max(getValues(model.predicted),na.rm=T)
plot(model.predicted)

writeRaster(final.ensemble.sum,filename=paste0(results.directory,"/Results/Correlative model/Suitability_",y.span[1],"_",y.span[length(y.span)],".tif"),format="GTiff",overwrite=T)

final.ensemble.sum <- calc(stack(ensemble),sum)
plot(final.ensemble.sum)
final.ensemble.sum[final.ensemble.sum < 1] <- 0
final.ensemble.sum[final.ensemble.sum >= 1] <- 1
plot(final.ensemble.sum)

writeRaster(final.ensemble.sum,filename=paste0(results.directory,"/Results/Correlative model/Suitability_2016.tif"),format="GTiff",overwrite=T)

plot(final.ensemble.sum)
presences <- read.occurrence.records(presence.records.file,exclude.na=FALSE,exclude.duplicated=FALSE,exclude.low.resolution=FALSE,exclude.out.of.range=FALSE,presence.records.mask.file=NULL,extra.columns=c("Year.of.record","Year.of.comparison","Occurrence.Data"))
presences[,1:2] <- relocate.coordinates.na(presences=presences[,1:2],rasters=shape,bathymetry=NULL,maximum.distance=500,use.species.depth=FALSE)
presences <- presences[which(!is.na(presences$Occurrence.Data)),]
presences.i <- presences[presences$Year.of.record == 2016,]
points(presences.i[presences.i$Occurrence.Data == 1,],col="Black")
points(presences.i[presences.i$Occurrence.Data == 0,],col="Gray")

write.csv(presences.i,file=paste0(results.directory,"/Results/Correlative model/Records_2016.csv"))

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
        if( ( ! FALSE %in% ( r.l <= disp ) ) & ( reconstruction.before[l,ncol(reconstruction.before)] == 1 )  ) { predicted[l] <- 1 }
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

# ----------------------------------------------------------------------------------------------------------------------

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

# --------------------------------------------------------------------------------------------------------

time.series <- data.frame(matrix(NA,ncol=length(vars)+1,nrow=length(2000:2016)))
colnames(time.series) <- c(vars,"Suitable.Habitats")

for( y in 1:length(2000:2016) ) {
  for( v in 1:length(vars)) {
    vals <- getValues(subset( get(vars[v]) ,y))
    time.series[y,v] <- mean(vals[!is.na(vals)])
    if( y == 17) { time.series[y,v] <- mean(vals[!is.na(vals)]) - 0.32 }
  }
}

time.series[,"Suitable.Habitats"] <- reconstruction.area

# --------------------------------------------

# 14 10
# TimeSeriesNitrate

plot(2000:2016,time.series$Suitable.Habitats / 1000, ylab="Suitable habitats (10^3 km2)" , xlab="Time (year)") # x 
lines(2000:2016,time.series$Suitable.Habitats / 1000,lty=2, col="grey")

plot(2000:2016,time.series$Ocean.temperature.Max , ylab="Ocean temperature max. (ºC)" , xlab="Time (year)")
lines(2000:2016,time.series$Ocean.temperature.Max,lty=2, col="grey")

(2000:2016)[sort(time.series$Suitable.Habitats,index.return=TRUE)$ix]
(2000:2016)[sort(time.series$Ocean.temperature.Max,index.return=TRUE,decreasing = TRUE)$ix]

sort(time.series$Suitable.Habitats,index.return=TRUE)
(2000:2016)[sort(time.series$Ocean.temperature.Max,index.return=TRUE,decreasing = TRUE)$ix]

# --------------------------------------------
