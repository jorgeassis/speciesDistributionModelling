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

## -----

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
closeAllConnections()

# setwd("~/Dropbox/Manuscripts/_ Under Revision/Global patterns of wild seaweed productivity")

source("mainFunctions.R")
library(foreach)
library(doParallel)

source("GitGIS/Dependencies/mainFunctions.R")
library(ggnewscale)
library(sf)
library(h3js) # devtools::install_github("saurfang/h3js")
library(h3jsr) # remotes::install_github("obrl-soil/h3jsr")

theme_map <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0.1),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.border = element_blank()
  )


dataLayersDirectory <- "Data/Climate/Present/"
detectCores(all.tests = TRUE, logical = TRUE)
nCores <- 16

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

resultsName <- "Subtidal Browns"
maskDistribution <- "Data/PAModels/allButSeagrasses_Diversity_Present_Reachable.tif"
maskDistribution <- raster(maskDistribution)
maskDistribution[!is.na(maskDistribution)] <- 1

Habitat <- c("Marine Forest")
Level <- c("Subtidal")

dataLayers <- c("LightAtBottom_Pred_Max.tif","OceanTemperature Surface Pred Min.tif","OceanTemperature Surface Pred Max.tif","Salinity Surface Pred Mean.tif","Nitrate Surface Pred Mean.tif","CoastalExposureW_Pred_Max.tif")
dataLayersName <- c("Light","TempMin","TempMax","Salinity","Nitrate","WaveEnergy") # "Ice" "SeaIceCover Surface Pred Mean.tif",
dataLayersMonotonocity <- c(+1,+1,-1,+1,+1,+1)

data.frame(dataLayersName,dataLayersMonotonocity)

dataset <- read.csv("Data/Productivity/Final NPP dataset.csv", stringsAsFactors=FALSE, fileEncoding="latin1", sep=",")
dataset <- dataset[which(dataset$Vegetation.type %in% Habitat & dataset$Level == Level & ! is.na(dataset$Long_converted)),]
nrow(dataset)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

resultsName <- "Intertidal Marine forest"

maskDistribution <- "Data/PAModels/allButSeagrasses_Diversity_Present_Reachable.tif"
maskDistributionIntertidal <- "Data/CoastLine.tif"
maskDistribution <- mask(raster(maskDistribution),raster(maskDistributionIntertidal))
maskDistribution[!is.na(maskDistribution)] <- 1

Habitat <- c("Marine Forest")
Level <- c("Intertidal")

dataLayers <- c("OceanTemperature Surface Pred Min.tif","OceanTemperature Surface Pred Max.tif","Salinity Surface Pred Mean.tif","Nitrate Surface Pred Mean.tif","CoastalExposureW_Pred_Max.tif") # "SeaIceCover Surface Pred Mean.tif",
dataLayersName <- c("TempMin","TempMax","Salinity","Nitrate","CoastalExposure") # "Ice",
dataLayersMonotonocity <- c(+1,-1,+1,+1,0)
data.frame(dataLayersName,dataLayersMonotonocity)

dataset <- read.csv("Data/Productivity/Final NPP dataset.csv", stringsAsFactors=FALSE, fileEncoding="latin1", sep=",")
dataset <- dataset[which(dataset$Vegetation.type %in% Habitat & dataset$Level == Level & ! is.na(dataset$Long_converted)),]
nrow(dataset)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

resultsName <- "Subtidal Red algae"
maskDistribution <- "Data/PAModels/LT2 RedAlgaeBenthicAggregated.tif"
maskDistribution <- raster(maskDistribution)

Habitat <- c("Red algae")
Level <- c("Subtidal")

library(gdalUtils)
maskDistribution <- raster::resample(maskDistribution,subset(rasterLayers,1))
maskDistribution[maskDistribution != 1] <- NA

dataLayers <- c("LightAtBottom_Pred_Max.tif","OceanTemperature Surface Pred Min.tif","OceanTemperature Surface Pred Max.tif","Salinity Surface Pred Mean.tif","Nitrate Surface Pred Mean.tif","CoastalExposureW_Pred_Max.tif")
dataLayersName <- c("Light","TempMin","TempMax","Salinity","Nitrate","CoastalExposure") # 
dataLayersMonotonocity <- c(+1,+1,-1,+1,+1,+1)
data.frame(dataLayersName,dataLayersMonotonocity)

dataset <- read.csv("Data/Productivity/Final NPP dataset.csv", stringsAsFactors=FALSE, fileEncoding="latin1", sep=",")
dataset <- dataset[which(dataset$Vegetation.type %in% Habitat & dataset$Level == Level & ! is.na(dataset$Long_converted)),]
nrow(dataset)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

resultsName <- "Intertidal Red algae"
maskDistribution <- "Data/PAModels/LT2 RedAlgaeBenthicAggregated.tif"
maskDistribution <- raster(maskDistribution)

Habitat <- c("Red algae")
Level <- c("Intertidal")

library(gdalUtils)
maskDistribution <- raster::resample(maskDistribution,subset(rasterLayers,1))
maskDistribution[maskDistribution != 1] <- NA

dataLayers <- c("SeaIceCover Surface Pred Mean.tif","OceanTemperature Surface Pred Min.tif","OceanTemperature Surface Pred Max.tif","Salinity Surface Pred Mean.tif","Nitrate Surface Pred Mean.tif","CoastalExposureW_Pred_Max.tif")
dataLayersName <- c("Ice","TempMin","TempMax","Salinity","Nitrate","CoastalExposure") # 
dataLayersMonotonocity <- c(-1,+1,-1,+1,+1,0)
data.frame(dataLayersName,dataLayersMonotonocity)

dataset <- read.csv("Data/Productivity/Data 6 - Turf split.csv", stringsAsFactors=FALSE, fileEncoding="latin1", sep=",")
dataset <- dataset[which(dataset$Vegetation_category %in% Habitat & dataset$Level %in% Level & ! is.na(dataset$Long_converted)),]
nrow(dataset)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

resultsName <- "Combined Red algae"
maskDistribution <- "Data/PAModels/LT2 RedAlgaeBenthicAggregated.tif"
maskDistribution <- raster(maskDistribution)

Habitat <- c("Red algae")
Level <- c("Subtidal", "Intertidal")

library(gdalUtils)
maskDistribution <- raster::resample(maskDistribution,subset(rasterLayers,1))
maskDistribution[maskDistribution != 1] <- NA

dataLayers <- c("LightAtBottom_Pred_Max.tif","OceanTemperature Surface Pred Min.tif","OceanTemperature Surface Pred Max.tif","Salinity Surface Pred Mean.tif","Nitrate Surface Pred Mean.tif","CoastalExposureW_Pred_Max.tif")
dataLayersName <- c("Light","TempMin","TempMax","Salinity","Nitrate","CoastalExposure") # 
dataLayersMonotonocity <- c(+1,+1,-1,+1,+1,+1)
data.frame(dataLayersName,dataLayersMonotonocity)

dataset <- read.csv("Data/Productivity/Data 6 - Turf split.csv", stringsAsFactors=FALSE, fileEncoding="latin1", sep=",")
dataset <- dataset[which(dataset$Vegetation_category %in% Habitat & dataset$Level %in% Level & ! is.na(dataset$Long_converted)),]
nrow(dataset)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

resultsName <- "Subtidal Turfs"
maskDistribution <- "Data/PAModels/LT2 RedAlgaeBenthicAggregated.tif"
maskDistribution <- raster(maskDistribution)

Habitat <- c("Algal turfs")
Level <- c("Subtidal")

library(gdalUtils)
maskDistribution <- raster::resample(maskDistribution,subset(rasterLayers,1))
maskDistribution[maskDistribution != 1] <- NA

dataLayers <- c("LightAtBottom_Pred_Max.tif","OceanTemperature Surface Pred Min.tif","OceanTemperature Surface Pred Max.tif","Salinity Surface Pred Mean.tif","Nitrate Surface Pred Mean.tif","CoastalExposureW_Pred_Max.tif")
dataLayersName <- c("Light","TempMin","TempMax","Salinity","Nitrate","CoastalExposure") # 
dataLayersMonotonocity <- c(+1,+1,-1,+1,+1,+1)
data.frame(dataLayersName,dataLayersMonotonocity)

dataset <- read.csv("Data/Productivity/Final NPP dataset.csv", stringsAsFactors=FALSE, fileEncoding="latin1", sep=",")
dataset <- dataset[which(dataset$Vegetation.type %in% Habitat & dataset$Level == Level & ! is.na(dataset$Long_converted)),]
nrow(dataset)

## -----------------------------------------------------
## -----------------------------------------------------
## -----------------------------------------------------

rasterLayers <- list.files(dataLayersDirectory,pattern = "tif",recursive = TRUE, full.names = TRUE)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
names(rasterLayers) <- dataLayersName
shape <- subset(rasterLayers,1)
datasetCoords <- data.frame(Lon=dataset$Long_converted,Lat=dataset$Lat_converted)
datasetCoords <- relocateNACoords(datasetCoords,rasterLayers)

# ------------------

scaled_logit <- function(x){ x } # log(x+1)
inv_scaled_logit <- function(x){ x  } # exp(x)-1
modelDataset <- data.frame(resp=scaled_logit(dataset$Avg_ann_prod_kg_C_m2_y),datasetCoords,method=dataset$Prod_method_general)

# ------------------
# ------------------

datasetEnvironment <- raster::extract(rasterLayers,data.frame(x=modelDataset$Lon,y=modelDataset$Lat))
modelDataset <- data.frame(modelDataset,datasetEnvironment)
modelDataset$method <- as.factor(modelDataset$method)
dataLayersMonotonocity <- c(0,dataLayersMonotonocity)
modelDataset <- modelDataset[,-which(colnames(modelDataset) %in% c("Lon","Lat"))]

set.seed(1)
kFolds <- kfold(1:nrow(modelDataset),k=6)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

brtLearning <- c(0.01,0.001,0.0001)
brtTreeDepth <- c(1,2,3,4,5,6)
brtBagFraction <- seq(0.1,0.9,0.1) # c(seq(0.1,0.9,0.1),2,3)
brtnTrees <- seq(50,1000,50)

comb = expand.grid(cv.k = 1:length(unique(kFolds)), learning.complex=brtLearning,tree.depth=brtTreeDepth , bag=brtBagFraction,trees=brtnTrees )

cl <- parallel::makeCluster(nCores)
registerDoParallel(cl)

cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .export=c("Dsquared"), .packages = c("gbm","dismo","ENMeval")) %dopar% {
  
  cv <- comb[c,1]
  l.rate <- comb[c,2]
  tree.c <- comb[c,3]
  bag <- comb[c,4]
  trees <- comb[c,5]
  
  train.dataset <- modelDataset
  test.dataset <- modelDataset
  
  # train.dataset <- modelDataset[kFolds != cv,]
  # test.dataset <- modelDataset[kFolds == cv,]
  
  model <- NULL
  set.seed(1)
  tryCatch(       
    model <- gbm(formula = resp ~ ., distribution = "gaussian",
                 data = train.dataset, var.monotone = dataLayersMonotonocity, n.trees = trees,
                 interaction.depth = tree.c, n.minobsinnode = 10, shrinkage = l.rate,
                 bag.fraction = bag, train.fraction = 1, cv.folds = 0,
                 keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL,
                 n.cores = NULL) , error=function(e) { model <<- NULL })
  
  if(!is.null(model)) {
    
    num.tress <- trees
    observed <- inv_scaled_logit(test.dataset[,1])
    predicted <- inv_scaled_logit(predict( model , test.dataset[,-1] , n.trees=num.tress,type="response"))
    
    model.deviance <- Dsquared(glm(predicted~observed,family=gaussian()))
    
    predicted.accuracy <- data.frame( cv.round=cv,tree.c=tree.c,l.rate=l.rate,bag=bag,trees=trees,accuracy=model.deviance)
    return(predicted.accuracy)
  }
  
}

stopCluster(cl); rm(cl) ; gc(reset=TRUE)

# ------------

best.model <- aggregate(list(cvIndex=cv.accuracy[,"accuracy"]), by = list(bag=cv.accuracy$bag,tree.c=cv.accuracy$tree.c,l.rate=cv.accuracy$l.rate,trees=cv.accuracy$trees), mean)
best.model[which.max(best.model$cvIndex),]

tree.c <- best.model[which.max(best.model$cvIndex),"tree.c"]
l.rate <- best.model[which.max(best.model$cvIndex),"l.rate"]
bag <- best.model[which.max(best.model$cvIndex),"bag"]
trees <- best.model[which.max(best.model$cvIndex),"trees"]

set.seed(1)
model <- gbm(formula = resp ~ ., distribution = "gaussian",
             data = modelDataset, var.monotone = dataLayersMonotonocity, n.trees = trees,
             interaction.depth = tree.c, n.minobsinnode = 1, shrinkage = l.rate,
             bag.fraction = bag, train.fraction = 1, cv.folds = 0,
             keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL,
             n.cores = NULL)

observed <- modelDataset[,1]
predicted <- predict( model , modelDataset[,-1] , n.trees=trees,type="response")
Dsquared(glm(predicted~observed,family=gaussian()))

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=12) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) ,
                   axis.text.x =element_text(size=11, margin = margin(t = 8, r = 0, b = 0, l = 0)),
                   axis.text.y =element_text(size=11, margin = margin(t = 0, r = 8, b = 0, l = 0)))

dataPlot <- data.frame(observed=observed,predicted=predicted)
plot1 <- ggplot(dataPlot, aes(x=observed, y=predicted)) +
  geom_point(size=2,shape=21, color="black", fill="black",alpha = 0.25) +
  geom_smooth(method=lm , color="black", fill="#B5CAE5", se=TRUE,size=0.5) +
  ylab("Predicted (productivity)") + xlab("Observed (productivity)") + mainTheme +
  annotate(geom="text", x=min(dataPlot$observed), y=max(dataPlot$predicted)-0.2, label=paste0("Dev. explained: ", round(Dsquared(glm(predicted~observed,family=gaussian())),digits=3),"\n","Mean error: ",round(mean(predicted-observed),digits=4),"\n","n: ",nrow(dataset)),size=4.5,family="Helvetica", color = "#22211d",hjust = 0)
plot1  

pdf(paste0("Results/",resultsName,"Fit.pdf"), width = 10, height = 9 )
plot1
dev.off()
write.csv(data.frame(observed=observed,predicted=predicted),file=paste0("Results/",resultsName,"Fit.csv"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Variable Contribution

sum.brt <- summary(model) ; colnames(sum.brt) <- c("Predictor","Contribution")
sum.brt <- sum.brt[sort(sum.brt$Contribution,index.return=TRUE)$ix,]

relativeImportancePlot <- ggplot(sum.brt) +
  geom_bar( aes(x= reorder(Predictor, Contribution) , y=Contribution), stat="identity", fill="black", alpha=0.5) +
  coord_flip() + theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
    panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
  ) + labs(x = "Predictor") + 
  labs(y = "Relative importance (%)") + geom_hline(aes(yintercept=5 ),color="Black", linetype="dashed", size=0.3) +
  annotate("text", y = 5 + 1 , x = 1 , label = "5%" , hjust = 0) + geom_hline(aes(yintercept=0 ),color="Gray", size=0.3)

pdf(paste0("Results/",resultsName,"RelContrib.pdf"), width = 12, height = 9 )
relativeImportancePlot
dev.off()
write.csv(sum.brt,file=paste0("Results/",resultsName,"RelContrib.csv"))

# -------------

sortVars <- as.character(relativeImportancePlot$data$Predictor[sort(relativeImportancePlot$data$Contribution,decreasing = T, index.return=T)$ix])

pdf(paste0("Results/",resultsName,"ParPlots.pdf"), width = 9, height = 12 )
par(mfrow = c(4,2))
for( var.name in sortVars ) {
  var.n <- which(model$var.names == var.name)
  modelPlot(model,rasterLayers,distribution="gaussian",variable.to.plot=var.n,print.limiting=FALSE,auto.limiting=FALSE,distribution.threshold=NULL,distribution.predicted=NULL,val.limiting=NULL,export.data.to.plot=FALSE,plot.x.lim=NULL,plot.y.lim=NULL)
}
dev.off()
par(mfrow=c(1,1))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Make predictions

maskCoastline <- raster("Data/CoastLine.tif")
predictionName <- "Present" # 2100RCP85
dataLayersDirectory <- "Data/Climate/Present/" # RCP26 RCP85 Present

rasterLayers <- list.files(dataLayersDirectory,pattern = "tif",recursive = TRUE, full.names = TRUE)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
names(rasterLayers) <- dataLayersName
shape <- subset(rasterLayers,1)

maskDistribution[maskDistribution != 1] <- NA
maskDistribution <- raster::resample(maskDistribution,rasterLayers)

rasterLayersPredict <- mask(rasterLayers,maskDistribution)
rasterLayersPredictDF <- as.data.frame(rasterLayersPredict,xy=T)
rasterLayersPredictDF <- rasterLayersPredictDF[!is.na(rasterLayersPredictDF$TempMax),]
rasterLayersPredictDFCoords <- rasterLayersPredictDF[,c("x","y")]

methodsInModel <- unique(modelDataset$method)
predictionDF <- matrix(NA,ncol=length(methodsInModel),nrow=nrow(rasterLayersPredictDFCoords))
predictionDF.conf <- matrix(NA,ncol=length(methodsInModel),nrow=nrow(rasterLayersPredictDFCoords))

for( m in 1:length(methodsInModel)) {
  
  m.name <- methodsInModel[m]
  rasterLayersPredictDF.i <- data.frame(rasterLayersPredictDF,method=m.name)
  predictionDF[,m] <- inv_scaled_logit( predict( model , rasterLayersPredictDF.i , n.trees=trees,type="response") )
  
  model.confidence.data <- data.frame(observed=modelDataset[,-1],predicted=predict( model , modelDataset[,-1] , n.trees=trees,type="response"))
  model.confidence <- lm( observed~predicted, model.confidence.data)
  
  prediction.conf.90 <- predict( model.confidence , data.frame(predicted=predictionDF[,m]) , interval="confidence")
  predictionDF.conf[,m] <- as.numeric(prediction.conf.90[,3] - prediction.conf.90[,2])
  
}

prediction <- apply(predictionDF,1,mean)
predictionMap <- rasterFromXYZ(cbind(rasterLayersPredictDFCoords,prediction))
predictionMap[predictionMap < 0] <- 0
predictionMap <- extend(predictionMap,maskCoastline)
predictionMap <- mask(predictionMap,maskCoastline)
writeRaster(predictionMap,filename=paste0("Results/",resultsName,predictionName,".tif"),format="GTiff",overwrite=T)

prediction.sd <- apply(predictionDF,1,sd) / length(methodsInModel)
predictionMap.sd <- rasterFromXYZ(cbind(rasterLayersPredictDFCoords,prediction.sd))
predictionMap.sd[predictionMap.sd < 0] <- 0
predictionMap.sd <- extend(predictionMap.sd,maskCoastline)
predictionMap.sd <- mask(predictionMap.sd,maskCoastline)
writeRaster(predictionMap.sd,filename=paste0("Results/",resultsName,predictionName,"SD Methods.tif"),format="GTiff",overwrite=T)

prediction.conf <- apply(predictionDF.conf,1,mean); head(prediction.conf)
predictionMap.conf <- rasterFromXYZ(cbind(rasterLayersPredictDFCoords,prediction.conf))
predictionMap.conf[predictionMap.conf < 0] <- 0
predictionMap.conf <- extend(predictionMap.conf,maskCoastline)
predictionMap.conf <- mask(predictionMap.conf,maskCoastline)
writeRaster(predictionMap.conf,filename=paste0("Results/",resultsName,predictionName,"Confidence.tif"),format="GTiff",overwrite=T)

predictionMaps <- list(predictionMap,predictionMap.sd,predictionMap.conf)
predictionMapsNames <- c("","SD Methods","Range of confidence")

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

mapExtent = c(xmin = -180, ymin = -90, xmax = 180, ymax = 90)
worldMap <- ne_countries(scale = 10, returnclass = "sf")
worldMapCoordRef <- crs(worldMap)

mainGlobalMap <- ggplot() + 
  geom_sf(data = worldMap,fill="#E6E6E6", colour = "#E6E6E6" , size=0.1 ) +
  theme(axis.ticks=element_blank()) + ylab("") + xlab("") +
  theme_map + coord_sf(crs = worldMapCoordRef)

## ------------
## ------------

for( map in 1:length(predictionMaps)) {
  
  rasterMap <- predictionMaps[map][[1]]
  fileName <- paste0("Results/",resultsName," ",predictionName," ",predictionMapsNames[map],".pdf")
  
  ## --------------
  
  # https://github.com/uber/h3/blob/master/docs/core-library/restable.md
  resolutionH3 <- 3
  
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
  rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
  rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex)
  rasterMapDF.polygons <- st_sf(rasterMapDF.polygons)
  
  scaled_logit <- function(x){ x } # log(x+1)
  inv_scaled_logit <- function(x){ x } # exp(x)-1
  
  # rasterMapDF$val[rasterMapDF$val > max(modelDataset$resp) ] <- max(modelDataset$resp)
  rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
  rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
  rasterMapDF.polygons <- st_wrap_dateline(rasterMapDF.polygons, options = "WRAPDATELINE=YES", quiet = TRUE)
  rasterMapDF.polygons <- sf:::as_Spatial(rasterMapDF.polygons)
  rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF.polygons$value))
  rasterMapDF.polygons$hex <- as.character(rasterMapDF.polygons$hex)
  rasterMapDF.polygons <- fortify(rasterMapDF.polygons)
  rasterMapDF.polygons <- rasterMapDF.polygons[rasterMapDF.polygons$lat >= extent(worldMap)[3],]
  rasterMapDF.polygons$value <- sapply(as.numeric(rasterMapDF.polygons$id), function(x) { rasterMapDF[x,"val"] })
  rasterMapDF.polygons$value[rasterMapDF.polygons$value < 0] <- 0
  
  minLegend <- min(rasterMapDF.polygons$value); minLegend
  maxLegend <- max(rasterMapDF.polygons$value); maxLegend
  
  if( Habitat == "Algal turfs" ) { rasterMapDF.polygons[rasterMapDF.polygons$lat >= 57.5 | rasterMapDF.polygons$lat <= -37.5 ,"value" ] <- NA }
  
  if(map == 3) {
    
    minLegend <- 0
    maxLegend <- 0.9
    myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
    
  }
  
  if(map == 2) {
    
    minLegend <- 0
    maxLegend <- 1.6
    myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
    
  }
  
  if(map == 1) {
    
    minLegend <- 0
    maxLegend <- 5
    myColors <- c("#00A8FC","#C6FCB6","#ECE15F","#F47D75","#B84F9D","#5F52A3")
    
  }
  
  #----------------------
  
  modelDataset.i <- modelDataset
  modelDataset.i <- data.frame(modelDataset.i,Lon=dataset$Long_converted,Lat=dataset$Lat_converted)
  modelDataset.i$resp[modelDataset.i$resp >= maxLegend ] <- maxLegend
  
  #----------------------
  
  plot1 <- mainGlobalMap +
    
    scale_colour_gradientn(colours = myColors, na.value = NA, values = seq(0,1,length.out=5), n.breaks = 6,labels=seq((minLegend),(maxLegend), length.out=6), guide = "colourbar", aesthetics = "fill", limits=c(minLegend,maxLegend) ) +
    geom_polygon(data = rasterMapDF.polygons, aes(x = long, y = lat,group=id, fill=value), colour = NA, size = 0) +
    
    new_scale("fill") +
    # scale_colour_gradientn(colours = myColors, values = seq(0,1,length.out=5), n.breaks = 6,labels=round(seq((min(modelDataset$resp)),(max(modelDataset$resp)), length.out=6), digits=2), guide = "colourbar", aesthetics = "fill", limits=c(minLegend,maxLegend) ) +
    
    scale_colour_gradientn(colours = myColors, na.value = NA, values = seq(0,1,length.out=6), n.breaks = 6,labels=round(seq((minLegend),(maxLegend), length.out=6), digits=2), guide = FALSE, aesthetics = "colour", limits=c(minLegend,6) ) +
    
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm'),
          legend.background = element_rect(fill = "#FFFFFF", color = NA)) + theme(legend.title=element_blank())
  
  if(map == 1) {
    
    plot1 <- plot1 + geom_point(data=modelDataset.i, aes(x=Lon, y=Lat,  colour=resp), size=3.5,colour="black", alpha=1) + # create a large black dot to serve as the outline
      geom_point(data=modelDataset.i, aes(x=Lon, y=Lat,  colour=resp), size=2.5, alpha=1)
    #scale_size_continuous(range = c(2,5))+
    
  }
  
  pdf(file=fileName,width=14,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
}

# ----------------------------
# ----------------------------

# Marine EcoR Estimates

ecoR <- shapefile("Data/marine_ecoregions.shp")

rasterMap <- predictionMaps[1][[1]]

# ---------

resRealm <- data.frame(Realm=unique(ecoR$REALM),minProductivity=NA,maxProductivity=NA,averageProductivity=NA,confidence5Productivity=NA,confidence95Productivity=NA)

for(i in 1:nrow(resRealm)) {
  ecoR.i <- ecoR[ecoR$REALM == unique(ecoR$REALM)[i],]
  rasterMap.i <- crop(rasterMap,ecoR.i)
  rasterMap.i <- mask(rasterMap.i,ecoR.i)
  resRealm[i,"minProductivity"] <- (cellStats(rasterMap.i,min))
  resRealm[i,"maxProductivity"] <- (cellStats(rasterMap.i,max))
  resRealm[i,"averageProductivity"] <- (cellStats(rasterMap.i,mean))
  
  rasterMap.i <- crop(predictionMap.conf,ecoR.i)
  rasterMap.i <- mask(rasterMap.i,ecoR.i)
  
  resRealm[i,"confidence5Productivity"] <- resRealm[i,"averageProductivity"] - (cellStats(rasterMap.i,mean))
  resRealm[i,"confidence95Productivity"] <- resRealm[i,"averageProductivity"] + (cellStats(rasterMap.i,mean))
  
}

resRealm[resRealm == Inf] <- 0
resRealm[resRealm == -Inf] <- 0
resRealm[ is.na(resRealm)] <- 0
write.csv(resRealm,file=gsub(".pdf"," Realm.csv",fileName))

# ---------

resProvince <- data.frame(Province=unique(ecoR$PROVINCE),minProductivity=NA,maxProductivity=NA,averageProductivity=NA,confidence5Productivity=NA,confidence95Productivity=NA)

for(i in 1:nrow(resProvince)) {
  
  ecoR.i <- ecoR[ecoR$PROVINCE == unique(ecoR$PROVINCE)[i],]
  error <- FALSE
  tryCatch(
    expr = {
      rasterMap.i <- crop(rasterMap,ecoR.i)
      rasterMap.i <- mask(rasterMap.i,ecoR.i)
    },
    error = function(e){
      error <<- FALSE
    }
  )
  
  if(error) { next }
  
  resProvince[i,"minProductivity"] <- (cellStats(rasterMap.i,min))
  resProvince[i,"maxProductivity"] <- (cellStats(rasterMap.i,max))
  resProvince[i,"averageProductivity"] <- (cellStats(rasterMap.i,mean))
  
  rasterMap.i <- crop(predictionMap.conf,ecoR.i)
  rasterMap.i <- mask(rasterMap.i,ecoR.i)
  
  resProvince[i,"confidence5Productivity"] <- resProvince[i,"averageProductivity"] - (cellStats(rasterMap.i,mean))
  resProvince[i,"confidence95Productivity"] <- resProvince[i,"averageProductivity"] + (cellStats(rasterMap.i,mean))
  
}

resProvince[resProvince == Inf] <- 0
resProvince[resProvince == -Inf] <- 0
resProvince[ is.na(resProvince)] <- 0
write.csv(resProvince,file=gsub(".pdf"," Province.csv",fileName))

# ---------

resEcoregion <- data.frame(Ecoregion=unique(ecoR$ECOREGION),minProductivity=NA,maxProductivity=NA,averageProductivity=NA,confidence5Productivity=NA,confidence95Productivity=NA)

for(i in 1:nrow(resEcoregion)) {
  ecoR.i <- ecoR[ecoR$ECOREGION == unique(ecoR$ECOREGION)[i],]
  error <- FALSE
  tryCatch(
    expr = {
      rasterMap.i <- crop(rasterMap,ecoR.i)
      rasterMap.i <- mask(rasterMap.i,ecoR.i)
    },
    error = function(e){
      error <<- FALSE
    }
  )
  
  if(error) { next }
  
  resEcoregion[i,"minProductivity"] <- (cellStats(rasterMap.i,min))
  resEcoregion[i,"maxProductivity"] <- (cellStats(rasterMap.i,max))
  resEcoregion[i,"averageProductivity"] <- (cellStats(rasterMap.i,mean))
  
  rasterMap.i <- crop(predictionMap.conf,ecoR.i)
  rasterMap.i <- mask(rasterMap.i,ecoR.i)
  
  resEcoregion[i,"confidence5Productivity"] <- resEcoregion[i,"averageProductivity"] - (cellStats(rasterMap.i,mean))
  resEcoregion[i,"confidence95Productivity"] <- resEcoregion[i,"averageProductivity"] + (cellStats(rasterMap.i,mean))
  
}

resEcoregion[resEcoregion == Inf] <- 0
resEcoregion[resEcoregion == -Inf] <- 0
resEcoregion[ is.na(resEcoregion)] <- 0
write.csv(resEcoregion,file=gsub(".pdf"," Ecoregion.csv",fileName))

# ----------------------------

# 0-10 maxT (summer isoterm) - Polar regions
# 10-15 - Cold Temperate
# 15-25- Warm temperate
# >25 - Tropical

maxTemp <- subset(rasterLayers,3)
thremalregion <- data.frame(Ecoregion=c("Polar regions","Cold Temperate","Warm temperate","Tropical"),minProductivity=NA,maxProductivity=NA,averageProductivity=NA,confidence5Productivity=NA,confidence95Productivity=NA)

for(i in 1:nrow(thremalregion)) {
  
  maxTemp.i <- maxTemp
  
  if( i == 1) { maxTemp.i[maxTemp.i > 10] <- NA; maxTemp.i[!is.na(maxTemp.i)] <- 1 }
  if( i == 2) { maxTemp.i[maxTemp.i < 10] <- NA; maxTemp.i[maxTemp.i > 15] <- NA;  maxTemp.i[!is.na(maxTemp.i)] <- 1 }
  if( i == 3) { maxTemp.i[maxTemp.i < 15] <- NA; maxTemp.i[maxTemp.i > 25] <- NA;  maxTemp.i[!is.na(maxTemp.i)] <- 1 }
  if( i == 4) { maxTemp.i[maxTemp.i < 25] <- NA; maxTemp.i[!is.na(maxTemp.i)] <- 1 }
  
  rasterMap.i <- crop(rasterMap,maxTemp.i)
  rasterMap.i <- extend(rasterMap.i,maxTemp.i)
  rasterMap.i <- mask(rasterMap.i,maxTemp.i)
  
  thremalregion[i,"minProductivity"] <- (cellStats(rasterMap.i,min))
  thremalregion[i,"maxProductivity"] <- (cellStats(rasterMap.i,max))
  thremalregion[i,"averageProductivity"] <- (cellStats(rasterMap.i,mean))
  
  rasterMap.i <- crop(predictionMap.conf,maxTemp.i)
  rasterMap.i <- mask(rasterMap.i,maxTemp.i)
  
  thremalregion[i,"confidence5Productivity"] <- thremalregion[i,"averageProductivity"] - (cellStats(rasterMap.i,mean))
  thremalregion[i,"confidence95Productivity"] <- thremalregion[i,"averageProductivity"] + (cellStats(rasterMap.i,mean))
  
}

thremalregion[thremalregion == Inf] <- 0
thremalregion[thremalregion == -Inf] <- 0
thremalregion[ is.na(thremalregion)] <- 0
write.csv(thremalregion,file=gsub(".pdf"," Thermalregion.csv",fileName))

# ----------------------------
