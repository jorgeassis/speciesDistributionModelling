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

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("Dependencies/mainFunctions.R")
source("0. Config.R")

# ---------------------
# ---------------------

if( length(list.dirs(resultsDirectory)) == 0 ) { dir.create(resultsDirectory) }
if( length(list.dirs(paste0(resultsDirectory,"/RData/"))) == 0 ) { dir.create(paste0(resultsDirectory,"/RData/")) }

# ------------------------------------------------------------------------------------
# Occurrence records

occurrenceRecords <- read.table(dataRecordsFile,sep=";",header=T)
occurrenceRecords <- data.frame(Lon=occurrenceRecords$Lon,Lat=occurrenceRecords$Lat)

# ---------------------------------
# ---------------------------------

shape <- listAllFiles(dataLayersDirectory,dataLayersFileType)[1]
shape <- raster(shape)

library(fasterize)

occurrenceRecords <- shapefile("../Data/Records/SimplifiedSingle2")
occurrenceRecordsSf <- st_as_sf(occurrenceRecords)
occurrenceRecordsRaster <- fasterize( occurrenceRecordsSf, shape, field = NULL, background = NA_real_)
occurrenceRecordsRasterPts <- Which(occurrenceRecordsRaster==1,cells=TRUE)
occurrenceRecordsRasterPts <- xyFromCell(occurrenceRecordsRaster,occurrenceRecordsRasterPts)

occurrenceRecords <- occurrenceRecordsRasterPts
occurrenceRecords <- data.frame(occurrenceRecords)
colnames(occurrenceRecords) <- c("Lon","Lat")
save(occurrenceRecords,file=paste0( resultsDirectory, "/RData/occurrenceRecordsInit.RData"))

# occurrenceRecordsFile <- read.csv(dataRecordsFile,se=";")

# ------------------
# Subset records  outside known distributions [1]

# occurrenceRecords <- occurrenceRecords[occurrenceRecords$Lat > 55 ,]
# occurrenceRecords <- occurrenceRecords[occurrenceRecords$Lat < 60 ,]
# occurrenceRecords <- occurrenceRecords[occurrenceRecords$Lon > 10 ,]
# occurrenceRecords <- occurrenceRecords[occurrenceRecords$Lon < 30 ,]

plotMap(occurrenceRecords,4,"black")

# ------------------
# Remove records outside known distributions [2]

# defineRegion(occurrenceRecords,"Lon","Lat")
# occurrenceRecords <- selectRecords(occurrenceRecords,"Lon","Lat")

plotMap(occurrenceRecords,4,"black")

# ------------------------------------------------------------------------------------
# Environmental layers

rasterLayers <- listAllFiles(dataLayersDirectory,dataLayersFileType)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
rasterLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer,minDepth,maxDepth,intertidalRegion)
names(rasterLayers) <- dataLayersName
names(rasterLayers) 

# ------------------------------------------------------------------------------------
# Drop no variation layers

cave <- function(x) { length(unique(x)) / length(x) }
randomLocations <- Which(!is.na(subset(rasterLayers,1)),cells=TRUE)
randomLocations <- xyFromCell( subset(rasterLayers,1) , sample(randomLocations,min(length(randomLocations),1000),replace=FALSE))
varRasterLayers <- which( apply( raster::extract(rasterLayers,randomLocations) ,2,cave) > 0.025)
rasterLayers <- subset(rasterLayers,varRasterLayers)

# ------------------------------------------------------------------------------------
# Relocate records to closest cell (those falling on land / unlikely depth)

occurrenceRecords <- relocateNACoords(occurrenceRecords,rasterLayers,relocateType,relocateSpeciesDistance,relocateSpeciesDepth,nCores)

# ------------------
# Remove deep records outside known depth range

# occurrenceRecordsDepths <- data.frame(Depth=abs( raster::extract(raster(bathymetryDataLayer),occurrenceRecords) ))
# occurrenceRecordsToKeep <- which( occurrenceRecordsDepths[,1] <= maxDepth & ! is.na(occurrenceRecordsDepths[,1]) )
# occurrenceRecords <- occurrenceRecords[ occurrenceRecordsToKeep,]

# ------------------------------------------------------------------------------------
# Determine autocorrelation 

distanceUncorr <- data.frame(Predictor=names(rasterLayers),Distance=NA)
for( i in 1:length(names(rasterLayers))) {
  distanceUncorr[i,2] <- spatialAutocorrelation(occurrenceRecords=occurrenceRecords,subset(rasterLayers,i),autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationSignif)
}

distanceUncorrPlot <- ggplot(distanceUncorr[sort(distanceUncorr[,2],decreasing = TRUE,index.return=TRUE)$ix,]) +
  geom_bar( aes(x= reorder(Predictor, Distance) , y=Distance), stat="identity", fill="black", alpha=0.5) +
  coord_flip() + theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
    panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
  ) + labs(x = "Predictor") + 
  labs(y = "Spatial correlation (km)") + geom_hline(aes(yintercept=round(mean(distanceUncorr[,2]),digits=2)   ),color="Black", linetype="dashed", size=0.3) +
  annotate("text", y = round(mean(distanceUncorr[,2]),digits=2) + 2 , x = 1 , label = paste0( round(mean(distanceUncorr[,2]),digits=2)," Km") , hjust = 0)

distanceUncorrPlot

pdf(file = paste0(resultsDirectory,"/SpatialAutocorrelation.pdf"), width=12, height=8 )
distanceUncorrPlot
dev.off()

meanCorrDistance <- mean(distanceUncorr[,2])
maxCorrDistance <- max(distanceUncorr[,2])

occurrenceRecordsSubsters <- unique(cellFromXY(subset(rasterLayers,1), occurrenceRecords))
occurrenceRecords <- xyFromCell(subset(rasterLayers,1), occurrenceRecordsSubsters)
occurrenceRecords <- spatialThinning(occurrenceRecords,rasterLayers,meanCorrDistance,verbose=FALSE)
nrow(occurrenceRecords)
plotMap(occurrenceRecords,3,"black")

save(meanCorrDistance,maxCorrDistance,file=paste0(resultsDirectory,"/RData/","corrDistances.RData"))
gc()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Generate Pseudo-absence

polygonPAType <- "exclude" # "exclude" "include"
polygonPA <- drawPolygon(rasterLayers,occurrenceRecords,"pa")
absences <- pseudoAbsences(occurrenceRecords,rasterLayers,patchContinuity=FALSE,polygonPA,polygonPAType,paMindistance,paType,paBiasKernelSurface,paBiasKernelSurfaceProb,paRatio,paEnvironmentStratification,plotFigure=TRUE)

# ------------------------------------------------------------------------------------
# Generate speciesData

speciesData <- data.frame(PA=c(rep(1,nrow(occurrenceRecords)),rep(0,nrow(absences))) , rbind(occurrenceRecords,absences)  )

# ------------------------------------------------------------------------------------
# Save records

write.table(speciesData,file=paste0(substring( dataRecordsFile , 1 , as.numeric(gregexpr("/",dataRecordsFile)[[1]])[length(as.numeric(gregexpr("/",dataRecordsFile)[[1]]))] ) , "speciesData.csv") , quote = FALSE, row.names = FALSE , sep=";")
write.table(occurrenceRecords,file=gsub(".csv"," Processed.csv",dataRecordsFile) , quote = FALSE, row.names = FALSE , sep=";")
write.table(absences,file=gsub(".csv"," pseudoAbsences.csv",dataRecordsFile) , quote = FALSE, row.names = FALSE , sep=";")

# ------------------------------------------------------------------------------------
# Maps and figures

worldMap <- ne_download(scale = "medium", returnclass = "sf")
worldMap.i <- st_crop(worldMap, extent(c(min(speciesData$Lon),max(speciesData$Lon),min(speciesData$Lat),max(speciesData$Lat)))+c(-5,5,-5,5))

occurrenceRecordsMap <- ggplot(data = worldMap.i) + geom_sf(color = NA, fill = "Gray") +
  geom_point(data = occurrenceRecords, aes(x = Lon, y = Lat), size = 0.5, color="red") + xlab("Longitude") + ylab("Latitude")

pdf(file = paste0(resultsDirectory,"/occurrenceRecordsMapGrid.pdf"), width=12 )
occurrenceRecordsMap
dev.off()

# --------------------------

occurrenceRecordsMap <- ggplot(data = worldMap.i) + geom_sf(color = NA, fill = "Gray") +
  geom_point(data = occurrenceRecords, aes(x = Lon, y = Lat), size = 0.5, color="red") + xlab("Longitude") + ylab("Latitude")

pdf(file = paste0(resultsDirectory,"/occurrenceRecordsMap.pdf"), width=12 )
occurrenceRecordsMap
dev.off()

# --------------------------

occurrenceRecordsDepths <- data.frame(Depth=abs( raster::extract(raster(bathymetryDataLayer),occurrenceRecords) ))
occurrenceRecordsDepthsPlot <- ggplot(occurrenceRecordsDepths, aes(x=Depth)) + geom_histogram(color="black", fill="#CACACA", bins=ifelse(nrow(occurrenceRecordsDepths) > 50,50,nrow(occurrenceRecordsDepths))) + geom_vline(aes(xintercept=quantile(Depth,probs=0.95,na.rm=T)  ),color="Black", linetype="dashed", size=0.3) + xlab("Depth (m)") + ylab("Number of records") 
occurrenceRecordsDepthsPlot <- occurrenceRecordsDepthsPlot + annotate("text", y = max(ggplot_build(occurrenceRecordsDepthsPlot)$data[[1]]$count) , x = 2 + quantile(occurrenceRecordsDepths$Depth,probs=0.95,na.rm=T) , label = paste0(round(quantile(occurrenceRecordsDepths$Depth,probs=0.95,na.rm=T),digits=2)," m"), hjust = 0)

pdf(file = paste0(resultsDirectory,"/occurrenceRecordsDepthsPlot.pdf"), width=12, height=8 )
occurrenceRecordsDepthsPlot
dev.off()

# --------------------------

occurrenceRecordsMap <- ggplot(data = worldMap.i) + geom_sf(color = NA, fill = "Gray") +
  geom_point(data = as.data.frame(absences), aes(x = Lon, y = Lat), size = 0.5, color="Black") + xlab("Longitude") + ylab("Latitude") +
  geom_point(data = occurrenceRecords, aes(x = Lon, y = Lat), size = 0.5, color="red") 

pdf(file = paste0(resultsDirectory,"/absenceRecordsMapGrid.pdf"), width=12 )
occurrenceRecordsMap
dev.off()

# --------------------------

gc()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
