## -------------------------------------------------------------------------------
## --------------------------------------------------
## --------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("Dependencies/mainFunctions.R")

nCores <- 12

tempFolder <- "/Volumes/Jellyfish/Temp"

modelsDirectory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Projected climate changes threaten marine forests biodiversity at global scales/Results/"
resultsFolder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Projected climate changes threaten marine forests biodiversity at global scales/Results/"

## ----------------

resultsFolder <- paste0(resultsFolder,"/Summary/")
if(! dir.exists(paste0(resultsFolder))) { dir.create(paste0(resultsFolder), recursive = T) }

## ----------------
## Get layers

layersToCalcSubseter <- NULL # "speciesRichness"
layersToCalcDiscard <- c("Change","Reachable") # NULL

layersToCalc <- list.files(modelsDirectory, full.names = T, recursive = T, pattern = ".tif")
if(!is.null(layersToCalcSubseter)) { layersToCalc <- layersToCalc[ which(apply(sapply(layersToCalcSubseter, function(x) { (grepl(x,layersToCalc)) }  ),1,sum) == length(layersToCalcSubseter)) ]  }
if(!is.null(layersToCalcDiscard)) { layersToCalc <- layersToCalc[ which(apply(sapply(layersToCalcDiscard, function(x) { (! grepl(x,layersToCalc)) }  ),1,sum) == length(layersToCalcDiscard)) ] }

layersToCalcNames <- gsub(modelsDirectory,"",layersToCalc)
layersToCalcNames <- gsub("/Rasters","",layersToCalcNames)
layersToCalcNames <- gsub("/","",layersToCalcNames)
layersToCalcNames <- gsub(".tif","",layersToCalcNames)
layersToCalcNames <- gsub("Stacking","",layersToCalcNames)
layersToCalcNames

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
## Global summary

resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*3)+1,nrow=1))
colnames(resultsDF) <- c("Region",paste0(layersToCalcNames,"MaxNumber"),paste0(layersToCalcNames,"MeanNumber"),paste0(layersToCalcNames,"Area"))
resultsDF[,1] <- "Global"

## ----------------

areaLayer <- raster::area(raster(layersToCalc[1]))

for(i in 1:length(layersToCalcNames)) {
  
  rasterLayer <- raster(layersToCalc[i]) 
  rasterLayerBinomial <- rasterLayer
  rasterLayerBinomial[rasterLayerBinomial > 1] <- 1
  rasterLayerArea <- rasterLayerBinomial * areaLayer

  rasterLayerVals <- getValues(rasterLayer)
  rasterLayerValsArea <- getValues(rasterLayerArea)
  
  resultsDF[1,i+1] <- max( unlist(rasterLayerVals),na.rm=T)
  resultsDF[1,i+1+length(layersToCalcNames)] <- mean( unlist(rasterLayerVals),na.rm=T)
  resultsDF[1,i+1+(length(layersToCalcNames)*2)] <- sum( unlist(rasterLayerValsArea) ,na.rm=T )
  
}

# Remove nonsense results and sort
resultsDF <- resultsDF[,c(-2,-3,-10,-11,-15,-16,-37,-36)]
resultsDF <- resultsDF[,c(1,sort(colnames(resultsDF[,-1]),index.return=T)$ix+1)]

write.csv(resultsDF, file=paste0(resultsFolder,"/summaryGlobal.csv"), row.names = FALSE)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
## Per Polygons

resultsName <- "MarineProvinces"

polygon <- "/Volumes/Jellyfish/Dropbox/Data/Spatial information/Shapefiles/Marine ecoregions/marine_ecoregions.shp"
polygonFeature <- "REALM" # NULL ECOREGION PROVINCE REALM

# polygon <- drawPolygon(mainLayer)

## ----------------

if(class(polygon) == "character") { polygon <- shapefile(polygon) }

if( ! is.null(polygonFeature) ) { polygonFeatureNames <- polygon@data[which(names(polygon@data) == polygonFeature)][,1] }

resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*3)+1,nrow=length(unique(polygonFeatureNames))))
colnames(resultsDF) <- c(resultsName,paste0(layersToCalcNames,"MaxNumber"),paste0(layersToCalcNames,"MeanNumber"),paste0(layersToCalcNames,"Area"))
resultsDF[,1] <- unique(polygonFeatureNames)

## ----------------

areaLayer <- raster::area(raster(layersToCalc[1]))

for(i in 1:length(layersToCalcNames)) {
  
  rasterLayer <- raster(layersToCalc[i]) 
  rasterLayerBinomial <- rasterLayer
  rasterLayerBinomial[rasterLayerBinomial > 1] <- 1
  rasterLayerArea <- rasterLayerBinomial * areaLayer
  
  for( p.i in 1:length(unique(polygonFeatureNames)) ){
    
    polygon.i <- polygon[which(polygon@data[which(names(polygon@data) == polygonFeature)][,1] == unique(polygonFeatureNames)[p.i]),]
    rasterLayerVals <- extract(rasterLayer,polygon.i)
    rasterLayerValsArea <- extract(rasterLayerArea,polygon.i)
    
    if( ! is.null(rasterLayerVals) ) { 
      
      resultsDF[p.i,i+1] <- max( unlist(rasterLayerVals),na.rm=T)
      resultsDF[p.i,i+1+length(layersToCalcNames)] <- mean( unlist(rasterLayerVals),na.rm=T)
      resultsDF[p.i,i+1+(length(layersToCalcNames)*2)] <- sum( unlist(rasterLayerValsArea) ,na.rm=T )
      
    }
  }
}

resultsDF[resultsDF == -Inf] <- NA

## ----------------

polygonExport <- polygon

if( nrow(resultsDF) == length(polygonExport) ) {
  
  row.names(resultsDF) <- resultsDF[,1]
  polygonExport <- SpatialPolygonsDataFrame(polygonExport, resultsDF)
  
}

if( nrow(resultsDF) != length(polygonExport) ) {
  
  polygonExport <- gUnaryUnion(polygonExport, id = polygonExport@data[which(names(polygonExport@data) == polygonFeature)][,1] )
  row.names(resultsDF) <- resultsDF[,1]
  polygonExport <- SpatialPolygonsDataFrame(polygonExport, resultsDF)
  
}

row.names(resultsDF) <- NULL
write.csv(resultsDF, file=paste0(resultsFolder,"/summaryPerPolygon",resultsName,".csv"), row.names = FALSE)
shapefile(polygonExport, filename=paste0(resultsFolder,"/summaryPerPolygon",resultsName,".shp"), overwrite=TRUE)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
## Per latitude

areaLayer <- raster::area(raster(layersToCalc[1]))
regionOfInterest <- extent(raster(layersToCalc[1]))
perLatitude <- data.frame()
bins <- seq(extent(regionOfInterest)[3],extent(regionOfInterest)[4], by=1)
bins <- data.frame(binFrom=bins[-length(bins)],binTo=bins[-1])

resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*3),nrow=nrow(bins)))
colnames(resultsDF) <- c(paste0(layersToCalcNames,"Max"),paste0(layersToCalcNames,"Mean"),paste0(layersToCalcNames,"Area"))
resultsDF <- data.frame(bins,resultsDF)

for(i in 1:length(layersToCalcNames)) {
  
  rasterLayer <- raster(layersToCalc[i]) 
  rasterLayerBinomial <- rasterLayer
  rasterLayerBinomial[rasterLayerBinomial > 1] <- 1
  rasterLayerArea <- rasterLayerBinomial * areaLayer
  
  for(b in 1:nrow(resultsDF)) {

    rasterLayer.bin <- crop(rasterLayer,extent(extent(regionOfInterest)[1],extent(regionOfInterest)[2],bins[b,1],bins[b,2]))
    rasterLayerArea.bin <- crop(rasterLayerArea,extent(extent(regionOfInterest)[1],extent(regionOfInterest)[2],bins[b,1],bins[b,2]))
    
    resultsDF[b,i+2] <- max( getValues(rasterLayer.bin) ,na.rm=T)
    resultsDF[b,i+2+length(layersToCalcNames)] <- mean( getValues(rasterLayer.bin) ,na.rm=T)
    resultsDF[b,i+2+(length(layersToCalcNames)*2)] <- sum( getValues(rasterLayerArea.bin) ,na.rm=T)

} }

row.names(resultsDF) <- NULL
write.csv(resultsDF, file=paste0(resultsFolder,"/summaryPerLatitudinalBin.csv"), row.names = FALSE)

fig1 <- ggplot() +
  geom_line( data=resultsDF,aes(x=binFrom, y=SpeciesRichnessPresentArea) , color="gray", size=0.4) +
  theme_minimal() +
  xlab("Latitude") + 
  ylab("SpeciesRichnessPresentArea") + geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.7) + coord_flip()

pdf(paste0(dumpFolder,"/ ",predictors[pred]," ",resultsName," Latitude.pdf"), width = 10 )
print(fig1)
dev.off()


