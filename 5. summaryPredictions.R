## -------------------------------------------------------------------------------
## --------------------------------------------------
## --------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("Dependencies/mainFunctions.R")

nCores <- 12

tempFolder <- "/Volumes/Jellyfish/Temp"

modelsDirectory <- "/Volumes/Jellyfish/Dropbox/Data/Distribution Models/Global distribution of marine forests/Results [Full models]/_ Per Species/"
resultsFolder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Projected climate changes threaten marine forests biodiversity at global scales/Results/"

## ----------------

resultsFolder <- paste0(resultsFolder,"/Summary/")
if(! dir.exists(paste0(resultsFolder))) { dir.create(paste0(resultsFolder), recursive = T) }

## ----------------
## Get layers

speciesToCalc <- list.files(modelsDirectory, full.names = F)
speciesToCalc <- speciesToCalc[sapply(1:length(speciesToCalc) ,  function(x) { length(list.files( paste0(modelsDirectory,speciesToCalc[x]) , full.names = T, recursive = T, pattern = "Ensemble_Present.tif")) } ) > 0]

layersToCalc <- list.files( paste0(modelsDirectory,speciesToCalc[1]) , full.names = T, recursive = T, pattern = ".tif")
layersToCalc <- layersToCalc[grepl("/Maps/Global/",layersToCalc)]
layersToCalc <- layersToCalc[!grepl("/_ Backup/",layersToCalc)]
layersToCalc <- layersToCalc[!grepl("Reachable",layersToCalc)]
layersToCalc <- layersToCalc[grepl("Reclass",layersToCalc)]
layersToCalc <- layersToCalc[!grepl("GainLoss",layersToCalc)]

layersToCalcNames <- sapply(1:length(layersToCalc) ,  function(x) { substr(layersToCalc[x], unlist(gregexpr("/",layersToCalc[x]))[length(unlist(gregexpr("/",layersToCalc[x])))] + 1 , nchar(layersToCalc[x])) } )
layersToCalcNames <- gsub(".tif","",layersToCalcNames)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
## Global summary for numerous species

resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames))+1,nrow=length(speciesToCalc)))
colnames(resultsDF) <- c("Species",paste0(layersToCalcNames,"Area"))
resultsDF[,1] <- speciesToCalc

resultsDF$LatChangeRCP26 <- NA
resultsDF$LatChangeRCP85 <- NA

## ----------------

areaLayer <- raster::area(raster(layersToCalc[1]))

for( sp in speciesToCalc ) {
  
    sp.loc <- which( speciesToCalc == sp)
    
    layersToCalcSP <- list.files( paste0(modelsDirectory,sp) , full.names = T, recursive = T, pattern = ".tif")
    layersToCalcSP <- layersToCalcSP[unlist(sapply( layersToCalcNames , function(x) { which(grepl(paste0("Global/",x),layersToCalcSP)) }))]
    
    for(i in 1:length(layersToCalcSP)) {
      
      rasterLayer <- raster(layersToCalcSP[i]) 
      rasterLayerBinomial <- rasterLayer
      rasterLayerBinomial[rasterLayerBinomial > 1] <- 1
      rasterLayerArea <- rasterLayerBinomial * areaLayer
      rasterLayerValsArea <- getValues(rasterLayerArea)
      resultsDF[sp.loc,i+1] <- sum( unlist(rasterLayerValsArea) ,na.rm=T )
      
    }
  
    # Latitudinal change in km
  
    rasterLayerPresent <- raster(layersToCalcSP[1]) 
    rasterLayerPresent[rasterLayerPresent == 0] <- NA
    rasterLayerPresentLat <- xyFromCell(rasterLayerPresent,Which(rasterLayerPresent == 1, cell=T))[,2]
    rasterLayerPresentLat <- (min(abs(rasterLayerPresentLat))) + ( (max(abs(rasterLayerPresentLat)) - min(abs(rasterLayerPresentLat)) ) / 2)
    
    rasterLayer <- raster(layersToCalcSP[2])
    rasterLayer[rasterLayer == 0] <- NA
    rasterLayerLat <- xyFromCell(rasterLayer,Which(rasterLayer == 1, cell=T))[,2]
    rasterLayerLat <- (min(abs(rasterLayerLat))) + ( (max(abs(rasterLayerLat)) - min(abs(rasterLayerLat)) ) / 2)
    resultsDF[sp.loc,"LatChangeRCP26"] <- ( max(c(rasterLayerPresentLat,rasterLayerLat)) - min(c(rasterLayerPresentLat,rasterLayerLat)) ) * 110.574 # km
    
    rasterLayer <- raster(layersToCalcSP[3])
    rasterLayer[rasterLayer == 0] <- NA
    rasterLayerLat <- xyFromCell(rasterLayer,Which(rasterLayer == 1, cell=T))[,2]
    rasterLayerLat <- (min(abs(rasterLayerLat))) + ( (max(abs(rasterLayerLat)) - min(abs(rasterLayerLat)) ) / 2)
    resultsDF[sp.loc,"LatChangeRCP85"] <- ( max(c(rasterLayerPresentLat,rasterLayerLat)) - min(c(rasterLayerPresentLat,rasterLayerLat)) ) * 110.574 # km

}

write.csv(resultsDF, file=paste0(resultsFolder,"/summaryPerSpecies.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Summary of predictions

# Review 

areaPredicted <- raster::area(finalEnsembleReclassM) * finalEnsembleReclassM
sum(getValues(areaPredicted),na.rm=TRUE) # km2

bathyPredicted <- finalEnsembleReclassM
bathyPredicted[bathyPredicted == 0] <- NA
bathyPredicted <- crop(raster(bathymetryDataLayer),bathyPredicted) * bathyPredicted * (-1)

bathyPredicted <- getValues(bathyPredicted)
bathyPredicted<- bathyPredicted[!is.na(bathyPredicted)] 
bathyPredicted<- bathyPredicted[bathyPredicted != 0] 

min(bathyPredicted)
mean(bathyPredicted)
sd(bathyPredicted)
max(bathyPredicted)
quantile(bathyPredicted,probs=0.95)
quantile(bathyPredicted,probs=0.025)

bathyPredictedHist <- hist(bathyPredicted, breaks=100)
bathyPredictedHist <- data.frame(depth=bathyPredictedHist$breaks[-length(bathyPredictedHist$breaks)],Frequency=bathyPredictedHist$counts)
bathyPredictedHist[,2] <- bathyPredictedHist[,2] / sum(bathyPredictedHist[,2])

bathyPredictedPlot <- ggplot(bathyPredictedHist) + ylim(0, 0.04) +
  geom_bar( aes(x=depth, y=Frequency), stat="identity", fill="black", alpha=0.5) +
  theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
    panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF") ) + 
  labs(x = "Depth (m)") + labs(y = "Suitable Habitats (proportion)") + 
  geom_vline(xintercept = quantile(bathyPredicted,probs=0.95), linetype="dashed", color = "#575757", size=0.7) +
  # geom_vline(xintercept = quantile(bathyPredicted,probs=0.5), linetype="dashed", color = "#575757", size=0.7) +
  geom_vline(xintercept = mean(bathyPredicted), linetype="solid", color = "black", size=0.7)

bathyPredictedPlot

pdf(file = paste0(resultsDirectory,"/depthDistribution.",predictName,".pdf"), onefile=FALSE, width=10, height=7 )
bathyPredictedPlot
dev.off()

# ----------------------------
# Depth shifts

predictName <- "RCP85" # Present MH LGM LIG RCP26 RCP85
finalEnsembleReclassM <- raster(paste0(resultsDirectory,"/lossGain.",predictName,".tif"))

# Loss 1

loss <- finalEnsembleReclassM
loss[loss != 1 ] <- NA
sum(getValues( raster::area(loss) * loss),na.rm=TRUE) # km2

# Gain 2
gain <- finalEnsembleReclassM
gain[gain != 2 ] <- NA
sum(getValues( raster::area(gain) * gain),na.rm=TRUE) # km2

# -------

bathyPredicted <- gain
bathyPredicted[bathyPredicted == 0] <- NA
bathyPredicted <- crop(raster(bathymetryDataLayer),bathyPredicted) * bathyPredicted * (-1)

bathyPredicted <- getValues(bathyPredicted)
bathyPredicted<- bathyPredicted[!is.na(bathyPredicted)] 
bathyPredicted<- bathyPredicted[bathyPredicted != 0] 

mean(bathyPredicted)
sd(bathyPredicted)
max(bathyPredicted)
quantile(bathyPredicted,probs=0.95)

bathyPredictedHist <- hist(bathyPredicted, breaks=100)
bathyPredictedHist <- data.frame(depth=bathyPredictedHist$breaks[-length(bathyPredictedHist$breaks)],Frequency=bathyPredictedHist$counts)
bathyPredictedHist[,2] <- bathyPredictedHist[,2] / sum(bathyPredictedHist[,2])

bathyPredictedPlot <- ggplot(bathyPredictedHist) +
  geom_bar( aes(x=depth, y=Frequency), stat="identity", fill="black", alpha=0.5) +
  theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
    panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF") ) + 
  labs(x = "Depth (m)") + labs(y = "Suitable Habitats (proportion)") + 
  geom_vline(xintercept = quantile(bathyPredicted,probs=0.95), linetype="dashed", color = "#575757", size=0.7) +
  # geom_vline(xintercept = quantile(bathyPredicted,probs=0.5), linetype="dashed", color = "#575757", size=0.7) +
  geom_vline(xintercept = mean(bathyPredicted), linetype="solid", color = "black", size=0.7)

bathyPredictedPlot

pdf(file = paste0(resultsDirectory,"/depthDistribution.Gain.",predictName,".pdf"), onefile=FALSE, width=10, height=7 )
bathyPredictedPlot
dev.off()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# End of Code