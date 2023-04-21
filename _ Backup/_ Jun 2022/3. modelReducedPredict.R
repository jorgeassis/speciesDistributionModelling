## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##  
##  Centre of Marine Sciences [ ccmar.ualg.pt ]
##  Faro, Portugal
##
## ---------------------------------------------------------------------
##
##  Machine Learning Species Distribution Modelling [ Ver.301 ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("Dependencies/mainFunctions.R")
source("0. Config.R")
load(file=paste0(resultsDirectory,"/RData/","envFullModel.RData"))

# -----------------------------

rasterLayers <- listAllFiles(dataLayersDirectory,dataLayersFileType)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
names(rasterLayers) <- dataLayersName
rasterLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer,minDepth-minDepthBuffer,maxDepth+maxDepthBuffer,intertidalRegion)

# -----------------------------
# Drop no Variation Predictors

cave <- function(x) { length(unique(x)) / length(x) }
randomLocations <- Which(!is.na(subset(rasterLayers,1)),cells=TRUE)
randomLocations <- xyFromCell( subset(rasterLayers,1) , sample(randomLocations,min(length(randomLocations),1000),replace=FALSE))
varRasterLayers <- which( apply( extract(rasterLayers,randomLocations) ,2,cave) > 0.025)
rasterLayers <- subset(rasterLayers,varRasterLayers)

# -----------------------------
# Correct layers

# rasterLayers <- correctLayer(rasterLayers,"Salinity","he",28,28)

# -----------------------------
# Correlated Pairs

pairsThershold <- 0.9
pairsCorrelated <- correlatedPairs(rasterLayers,speciesData,pairsThershold,dataLayersMonotonocity)

allComb <- paste0("expand.grid(",paste0(sapply(names(rasterLayers),function(r){ paste0(r,"=c(0,1)") }),collapse = ", "),",stringsAsFactors = TRUE)")
allComb <- eval(parse(text=allComb))

if(nrow(pairsCorrelated) > 0) {
  for(i in 1:nrow(pairsCorrelated)) {
    allComb <- allComb[-which(allComb[,pairsCorrelated[i,1]] == 1 & allComb[,pairsCorrelated[i,2]] == 1),]
  }
}

allComb <- allComb[apply(allComb,1,sum) >= ncol(allComb) - nrow(pairsCorrelated) ,]
allComb

# ------------------------------------------------------------------------------------
# Produce models

ensembleBRT <- list()
ensembleMBOOST <- list()
ensembleMPNN <- list()
varsKeptBRT <- character(0)
varsKeptMBOOST <- character(0)
keepPredictorsThreshold <- 0

for( algorithmType in algorithms ) {
  
  modelType <- algorithmType
  
  for( m in 1:nrow(allComb)) {
    
    rasterLayersBkSimp <- rasterLayers
    rasterLayers <- subset(rasterLayers,which(names(rasterLayers) %in% colnames(allComb)[which(allComb[m,] == 1)]))
    
    if( nrow(allComb) == 1) { fullModel <- modelBRTFull }
    
    if( nrow(allComb) > 1 ) { 
      
      source("Dependencies/modelMe.R")
      fullModel <- modelFull
      
    }
    
    simplifyType <- "auc" # deviance
    relativeContribution <- relativeContributionBRT
    source("Dependencies/simplifyModel.R")
    
    pdf(file = paste0(resultsDirectory,"/simplifyModel",modelType,m,".pdf"), onefile=FALSE, width=8 )
    print(model.m$plot)
    dev.off()
    
    if( modelType == "BRT") { 
      ensembleBRT <- c(ensembleBRT,list(model.m$model))
      varsKeptBRT <- unique(c(varsKeptBRT,model.m$keptVariables))
    }
    
    if( modelType == "MBOOST") { 
      ensembleMBOOST <- c(ensembleMBOOST,list(model.m$model))
      varsKeptMBOOST <- unique(c(varsKeptMBOOST,model.m$keptVariables))
    }
    
    rasterLayers <- rasterLayersBkSimp
    
  }
}

# ------------------------------------------------------------------------------------
# Save models

paste0( varsKeptBRT ,collapse = ",")
save(ensembleBRT,file=paste0( resultsDirectory, "/RData/ensembleBRT.RData"))

paste0( varsKeptMBOOST ,collapse = ",")
save(ensembleMBOOST,file=paste0( resultsDirectory, "/RData/ensembleMBOOST.RData"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

dataLayersDirectory <- "../Data/Climate/LGM/"
predictName <- "LGM" # Present MH LGM LIG RCP26 RCP85

newIntertidalRegion <- intertidalRegion
coastLineDataLayer <- "Dependencies/Data/Rasters/BO2CoastLineLGM.tif" # BO2CoastLineLGM.tif BO2CoastLine.tif
regionBuffer <- 20
newMinDepthBuffer <- minDepthBuffer
newMaxDepthBuffer <- maxDepthBuffer

rasterLayers <- listAllFiles(dataLayersDirectory,dataLayersFileType)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
names(rasterLayers) <- dataLayersName
predictLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer,ifelse(minDepth-newMinDepthBuffer < 0 , 0 , minDepth-newMinDepthBuffer),maxDepth+newMaxDepthBuffer,newIntertidalRegion)
plot(subset(predictLayers,1))

# ------------------------------------

finalEnsemble <- list()
for( algorithmType in algorithms) {
  if(algorithmType == "BRT") { finalEnsemble <- c(finalEnsemble,list( sapply(1:length(ensembleBRT), function(x) { predictDistribution(predictLayers,ensembleBRT[[x]] , reclassToOne=FALSE) } ))) }
  if(algorithmType == "MBOOST") { finalEnsemble <- c(finalEnsemble,list( sapply(1:length(ensembleMBOOST), function(x) { predictDistribution(predictLayers,ensembleMBOOST[[x]] , reclassToOne=FALSE) } ))) }
  if(algorithmType == "MPNN") { stop("!") }
}

beginCluster(nCores)
finalEnsembleSD <- clusterR(stack(unlist(finalEnsemble)), calc, args=list(sd))
finalEnsemble <- clusterR(stack(unlist(finalEnsemble)), calc, args=list(median))
endCluster()

plot(finalEnsemble)

# --------------------

accuracyPredicted(finalEnsemble,speciesData,"area") # Use threshold dependent index
write.csv(t(accuracyPredicted(finalEnsemble,speciesData,"tss") ),paste0(resultsDirectory,"/summaryModelEnsemblePerformance.csv"), row.names = TRUE)

# -------------------

reclassThreshold <- 0.1414141 #
finalEnsembleReclass <- reclassifyPredicted(finalEnsemble,speciesData,method="directReclass",reclassThreshold)
plot(finalEnsembleReclass)

writeRaster(finalEnsemble,filename=paste0(resultsDirectory,"/Ensemble.",predictName,".tif"),format="GTiff",overwrite=T)
writeRaster(finalEnsembleReclass,filename=paste0(resultsDirectory,"/EnsembleReclass.",predictName,".tif"),format="GTiff",overwrite=T)

# -------------------
# -------------------

# Reachability

finalEnsembleReclass[finalEnsembleReclass == 0] <- NA
finalEnsembleReclassRegions <- raster::clump(finalEnsembleReclass)

regionsUsed <- extract(finalEnsembleReclassRegions,occurrenceRecords)
regionsUsed <- regionsUsed[!is.na(regionsUsed)]

finalEnsembleReclass <- finalEnsembleReclassRegions %in% regionsUsed
finalEnsembleReclass[finalEnsembleReclass == 0] <- NA
plot(finalEnsembleReclass, col="Black")
writeRaster(finalEnsembleReclass,filename=paste0(resultsDirectory,"/EnsembleReclassReachable_",predictName,".tif"),format="GTiff",overwrite=T)

# ---------------------
# Save Environment

save.image(file=paste0(resultsDirectory,"/RData/envReducedPredict.RData"))
# load(paste0(resultsDirectory,"envReducedPredict.RData"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# Ensembe tipping points

finalEnsembleReclassi <- finalEnsembleReclass
finalEnsembleReclassi[finalEnsembleReclassi >= reclassThreshold] <- 1
finalEnsembleReclassi[finalEnsembleReclassi < reclassThreshold] <- NA
ensembleTP <- data.frame()

for(i in 1:length(names(predictLayers))) {
  
  temp <- getValues(subset(predictLayers,i) * finalEnsembleReclassi)
  temp <- range(temp,na.rm=T)
  ensembleTP <- rbind(ensembleTP,data.frame(predicor=names(subset(predictLayers,i)) , as.character(temp)))
  
}
ensembleTP
write.csv(ensembleTP,paste0(resultsDirectory,"/summaryModelEnsembleTippingPoints.csv"), row.names = TRUE)

# End of Code
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
