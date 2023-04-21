## -----------------------

# Remove variables with permutation importance lower than x% only if testing
# TSS doesn't decrease

rasterLayersBk <- rasterLayers
modelDataBk <- modelData

keptVariables <- names(subset(rasterLayers, which(names(rasterLayers) %in% relativeContributionDF[which(relativeContributionDF$relImportance > reduceModelcomplexityThreshold),"Predictor"]) ))
if(length(keptVariables) == 1) { keptVariables <- as.character(relativeContributionDF[ sort(relativeContributionDF$relImportance, index.return=T, decreasing=T)$ix,1])[1:2] }
removedVariables <- names(rasterLayers)[which(! names(rasterLayers) %in% keptVariables)]

if( "layer" %in% keptVariables) { stop("Error :: 914") }

if( length(removedVariables) > 0) {
  
  rasterLayers <- subset(rasterLayers,which(names(rasterLayers) %in% keptVariables))
  modelData <- prepareModelData(occurrenceRecords, pseudoAbsences, rasterLayers)
  
  source("Dependencies/modelMe.R", local = TRUE)
  
}

rasterLayers <- rasterLayersBk
modelData <- modelDataBk

## -----------------------

modelReduced <- list(model=model,keptVariables=keptVariables,removedVariables=removedVariables)

## -----------------------

gc(reset=TRUE)
