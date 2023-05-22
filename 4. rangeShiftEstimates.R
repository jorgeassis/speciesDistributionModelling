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

# -----------------------------
# Read data

modelDataList <- list.files(paste0(resultsDirectory,"/Predictions/"), pattern="RData", full.names = TRUE, recursive = TRUE)

if(   "Reclass" %in% exportType ) { modelDataList.i <- modelDataList[grepl("Reclass",modelDataList)] }
if( ! "Reclass" %in% exportType ) { modelDataList.i <- modelDataList[!grepl("Reclass",modelDataList)] }

if(   "Reachable" %in% exportType ) { modelDataList.i <- modelDataList.i[grepl("Reachable",modelDataList.i)] }
if( ! "Reachable" %in% exportType ) { modelDataList.i <- modelDataList.i[!grepl("Reachable",modelDataList.i)] }

if(   "Global" %in% exportType ) { modelDataList.i <- modelDataList.i[grepl("Global.RData",modelDataList.i)] }    
if( ! "Global" %in% exportType ) { modelDataList.i <- modelDataList.i[!grepl("Global.RData",modelDataList.i)] }    

baselinePredictionFile <- modelDataList.i[grepl("Baseline",modelDataList.i)]
baselinePrediction <- loadRData(baselinePredictionFile)

# ----------

shape <- baselinePrediction
shape[shape == 1] <- 0

for( scenario in scenariosToPredict[scenariosToPredict != "Baseline"] ) {

  scenarioPredictionFile <- gsub("Baseline",scenario,baselinePredictionFile)
  scenarioPrediction <- loadRData(scenarioPredictionFile)

  rangeShifts <- Gain <- Loss <- Refugia <- shape
  
  if( ! 1 %in% getValues(scenarioPrediction) ) { 
    
    Loss <- baselinePrediction
    rangeShifts <- baselinePrediction
    rangeShifts[rangeShifts == 1] <- -1
    
    save( Gain , file = gsub("ensemble","ensembleGain",scenarioPredictionFile), compress=TRUE, compression_level=6)
    save( Loss , file = gsub("ensemble","ensembleLoss",scenarioPredictionFile), compress=TRUE, compression_level=6)
    save( Refugia , file = gsub("ensemble","ensembleRefugia",scenarioPredictionFile) , compress=TRUE, compression_level=6)
    save( rangeShifts , file = gsub("ensemble","ensembleRangeShifts",scenarioPredictionFile) , compress=TRUE, compression_level=6)
    
    save( shape , file=gsub("ensemble","ensembleGainDistance",scenarioPredictionFile) , compress=TRUE, compression_level=6)
    save( shape , file=gsub("ensemble","ensembleRefugiaDistance",scenarioPredictionFile) , compress=TRUE, compression_level=6)
    
    rangeShifts <- data.frame(LatShift = 0, distanceCentroidShift=0, distanceToRefugiaMean=0, distanceToRefugiaMax=0, distanceRangeGainMean=0, distanceRangeGainMax=0)
    save( rangeShifts , file=paste0( resultsDirectory, "/SummaryModels/rangeShiftEstimates",ifelse( "Reachable" %in% exportType , "Reachable" , "unConstrained" ),scenario,".RData"))
    
    next
    
  }

  rangeShifts[Which(baselinePrediction == 0 , cells=T )] <- NA
  rangeShifts[Which(scenarioPrediction == 1 & baselinePrediction == 0 , cells=T )] <- 1
  rangeShifts[Which(scenarioPrediction == 0 & baselinePrediction == 1 , cells=T )] <- -1
  rangeShifts[Which(scenarioPrediction == 1 & baselinePrediction == 1 , cells=T )] <- 0
  
  Gain[Which(scenarioPrediction == 1 & baselinePrediction == 0 , cells=T )] <- 1
  Loss[Which(scenarioPrediction == 0 & baselinePrediction == 1 , cells=T )] <- 1
  Refugia[Which(scenarioPrediction == 1 & baselinePrediction == 1 , cells=T )] <- 1

  save( Gain , file = gsub("ensemble","ensembleGain",scenarioPredictionFile), compress=TRUE, compression_level=6)
  save( Loss , file = gsub("ensemble","ensembleLoss",scenarioPredictionFile), compress=TRUE, compression_level=6)
  save( Refugia , file = gsub("ensemble","ensembleRefugia",scenarioPredictionFile) , compress=TRUE, compression_level=6)
  save( rangeShifts , file = gsub("ensemble","ensembleRangeShifts",scenarioPredictionFile) , compress=TRUE, compression_level=6)
  
  # -------------
  # Range shift estimates
  
  rangeShifts <- data.frame(LatShift = abs(mean(range(xyFromCell(Gain,Which(Gain == 1, cells=T))[,2]))) - abs(mean(range(xyFromCell(Refugia,Which(Refugia == 1, cells=T))[,2]))))

  rangeShifts$distanceCentroidShift <- spDists( matrix( c(mean((range(xyFromCell(baselinePrediction,Which(baselinePrediction == 1, cells=T))[,1]))),mean((range(xyFromCell(baselinePrediction,Which(baselinePrediction == 1, cells=T))[,2])))), nrow=1),
                                                matrix( c(mean((range(xyFromCell(scenarioPrediction,Which(scenarioPrediction == 1, cells=T))[,1]))),mean((range(xyFromCell(scenarioPrediction,Which(scenarioPrediction == 1, cells=T))[,2])))), nrow=1),
                                                longlat = TRUE )

  
  distRefugia <- distExpansion <- shape
  distRefugia[] <- NA
  distExpansion[] <- NA
  
  # -------------
  # Distance to refugia 
  
  cellsbaselinePrediction <- Which(baselinePrediction == 1 , cells=TRUE)
  cellsRefugiaPrediction <- Which(Refugia == 1 , cells=TRUE)
  coordsRefugiaPrediction <- xyFromCell(distRefugia,cellsRefugiaPrediction)
  
  cellsZeroDistance <- intersect(cellsbaselinePrediction,cellsRefugiaPrediction)
  distRefugia[cellsZeroDistance] <- 0
  
  cellsNonZeroDistance <- setdiff(cellsbaselinePrediction,cellsRefugiaPrediction)
  cellsNonZeroDistanceCalc <- numeric(length(cellsNonZeroDistance))
  coordsRangeshifts <- xyFromCell(distRefugia,cellsNonZeroDistance)
  
  if( nrow(coordsRefugiaPrediction) > 0 & nrow(coordsRangeshifts) > 0) {
    distanceCalc <- get.knnx( coordsRefugiaPrediction , matrix(coordsRangeshifts, ncol=2), k= 1 , algorithm="kd_tree" )
    distanceCalc <- sapply(1:nrow(coordsRangeshifts), function(x) { d2km(distanceCalc$nn.dist[x], base.latitude = coordsRangeshifts[x,2]) })
    distRefugia[cellsNonZeroDistance] <- distanceCalc
  }
  
  save( distRefugia , file=gsub("ensemble","ensembleRefugiaDistance",scenarioPredictionFile) , compress=TRUE, compression_level=6)
  
  distRefugia[distRefugia == 0 ] <- NA
  rangeShifts$distanceToRefugiaMean <- cellStats(distRefugia,mean)
  rangeShifts$distanceToRefugiaMax <- cellStats(distRefugia,max)
  
  # Expansion region
  
  cellsbaselinePrediction <- Which(baselinePrediction == 1 , cells=TRUE)
  cellsGainPrediction <- Which(Gain == 1 , cells=TRUE)
  
  coordsBaselinePrediction <- xyFromCell(distExpansion,cellsbaselinePrediction)
  coordsGainPrediction <- xyFromCell(distExpansion,cellsGainPrediction)
  
  if( nrow(coordsBaselinePrediction) > 0 & nrow(coordsGainPrediction) > 0) {

    distanceCalc <- get.knnx( coordsBaselinePrediction , matrix(coordsGainPrediction, ncol=2), k= 1 , algorithm="kd_tree" )
    distanceCalc <- sapply(1:nrow(coordsGainPrediction), function(x) { d2km(distanceCalc$nn.dist[x], base.latitude = coordsGainPrediction[x,2]) })
    distExpansion[cellsGainPrediction] <- distanceCalc
  }
  
  save( distExpansion , file=gsub("ensemble","ensembleGainDistance",scenarioPredictionFile) , compress=TRUE, compression_level=6)
  
  distExpansion[distExpansion == 0 ] <- NA
  rangeShifts$distanceRangeGainMean <- cellStats(distExpansion,mean)
  rangeShifts$distanceRangeGainMax <- cellStats(distExpansion,max)
  
  save( rangeShifts , file=paste0( resultsDirectory, "/SummaryModels/rangeShiftEstimates",ifelse( "Reachable" %in% exportType , "Reachable" , "unConstrained" ),scenario,".RData"))
  
}
