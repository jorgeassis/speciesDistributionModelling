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

## ----------------------
## Relative contribution of predictors

algorithmsSP <- (c(algorithms,"ensemble"))[ sapply(c(algorithms,"ensemble"), function(x) { length(list.files(paste0(dataDirectory,"/SummaryModels/"), pattern=paste0(x,"Contrib"))) > 0 } )]
contributionDF <- data.frame(name=species)

for( algorithmType in algorithmsSP ) {
  contributionDF.i <- list.files(paste0(dataDirectory,"/SummaryModels/"), pattern=paste0(algorithmType,"Contrib"), full.names = TRUE)[1]
  contributionDF.i <- loadRData(contributionDF.i)
  contributionDF.c <- data.frame(matrix(NA,ncol=length(dataLayersName),nrow=1))
  colnames(contributionDF.c) <- paste0(dataLayersName,".",algorithmType)

  for( pred in contributionDF.i$Predictor ) {
    colMatch <- sapply( paste0(pred,".",algorithmType) , function(x) which(colnames(contributionDF.c) == x ) )
    if(length(colMatch[[1]]) == 0) { next }
    contributionDF.c[1, colMatch ] <- contributionDF.i[which(contributionDF.i$Predictor == pred),"relImportance"]
  }
  contributionDF.c[is.na(contributionDF.c)] <- 0
  contributionDF <- cbind(contributionDF,contributionDF.c)
}

## ----------------------
## Performance

algorithmsSP <- (c(algorithms,"ensemble"))
performanceDF <- data.frame(name=species)

for( algorithmType in algorithmsSP ) {
  
  if( algorithmType != "ensemble" ) {
    
    performance.i <- list.files(paste0(dataDirectory,"/Models/"), pattern=paste0(algorithmType), full.names = TRUE)
    performance.i <- performance.i[which(!grepl("reduced",performance.i))]
    performance.i <- loadRData(performance.i)
    performanceCV <- performance.i$performanceCV[-1]
    performanceCV.sd <- performance.i$performanceCV.sd[-1]
    performance <- performance.i$performance[-1]

    performance.i <- c(interleave(performanceCV.sd,performanceCV),unlist(performance))
    names(performance.i) <- c(interleave(paste0("sdCV.",names(performanceCV),".",algorithmType),paste0("meanCV.",names(performanceCV.sd),".",algorithmType)),paste0(names(performance),".",algorithmType))
    
    performanceDF <- cbind(performanceDF,t(data.frame(performance.i)))
    
  }
  
  if( algorithmType == "ensemble" ) {
    
    performance.i <- list.files(paste0(dataDirectory,"/SummaryModels/"), pattern=paste0(algorithmType,"Performance.RData"), full.names = TRUE)
    performance.i <- loadRData(paste0(dataDirectory,"/SummaryModels/",algorithmType,"Performance.RData"))
    performance.i <- performance.i[-1]
    colnames(performance.i) <- paste0("",names(performance.i),".","Ensemble")
    performanceDF <- cbind(performanceDF,performance.i)
    
  }

}
row.names(performanceDF) <- NULL

## ----

if( "ensemblePerformanceReachable.RData" %in% list.files(paste0(dataDirectory,"/SummaryModels/") ) ) {
  performance.r <- loadRData(paste0(dataDirectory,"/SummaryModels/ensemblePerformanceReachable.RData"))[,-1]
  names(performance.r) <- paste0(names(performance.r),".EnsembleReachable")
  performanceDF <- cbind(performanceDF,performance.r)
}

## ----------------------
## Tipping Points

algorithmsSP <- (c(algorithms,"Ensemble"))[ sapply(c(algorithms,"Ensemble"), function(x) { length(list.files(paste0(dataDirectory,"/SummaryModels/"), pattern=paste0(x,"TippingPoints"))) > 0 } )]
tippingPointsDF <- data.frame(name=species)

for( algorithmType in algorithmsSP ) {
  
  tippingPoints.i <- loadRData(paste0(dataDirectory,"/SummaryModels/",algorithmType,"TippingPoints.RData"))
  tippingPoints.c <- data.frame(matrix(NA,ncol=length(dataLayersName)*(ncol(tippingPoints.i)-1),nrow=1))
  colnames(tippingPoints.c) <- c(sapply(colnames(tippingPoints.i)[-1],function(x) paste0(dataLayersName,".",algorithmType,".",x) ))
    
  for( c in 2:ncol(tippingPoints.i)) {
    colMatch <- sapply( paste0(tippingPoints.i[,1],".",algorithmType,".",colnames(tippingPoints.i)[c]), function(x) which(colnames(tippingPoints.c) == x ) )
    tippingPoints.c[1, colMatch ] <- tippingPoints.i[,c]
  }
  
  tippingPointsDF <- cbind(tippingPointsDF,tippingPoints.c)

}

## ----------------------
## Depth Range

depthRangeShifts.1 <- loadRData(paste0(dataDirectory,"/SummaryModels/depthDataPredicted.RData"))
depthRangeShifts.2 <- loadRData(paste0(dataDirectory,"/SummaryModels/depthDataPredictedLossGainRefugia.RData"))
depthRangeShifts.2 <- depthRangeShifts.2[-4,-1]

depthRangeShifts <- cbind(depthRangeShifts.1,depthRangeShifts.2)
depthRange <- numeric(0)

for( c in 2:(ncol(depthRangeShifts))) {
  for(r in 1:nrow(depthRangeShifts)) {
    
    value <- depthRangeShifts[r,c]
    valueName <- paste0(tolower(depthRangeShifts[r,1]),names(depthRangeShifts)[c])
    names(value) <- valueName
    depthRange <- c(depthRange,value)
  }
}

depthRange <- as.data.frame(t(data.frame(depthRange)))
depthRange <- data.frame(name=species,depthRange)
rownames(depthRange) <- NULL
 
## ----------------------
## Range shifts

rangeShiftsEstimates <- numeric()

for( scenario in scenariosToPredict ) {

  rangeShiftsEstimates.i <- list.files( paste0(dataDirectory,"/SummaryModels/" ), recursive=T, full.names=T, pattern="rangeShiftEstimates")
  rangeShiftsEstimates.i <- rangeShiftsEstimates.i[grepl(scenario,rangeShiftsEstimates.i)]
  
  if(length(rangeShiftsEstimates.i) == 0) { next }
  
  rangeShiftsEstimates.i <- loadRData(rangeShiftsEstimates.i)
  names(rangeShiftsEstimates.i) <- paste0(names(rangeShiftsEstimates.i),".",scenario)
  rangeShiftsEstimates <- c(rangeShiftsEstimates,rangeShiftsEstimates.i)
  
}

rangeShiftsEstimatesDF <- data.frame(name=species,t(data.frame(unlist(rangeShiftsEstimates))))
row.names(rangeShiftsEstimatesDF) <- NULL
           