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

algorithmsSP <- (c(algorithms,"ensemble"))[ sapply(c(algorithms,"ensemble"), function(x) { length(list.files(paste0(dataDirectory,"/SummaryModels/"), pattern=paste0(x,"Performance"))) > 0 } )]
performanceDF <- data.frame(name=species)

for( algorithmType in algorithmsSP ) {
  performance.i <- list.files(paste0(dataDirectory,"/SummaryModels/"), pattern=paste0(algorithmType,"Performance.RData"), full.names = TRUE)
  performance.i <- loadRData(paste0(dataDirectory,"/SummaryModels/",algorithmType,"Performance.RData"))
  colnames(performance.i) <- paste0(colnames(performance.i),".",algorithmType)
  performanceDF <- cbind(performanceDF,performance.i)
}
row.names(performanceDF) <- NULL

## ----

if( "ensemblePerformanceReachable.RData" %in% list.files(paste0(dataDirectory,"/SummaryModels/") ) ) {
  performance.r <- loadRData(paste0(dataDirectory,"/SummaryModels/ensemblePerformanceReachable.RData"))[,-c(1)]
  names(performance.r) <- paste0(names(performance.r),".ensembleReachab")
  performanceDF <- cbind(performanceDF,performance.r)
}

## ----

algorithmsSP <- (c(algorithms,"ensemble"))[ sapply(c(algorithms,"ensemble"), function(x) { length(list.files(paste0(dataDirectory,"/Models/"), pattern=paste0(x,".RData"))) > 0 } )]
performanceDF.i <- data.frame(name=species)

for( algorithmType in algorithmsSP ) {
  performance.i <- list.files(paste0(dataDirectory,"/Models/"), pattern=paste0(algorithmType,".RData"), full.names = TRUE)
  performance.i <- performance.i[which(!grepl("reduced",performance.i))]
  performance.i <- loadRData(performance.i)$performance
  colnames(performance.i) <- paste0(colnames(performance.i),".",algorithmType)
  performanceDF.i <- cbind(performanceDF.i,performance.i)
}
row.names(performanceDF.i) <- NULL

performanceDF <- cbind(performanceDF.i,performanceDF)

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

depthRangeShifts <- loadRData(paste0(dataDirectory,"/SummaryModels/depthDataPredicted.RData"))
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
           