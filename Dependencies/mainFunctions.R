# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
# One Pipeline for Modelling the distribution of marine Species
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

packages.to.use <- c(  
  "datawizard"
  ,"DescTools"
      ,"credentials"
      ,"SDMtune"
      ,"modEvA"
      ,"blockCV"
      ,"ecospat" 
      ,"plyr"
      ,"raster"
      ,"sf"
      ,"maptools"
      ,"gridExtra"
      ,"xgboost"
      ,"devtools"
      ,"patchwork"
      ,"usdm"
      ,"sm"
      ,"sf"
      ,"SDMTools"
      , "bcp"
      , "spThin"
      , "ecodist"
      , "rnaturalearth"
      , "ggplot2"
      , "mboost"
      , "raster"
      , "dismo"
      , "gbm"
      , "sp"
      , "parallel"
      , "doParallel"
      , "biganalytics"
      , "rgeos"
      , "rgdal"
      , "ENMeval"
      , "remotes"
      ,"leaflet"
      , "monmlp"
      , "okara"
      , "FNN"
      #,"RStoolbox" 
      )

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  
 # if( ! package %in% rownames(installed.packages()) ) { stop("From folder !") }

  if( ! package %in% rownames(installed.packages()) & package == "okara" ) { remotes::install_github("omerkara/okara") }
  if( ! package %in% rownames(installed.packages()) & package == "RStoolbox" ) { install_github("bleutner/RStoolbox") }
  if( ! package %in% rownames(installed.packages()) & package == "USE" ) { devtools::install_github("danddr/USE") }
  if( ! package %in% rownames(installed.packages()) & package == "SDMTools" ) { install_version("SDMTools", "1.1-221") }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package , type = "source") }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

# ---------------------

exportPackages <- function() {
  
  if( ! "miniCRAN" %in% rownames(installed.packages()) ) { install.packages( "miniCRAN" ) }
  library(miniCRAN)
  
  makeRepo(packages.to.use, path = "Dependencies/Packages/", type = "source", Rversion = paste0(R.Version()$major,".",R.Version()$minor),quiet=TRUE)
  cat("Packages export:", "TRUE","\n")
  
}

# ---------------------

interleave <- function(x, y)
{
  m <- length(x)
  n <- length(y)
  xi <- yi <- 1
  len <- m + n
  err <- len %/% 2
  res <- vector()
  for (i in 1:len)
  {
    err <- err - m
    if (err < 0)
    {
      
      res[i] <- x[xi]
      xi <- xi + 1
      err <- err + len
    } else
    {
      res[i] <- y[yi]
      yi <- yi + 1
    }
  }
  res
}

# ---------------------

exportRequirements <- function() {
  
  if( file.exists("requirements.txt") ) { file.remove("requirements.txt")  }
  
  requirements <- data.frame()
  for(package in packages.to.use) {
    lineToWrite <- paste0(package,"==",packageVersion(package))
    requirements <- rbind(requirements,lineToWrite)
  }
  
  cat("Requirements file:", "TRUE","\n")
  write.table(requirements,file="requirements.txt",row.names = FALSE,col.names = FALSE, quote = FALSE)
  return(NULL)
  
}

# ---------------------

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## -----------------------

ensembleModels <- function(modelDataList,rasterLayers,type) {
  
  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1
  cells <- Which(!is.na(shape), cell=TRUE)
  matrixLayers <- as.data.frame(rasterLayers[cells])
  
  ensemble <- matrix(NA,ncol=length(modelDataList),nrow=length(cells))
  ensemble.performance <- numeric(0)
  
  for(m in 1:length(modelDataList)) {

    model.i <- loadRData(modelDataList[m])
    modelDataList.i <- model.i$models
    ensemble.performance <- c(ensemble.performance,model.i$performance$auc)
      
    prediction <- matrix(NA,ncol=length(modelDataList.i),nrow=length(cells))
    matrixLayers.i <- matrixLayers
    
    for( k in 1:length(modelDataList.i) ) {
      
      if( class(modelDataList.i[[k]]) == "xgb.Booster") { 
        matrixLayers.i = xgb.DMatrix(data = data.matrix(matrixLayers[modelDataList.i[[k]]$feature_names]), label = rep(0,nrow(matrixLayers))) 
      }
      
      prediction[,k] <- predict(modelDataList.i[[k]],matrixLayers.i, type="response")
    }
    ensemble[,m] <- apply(prediction,1,mean)
    
  }
  
  shape[] <- NA
  rast <- shape
  
  if( type == "weiAverage") { rast[cells] <- apply(ensemble, 1, weighted.mean, w=ensemble.performance) }
  if( type == "average") { rast[cells] <- apply(ensemble,1,mean) }
  
  rast[rast < 0] <- 0
  rast[rast > 1] <- 1
  
  rast.sd <- shape
  rast.sd[cells] <- apply(ensemble,1,sd)
  rast.sd[rast.sd < 0] <- 0
  rast.sd[rast.sd > 1] <- 1
  
  names(rast) <- "layer"
  names(rast.sd) <- "layer"
  
  res <- list(prediction=rast,prediciton.sd=rast.sd)
  return(res)

}

## -----------------------

drawPolygon <- function(r1) {
  
  plot(r1)
  poly <- spatstat.geom::clickpoly(add=TRUE)
  
  p = Polygon(cbind(poly$bdry[[1]]$x,poly$bdry[[1]]$y))
  ps = Polygons(list(p),1)
  sps = SpatialPolygons(list(ps))
  plot(sps,add=TRUE,col="#727272")
  
  crs(sps) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  sps <- as(sps, "SpatialPolygons")
  
  return( sps )
  
}

## -----------------------

predictedPerformance <- function(prediction,modelData,cvIndex,type="regular") {
  
  observed <- modelData@pa
  predicted <- raster::extract(prediction,modelData@coords)
  predicted.accuracy <- modelPerformance(observed,predicted,cvIndex)
  predicted.accuracy.threshold <- predicted.accuracy$threshold
  predicted.accuracy.boyce <- predicted.accuracy$boyce
  
  if( type == "threshold" ) {
    observed[which( observed == 0 & predicted >= predicted.accuracy$threshold )] <- 1
    predicted.accuracy <- modelPerformance(observed,predicted,cvIndex)
    predicted.accuracy$threshold <- predicted.accuracy.threshold
    predicted.accuracy$boyce <- predicted.accuracy.boyce
  }
  
  return(predicted.accuracy)
  
}

## -----------------------

reachableRange <- function(prediction,dispersalFactor,nonReachableCells=NULL) {
  
  if(!is.null(nonReachableCells)) { prediction[nonReachableCells] <- 0 }
  
  ensembleReclassReach <- aggregate(prediction,dispersalFactor, fun=max)
  ensembleReclassReach[ensembleReclassReach == 0] <- NA
  ensembleReclassRegions <- raster::clump(ensembleReclassReach)
  regionsUsed <- raster::extract(ensembleReclassRegions,occurrenceRecords)
  
  ensembleReclassReach <- ensembleReclassRegions %in% unique(regionsUsed[!is.na(regionsUsed)])
  ensembleReclass[ Which(ensembleReclass == 1, cell=T)[which(raster::extract(ensembleReclassReach,xyFromCell( ensembleReclass , Which(ensembleReclass == 1, cell=T))) == 0)] ] <- 0

  return(ensembleReclass)
  
}

## -----------------------

modelPerformance <- function(observed,predicted,index) {
  
    keepData <- which(!is.na(observed) & !is.na(predicted))
    observed <- observed[keepData]
    predicted <- predicted[keepData]
    
    predicted.accuracy <- SDMTools::accuracy( observed , predicted , threshold = 100 )
    predicted.accuracy$deviance <- modEvA::Dsquared(obs = observed, pred=predicted, family = "binomial")
    predicted.accuracy$aicc <- AIC(glm(observed~predicted))

    boyce.a <- predicted
    boyce.a <- boyce.a[!is.na(boyce.a)]
    boyce.p <- predicted[which(observed==1)]
    boyce.p <- boyce.p[!is.na(boyce.p)]
    
    if( length(boyce.p) == 1) { boyce.p <- c(boyce.p-0.01,boyce.p-0.005,boyce.p,boyce.p+0.005,boyce.p+0.01) }
    if( length(unique(boyce.p)) == 1) { boyce.p <- sapply(boyce.p, function(x) {   x + sample(c(0.001,-0.001), 1)}) }
    boyce.value <- ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE,rm.duplicate = TRUE)
    predicted.accuracy$boyce <- boyce.value$cor
  
    if(index == "boyce" ) { predicted.accuracy <- predicted.accuracy[which.min(abs(predicted.accuracy$threshold - min(  boyce.value$HS[ max( ifelse( length(which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))]) > 0 , which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))] , 1)) ]  ))),] }
    if(index == "tss" ) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
    if(index == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    if(index == "deviance") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    if(index == "minimumTrain") { predicted.accuracy <- predicted.accuracy[max(which(predicted.accuracy$sensitivity == 1)),] }
    if(index == "minimumTrain95") { predicted.accuracy <- predicted.accuracy[max(which(predicted.accuracy$sensitivity >= 0.95 )),] }
    
    predicted.accuracy$tss <- predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1
    names(predicted.accuracy)[2] <- "auc"
    
    ## -------

    return(predicted.accuracy)
  
}

## -----------------------

modelGridSearch <- function(modelType,hyperParam,modelData,cvFolds,monotonicity,cvIndex){
  
  combP <- expand.grid(hyperParam)
  combPResults <- data.frame()
  
  for( c in 1:nrow(combP) ) {
    
    performance <- data.frame()
    
    for( k in 1:ncol(cvFolds$train) ) {
      
      trainData <- data.frame(pa=modelData@pa[cvFolds$train[,k]],modelData@data[cvFolds$train[,k],])
      trainData <- trainData[c("pa",names(which(apply(trainData[,-1],2,function(x) { length(unique(x)) / length(x)}) > 0.05) ))]
      testData <- data.frame(pa=modelData@pa[cvFolds$test[,k]],modelData@data[cvFolds$test[,k],])
      testData <- testData[names(trainData)]
      testDataPA <- testData$pa
      
      monotonicity.i <- monotonicity[match( names(trainData[-1]) , names(monotonicity))]
      
      if(modelType == "BRT"){
        if( "n.trees" %in% names(combP)) { n.trees = combP[c,"n.trees"] } else { n.trees = 100 }
        if( "interaction.depth" %in% names(combP)) { interaction.depth = combP[c,"interaction.depth"] } else { interaction.depth = 2 }
        if( "shrinkage" %in% names(combP)) { shrinkage = combP[c,"shrinkage"] } else { shrinkage = 0.1 }
        model <- NULL
        tryCatch( model <- gbm(pa ~ ., distribution = "bernoulli", data=trainData, var.monotone=monotonicity.i, n.trees = n.trees, interaction.depth = interaction.depth, shrinkage = shrinkage, cv.folds = 0, n.minobsinnode = 5, verbose = FALSE), error=function(e) { Error <<- TRUE })
      }

      if(modelType == "MBOOST"){
        if( "df" %in% names(combP)) { df = combP[c,"df"] } else { df = 10 }
        if( "mstop" %in% names(combP)) { mstop = combP[c,"mstop"] } else { mstop = 50 }
        if( "shrinkage" %in% names(combP)) { shrinkage = combP[c,"shrinkage"] } else { shrinkage = 0.1 }
        trainData[,1] <- as.factor(trainData[,1])
        constr_mono <- as.character(monotonicity.i)
        constr_mono[constr_mono == "-1"] <- "decreasing"
        constr_mono[constr_mono == "1"] <- "increasing"
        rhs <- paste(c(paste("bmono(", colnames(trainData)[-1], ", constraint = \"", constr_mono,"\", df = ",df,")", sep = "")),collapse = " + ")
        fm_mono <- as.formula(paste(colnames(trainData)[1], " ~ ", rhs, collapse = ""))
        ctrl <- boost_control(mstop = mstop, trace = FALSE, nu=shrinkage,stopintern=TRUE)
        model <- NULL
        tryCatch( model <- mboost(fm_mono, data = trainData, control = ctrl, family = Binomial(type = "adaboost" , link = "logit" )) , error=function(e) { Error <<- TRUE })
      }
      
      if(modelType == "XGBOOST"){
        if( "max.depth" %in% names(combP)) { max.depth = combP[c,"max.depth"] } else { max.depth = 2 }
        if( "gamma" %in% names(combP)) { gamma = combP[c,"gamma"] } else { gamma = 0 }
        if( "nrounds" %in% names(combP)) { nrounds = combP[c,"nrounds"] } else { nrounds = 50 }
        trainData = xgb.DMatrix(data = data.matrix(trainData[,-1]), label = trainData[,1], nthread = 1)
        testData = xgb.DMatrix(data = data.matrix(testData[,-1]), label = testData[,1], nthread = 1)
        model <- NULL
        tryCatch( model <- xgboost(data = trainData, monotone_constraints=monotonicity.i, max_depth = max.depth, gamma=gamma, nrounds = nrounds , verbose = 0, nthread = 1, objective="binary:logistic") , error=function(e) { Error <<- TRUE })
      }
      
      if( is.null(model) ) { next }
      
      predTest <- predict(model,testData, type="response")
      performance <- rbind(performance,modelPerformance(testDataPA,predTest,cvIndex))
      
    }
    
    if( nrow(performance) > 0) {
    combPResults <- rbind(combPResults,data.frame(combP[c,],t(data.frame(apply(performance,2,mean)))))
    }
    
  }
  
  return(combPResults)
  
}

## -----------------------

modelTrainCV <- function(modelType,modelData,cvFolds,monotonicity,gridMatrix,cvIndex){
  
  combPResults <- data.frame()
  performance <- data.frame()
  prediction <- data.frame(matrix(NA,ncol=ncol(cvFolds$train),nrow=nrow(modelData@data)))
  models <- list()
  
  for( k in 1:ncol(cvFolds$train) ) {
    
    trainData <- data.frame(pa=modelData@pa[cvFolds$train[,k]],modelData@data[cvFolds$train[,k],])
    trainData <- trainData[c("pa",names(which(apply(trainData[,-1],2,function(x) { length(unique(x)) / length(x)}) > 0.05) ))]
    
    testData <- data.frame(pa=modelData@pa[cvFolds$test[,k]],modelData@data[cvFolds$test[,k],])
    testData <- testData[names(trainData)]
    testDataPA <- testData$pa
    fullData <- data.frame(pa=modelData@pa,modelData@data)
    fullData <- fullData[names(trainData)]
    
    monotonicity.i <- monotonicity[match( names(trainData[-1]) , names(monotonicity))]
    
    model <- NULL
    cvIndexSorting <- sort(gridMatrix[,cvIndex], decreasing = TRUE, index.return=T)$ix
    sortTrial <- 0
    
    while( is.null(model)) {
      
      sortTrial <- sortTrial + 1
      if((nrow(gridMatrix) == sortTrial) ) { break }
      
      if(modelType == "BRT"){
        hyperParam <- list(n.trees=gridMatrix[cvIndexSorting[sortTrial],"n.trees"],interaction.depth = gridMatrix[cvIndexSorting[sortTrial],"interaction.depth"], shrinkage = gridMatrix[cvIndexSorting[sortTrial],"shrinkage"])
        if( "n.trees" %in% names(hyperParam)) { n.trees = hyperParam$n.trees } else { n.trees = 100 }
        if( "interaction.depth" %in% names(hyperParam)) { interaction.depth = hyperParam$interaction.depth } else { interaction.depth = 2 }
        if( "shrinkage" %in% names(hyperParam)) { shrinkage = hyperParam$shrinkage } else { shrinkage = 0.1 }
        tryCatch( model <- gbm(pa ~ ., distribution = "bernoulli", data=trainData, var.monotone=monotonicity.i, n.trees = n.trees, interaction.depth = interaction.depth, shrinkage = shrinkage, cv.folds = 0, verbose = FALSE) , error=function(e) { Error <<- TRUE })
        
        if(is.null(model)) { tryCatch( model <- gbm(pa ~ ., distribution = "bernoulli", data=rbind(trainData,trainData), var.monotone=monotonicity.i, n.trees = n.trees, interaction.depth = interaction.depth, shrinkage = shrinkage, cv.folds = 0, verbose = FALSE) , error=function(e) { Error <<- TRUE }) }
        
      }
      
      if(modelType == "MBOOST"){
        hyperParam <- list(df=gridMatrix[cvIndexSorting[sortTrial],"df"],mstop = gridMatrix[cvIndexSorting[sortTrial],"mstop"], shrinkage = gridMatrix[cvIndexSorting[sortTrial],"shrinkage"])
        if( "df" %in% names(hyperParam)) { df =hyperParam$df } else { df = 10 }
        if( "mstop" %in% names(hyperParam)) { mstop = hyperParam$mstop } else { mstop = 50 }
        if( "shrinkage" %in% names(hyperParam)) { shrinkage = hyperParam$shrinkage } else { shrinkage = 0.1 }
        trainData[,1] <- as.factor(trainData[,1])
        constr_mono <- as.character(monotonicity.i)
        constr_mono[constr_mono == "-1"] <- "decreasing"
        constr_mono[constr_mono == "1"] <- "increasing"
        rhs <- paste(c(paste("bmono(", colnames(trainData)[-1], ", constraint = \"", constr_mono,"\", df = ",df,")", sep = "")),collapse = " + ")
        fm_mono <- as.formula(paste(colnames(trainData)[1], " ~ ", rhs, collapse = ""))
        ctrl <- boost_control(mstop = mstop, trace = FALSE, nu=shrinkage,stopintern=TRUE)
        tryCatch( model <- mboost(fm_mono, data = trainData, control = ctrl, family = Binomial(type = "adaboost" , link = "logit" )) , error=function(e) { Error <<- TRUE })
      }
      
      if(modelType == "XGBOOST"){
        hyperParam <- list(max.depth=gridMatrix[cvIndexSorting[sortTrial],"depth"],gamma = gridMatrix[cvIndexSorting[sortTrial],"gamma"], nrounds = gridMatrix[cvIndexSorting[sortTrial],"nrounds"])
        if( "max.depth" %in% names(hyperParam)) { max.depth = hyperParam$max.depth } else { max.depth = 2 }
        if( "gamma" %in% names(hyperParam)) { gamma = hyperParam$gamma } else { gamma = 0 }
        if( "nrounds" %in% names(hyperParam)) { nrounds = hyperParam$nrounds } else { nrounds = 50 }
        trainData = xgb.DMatrix(data = data.matrix(trainData[,-1]), label = trainData[,1])
        testData = xgb.DMatrix(data = data.matrix(testData[,-1]), label = testData[,1])
        fullData = xgb.DMatrix(data = data.matrix(fullData[,-1]), label = fullData[,1])
        tryCatch( model <- xgboost(data = trainData, monotone_constraints=monotonicity.i, max_depth = max.depth, gamma=gamma, nrounds = nrounds , verbose = 0, objective="binary:logistic") , error=function(e) { Error <<- TRUE })
      }
      
    }

    if( ! is.null(model) ) {
      
      models <- c(models,list(model))
      predTest <- predict(model,testData, type="response")
      prediction[,k] <- predict(model,fullData, type="response")
      performance <- rbind(performance,modelPerformance(testDataPA,predTest,cvIndex))
   
      combPResults <- list(models=models,
                           performanceCV=apply(performance,2,mean),
                           performanceCV.sd=apply(performance,2,sd),
                           performance=modelPerformance(modelData@pa,apply(prediction,1,mean,na.rm=T),cvIndex))
      
    }

  }
  
   return(combPResults)
  
}

## -----------------------

prepareModelData <- function(p,a,env) {
  
  m <- subset(env,1)
  
  p.i <- raster::extract(m,p)
  p <- p[which(!is.na(p.i)),]
  
  p <- xyFromCell(m,cellFromXY(m,p))
  p <- unique(p)
  
  a.i <- raster::extract(m,a)
  a <- a[which(!is.na(a.i)),]
  
  a <- xyFromCell(m,cellFromXY(m,a))
  a <- unique(a)

  tryCatch( modelData <- prepareSWD(species = "Model species", p = p, a = a, env = env), error=function(e) { modelData <<- prepareSWD(species = "Model species", p = p, a = a, env = terra::rast(env)) })
  return(modelData)
  
}

## -----------------------

getBlocks <- function(modelData,rasterLayers,corrDistance,k) {

  library(blockCV)
  library(SDMtune)
  
  occ <- modelData@coords
  occ <- SpatialPoints(coords = occ[, 1:2], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") )
  occ$occ <- modelData@pa
  
  sb <- NULL
  while(is.null(sb)) {
    tryCatch( sb <- blockCV::cv_spatial(x = occ, # sf or SpatialPoints of sample data
                               column = "occ", # the response column (binary or multi-class)
                               r = as(rasterLayers, "SpatRaster"), # a raster for background (optional)
                               size = corrDistance * 1000 , # size of the blocks in metres
                               k = k, # number of folds
                               hexagon = TRUE, # use hexagonal blocks - defualt
                               selection = "systematic", # random blocks-to-fold
                               iteration = 50, # find evenly dispersed folds
                               biomod2 = FALSE,
                               progress = FALSE) , error=function(e) { Error <<- TRUE })
  
    corrDistance <- corrDistance - (corrDistance * 0.1)
    if(corrDistance < 10) { stop("Error :: 915")}
    
  }
  
  train <- matrix(FALSE,ncol=k,nrow=length(occ))
  test <- matrix(FALSE,ncol=k,nrow=length(occ))
  for( i in 1:k) {
    test[ sb$folds_list[[i]][[2]] , i] <- TRUE
    train[ sb$folds_list[[i]][[1]]  , i] <- TRUE
  }
  
  train.vect.p <- which(apply(train[which(modelData@pa==1),],2,sum) > 0 )
  test.vect.p <- which(apply(test[which(modelData@pa==1),],2,sum) > 0 )
  foldsList <- intersect(train.vect.p,test.vect.p)
  
  if( length(which( ! 1:k %in% train.vect.p)) > 0 ) { 
    
    for( k.i in which( ! 1:k %in% train.vect.p) ) {
      train[which(modelData@pa==1),k.i] <- train[which(modelData@pa==1),sample(x=c(train.vect.p,train.vect.p),size=1)]
    }
    
    train.vect.p <- which(apply(train[which(modelData@pa==1),],2,sum) > 0 )
    test.vect.p <- which(apply(test[which(modelData@pa==1),],2,sum) > 0 )
    foldsList <- intersect(train.vect.p,test.vect.p)

  }
  
  if( length(which( ! 1:k %in% test.vect.p)) > 0 ) { 
    
    for( k.i in which( ! 1:k %in% test.vect.p) ) {
      test[which(modelData@pa==1),k.i] <- test[which(modelData@pa==1),sample(x=c(test.vect.p,test.vect.p),size=1)]
    }
    
    train.vect.p <- which(apply(train[which(modelData@pa==1),],2,sum) > 0 )
    test.vect.p <- which(apply(test[which(modelData@pa==1),],2,sum) > 0 )
    foldsList <- intersect(train.vect.p,test.vect.p)
    
  }
  
  figureBlocks <- sb$blocks$blocks
  figureBlocksK <- sb$blocks$folds
  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1

  figurePartitioning <- cv_plot(cv = sb,
                                x = occ, 
                                r = as(shape, "SpatRaster"),
                                points_alpha = 0.5,
                                num_plots=foldsList,
                                nrow = 2)

  folds <- list( figureBlocks=figureBlocks,figureBlocksK=figureBlocksK, figurePartitioning=figurePartitioning, train = train[,foldsList] , test = test[,foldsList] )

  return(folds)
  
}

## -----------------------

correctDateLineMap <- function(spatialObjectSF) {
  
  if(! TRUE %in% (class(spatialObjectSF) == "sf")) { stop("Object must be of class sf")}
  
  dggsOcean <- spatialObjectSF
  dggsOcean <- as(dggsOcean, "Spatial")
  listToRemove <- numeric()
  listToAdd <- list()
  
  for(i in 1:length(dggsOcean)) {
    
    if( extent(dggsOcean[i,])[1] < 0 & extent(dggsOcean[i,])[2] > 0 & max(extent(dggsOcean[i,])[1:2])-min(extent(dggsOcean[i,])[1:2]) > 180 ) { 
      
      listToRemove <- c(listToRemove,i)
      
      coords1 <- fortify(dggsOcean[i,])[,1:2]
      coords2 <- fortify(dggsOcean[i,])[,1:2]
      coords1[coords1[,1] < 0 ,1] <- 180
      coords2[coords2[,1] > 0 ,1] <- -180
      coords2 <- rbind(coords2,data.frame(long=-180,lat=max(coords2[,2])))
      coords2 <- rbind(coords2,data.frame(long=-180,lat=min(coords2[,2])))
      coords2 <- coords2[chull(coords2),] 
      
      coords1 <- spPolygons(as.matrix(coords1))
      coords2 <- spPolygons(as.matrix(coords2))
      
      listToAdd.i <- unionSpatialPolygons(union(coords1, coords2), c(1,1))
      listToAdd.i <- SpatialPolygonsDataFrame(listToAdd.i,as.data.frame(dggsOcean[i,]), match.ID = FALSE)
      
      if( length(listToAdd) >  0 ) { listToAdd <- bind(listToAdd, listToAdd.i) }
      if( length(listToAdd) == 0 ) { listToAdd <- listToAdd.i }
      
    }
    
    if( (extent(dggsOcean[i,])[1] > 0 & extent(dggsOcean[i,])[2] < 0) ) { stop(i) }
    
  }
  
  if( length(listToRemove) > 0) {
    
    dggsOcean <- dggsOcean[-listToRemove,]
    dggsOcean <- bind(dggsOcean,listToAdd)
    
  }
  
  dggsOcean <- st_as_sf(dggsOcean)
  
  return(dggsOcean)      
}

## -----------------------

plotMap <- function(coordinates,radius,color) {
  
  set.seed(42)
  
  m <- leaflet()
  m <- addTiles(m)
  # m <- addMarkers(m, lng=coordinates[,1], lat=coordinates[,2], popup=paste0( "Species record ") , icon = greenLeafIcon)
  
  m <- addCircleMarkers(m, lng=coordinates[,1], lat=coordinates[,2], 
                        popup=paste0( "Species record ") , 
                        radius = radius, color = color , 
                        stroke = FALSE, fillOpacity = 0.5 )
  
  return(m)
  
}

## -----------------------

processLayers <- function(rasterLayers,occurrenceRecords,regionBuffer,minDepth,maxDepth,intertidal) {
  
  if( ! is.null(minDepth) ) { if( minDepth == "NULL" ) { minDepth <- NULL } }
  if( ! is.null(maxDepth) ) { if( maxDepth == "NULL" ) { maxDepth <- NULL } }

  rastersBathymetryID <- NULL
  rastersIntertidalID <- NULL
  
  rasters <- rasterLayers
  
  if( ! is.null(regionBuffer)  ) {
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    colnames(sink.points.poly) <- dataRecordsNames
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    sink.points.poly.list <- list()
    for( l in 1:length(sink.points.poly) ) {
      sink.points.poly.list <- c(sink.points.poly.list, list(gBuffer( sink.points.poly[l,], width=regionBuffer )))
    }
    sink.points.poly <- do.call("rbind",c(args = sink.points.poly.list, makeUniqueIDs = TRUE))
    row.names(sink.points.poly) <- as.character(1:length(sink.points.poly))
    sink.points.poly <- SpatialPolygonsDataFrame(sink.points.poly, data.frame(val=as.character(1:length(sink.points.poly))))
    sink.points.poly$val <- 1
    sink.points.poly.base <- gUnaryUnion(sink.points.poly, id = sink.points.poly$val)
    
    rasters <- crop(rasters,sink.points.poly.base)
    rasters <- mask(rasters,sink.points.poly.base)
    
  }

  if( ! is.null(regionBufferPoly)  ) { 
    
    regionBufferPoly.shp <- shapefile(regionBufferPoly)
    regionBufferPoly.shp <- regionBufferPoly.shp[,regionBufferPolyLevel]
    names(regionBufferPoly.shp) <- "Level"
    regionBufferPoly.shp.DF <- as.data.frame(regionBufferPoly.shp)
    regionBufferPoly.shp <- gUnaryUnion(regionBufferPoly.shp, id = regionBufferPoly.shp@data$Level, checkValidity=2L)
    mathingVector <- unlist(sapply(1:length(regionBufferPoly.shp), function(x) which(regionBufferPoly.shp.DF$Level ==  slot(slot(regionBufferPoly.shp, "polygons")[[x]],"ID") )[1]))
    regionBufferPoly.shp.DF <- data.frame(Level=regionBufferPoly.shp.DF[mathingVector,])
    rownames(regionBufferPoly.shp.DF) <- unlist(sapply(1:length(regionBufferPoly.shp), function(x) slot(slot(regionBufferPoly.shp, "polygons")[[x]],"ID") ))
    regionBufferPoly.shp <- SpatialPolygonsDataFrame(regionBufferPoly.shp, regionBufferPoly.shp.DF )
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    colnames(sink.points.poly) <- dataRecordsNames
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    
    closestPoly <- apply(gDistance(sink.points.poly, regionBufferPoly.shp,byid=TRUE),2,which.min)
    closestPoly <- unique(regionBufferPoly.shp[unique(closestPoly),]$Level)
    
    closestPolyAdjacent <- gTouches(regionBufferPoly.shp, byid=TRUE)[which(row.names(gTouches(regionBufferPoly.shp, byid=TRUE)) %in% closestPoly),]
    closestPolyAdjacent <- names(which(apply(closestPolyAdjacent,2,sum) != 0))
    
    regionBufferPoly.shp <- regionBufferPoly.shp[regionBufferPoly.shp$Level %in% unique(c(closestPoly,closestPolyAdjacent)),]
    rasters <- crop(rasters,regionBufferPoly.shp)
    rasters <- mask(rasters,regionBufferPoly.shp)
    
  }
  
  ## -----------------------
  
  regions <- subset(rasters,1)
  regions[!is.na(regions) ] <- 1
  regions <- aggregate(regions,5, fun=max, na.rm=TRUE)
  regions <- clump(regions)
  regionsUsed <- unique(raster::extract(regions,occurrenceRecords))
  regionsUsed <- regionsUsed[ !is.na(regionsUsed) ]
  regions[ ! regions %in% regionsUsed ] <- NA
  rastersLoc <- xyFromCell( subset(rasters,1) , Which( !is.na(subset(rasters,1)), cell=T))
  rastersLocNA <- which(is.na(raster::extract(regions,rastersLoc)))
  rasters[ Which( !is.na(subset(rasters,1)), cell=T)[rastersLocNA] ] <- NA
  
  ## -----------------------
  
  if( ! is.null(intertidal) ) {
    
    intertidalLayer <- raster(intertidal)
    intertidalLayer <- raster::crop(intertidalLayer, rasters)
    rastersIntertidalID <- Which( is.na(intertidalLayer), cells=TRUE)
    
  }
  
  ## -----------------------
  
  if( ! is.null(minDepth) | ! is.null(maxDepth) ) {
    
    if( is.null(minDepth) ) { minDepth <- 0 }
    if( is.null(maxDepth) ) { maxDepth <- 99999 }
    
    if( length(bathymetryDataLayer) > 1) {
      minDepthRaster <- raster(bathymetryDataLayer[grepl("Min",bathymetryDataLayer)])
      maxDepthRaster <- raster(bathymetryDataLayer[grepl("Max",bathymetryDataLayer)])
      range <- c(minDepth, maxDepth) * (-1)
      shapeRaster <- minDepthRaster <= max(range) & maxDepthRaster >= min(range)
      bathymetry <- crop(shapeRaster, rasters)
      bathymetry <- mask(bathymetry, subset(rasters,1))
      bathymetry[bathymetry == 0] <- NA
      rastersBathymetryID <- Which( is.na(bathymetry), cells=TRUE)
      
    }

    if( length(bathymetryDataLayer) == 1) {
      
      minDepth <- minDepth * (-1)
      maxDepth <- maxDepth * (-1)
      
      rclmat <- data.frame(from=c(-99999,maxDepth,minDepth),to=c(maxDepth,minDepth,0),value=c(NA,1,NA))
      rclmat <- rclmat[rclmat$from != rclmat$to,]
      bathymetry <- raster(bathymetryDataLayer[1])
      bathymetry <- crop(bathymetry, rasters)
      bathymetry <- reclassify(bathymetry, as.matrix(rclmat))
      rastersBathymetryID <- Which( is.na(bathymetry), cells=TRUE)

    }
   
  }
  
  ## -----------------------
  
  rastersID <- unique(c(rastersIntertidalID,rastersBathymetryID))
  rasters[rastersID] <- NA
  
  if(class(rasters) == "RasterBrick") { rasters <- stack(rasters)}
  return(rasters)
  
}

## -----------------------

dropNoVariationLayers <- function(rasterLayers,records) {
  
  cave <- function(x) { length(unique(x)) / length(x) }
  varRasterLayers <- which( apply( raster::extract(rasterLayers,records) ,2,cave) > 0.025)
  rasterLayers <- subset(rasterLayers,varRasterLayers)
  return(rasterLayers)
  
}

## -----------------------

relocateNACoords <- function(occurrenceRecords,rasterLayers,relocateType,relocateSpeciesDistance,relocateSpeciesDepth) {
  
  occurrenceRecords.bk <- occurrenceRecords
  
  set.seed(42)
  
  if(relocateType == "nonDistance") {
    
    vals <- raster::extract(subset(rasterLayers,1),occurrenceRecords)
    occurrenceRecords <- occurrenceRecords[which(!is.na(vals)),]
    
  }
  
  if(relocateType == "distance") {
    
    if( relocateSpeciesDepth ) { 
      
      if( is.null(minDepth) ) { minDepth <- 0 }
      if( is.null(maxDepth) ) { maxDepth <- 99998 }
      
      minDepth.i <- minDepth * (-1)
      maxDepth.i <- maxDepth * (-1)
      
      rclmat <- data.frame(from=c(-99999,maxDepth.i,minDepth.i),to=c(maxDepth.i,minDepth.i,0),value=c(NA,1,NA))
      rclmat <- rclmat[rclmat$from != rclmat$to,]
      
      shape <- subset(rasterLayers,1)
      
      if( length(bathymetryDataLayer) == 1 ) { bathymetry <- raster(bathymetryDataLayer) }
      if( length(bathymetryDataLayer) >= 2 ) { bathymetry <- calc(stack(bathymetryDataLayer),mean) }
      
      bathymetry <- crop(bathymetry, shape )
      bathymetry <- mask(bathymetry,shape )
      
      shape <- reclassify(bathymetry, as.matrix(rclmat))
      
    }
    
    if( ! relocateSpeciesDepth ) { shape <- subset(rasterLayers,1) } 
    
    to.relocate <- unique(which(is.na(raster::extract(shape,occurrenceRecords))))
    coordinates.to.relocate <- matrix(occurrenceRecords[c(1,2),],ncol=2)
    
    if( nrow(coordinates.to.relocate) > 0 ) { 
      
      old.presences <- occurrenceRecords[ (1:nrow(occurrenceRecords))[! 1:nrow(occurrenceRecords) %in% to.relocate] ,]
      
      correct.points <- xyFromCell(shape,  Which(!is.na(shape), cells=TRUE) ) # result
      
      cat( paste0("Relocating ",length(to.relocate)," Points that were falling out of range"))
      cat( paste0("\n"))
      
      near.cells <- numeric(nrow(coordinates.to.relocate))
      
      for(p in 1:nrow(coordinates.to.relocate)) {
        
        near.cell.p <- spDistsN1( as.matrix(correct.points), as.matrix(coordinates.to.relocate[p,]),longlat=TRUE)
        
        if( near.cell.p[which.min(near.cell.p)] <= sqrt(sum(relocateSpeciesDistance^2,relocateSpeciesDistance^2)) ) {
          
          near.cell.p <- which.min(near.cell.p)
          
        } else {   near.cell.p <- NA }
        
        near.cells[p] <- near.cell.p
        
      }
      
      relocated <- which(!is.na(near.cells))
      
      if( length(relocated) > 0) {
        
        near.cells <- data.frame(Lon=correct.points[near.cells[relocated],1],Lat=correct.points[near.cells[relocated],2])
        colnames(near.cells) <- colnames(old.presences)
        occurrenceRecords <- rbind(old.presences,near.cells)
        
      }
      
      cat( paste0("Done relocating"))
      cat( paste0("\n"))
      
    }
    
    ## -----------------------
    
    if( nrow(coordinates.to.relocate) == 0) { 
      
      cat( paste0("None to Relocate"))
      cat( paste0("\n"))
      
    }
    
    ## -----------------------
    
  }
  
  toRemove <- raster::extract(shape,occurrenceRecords)
  toRemove <- which(is.na(toRemove))
  
  if( length(toRemove) > 0 ) { occurrenceRecords <- occurrenceRecords[-toRemove,] }
  
  if( nrow(occurrenceRecords) == 0 ) { occurrenceRecords <-occurrenceRecords.bk }
  

  return( occurrenceRecords )
  
}

## -----------------------

spatialAutocorrelation <- function(rasterLayers,autocorrelationClassDistance,autocorrelationMaxDistance) {
  
  rasterLayers.i <- as.data.frame(rasterLayers,xy=T,na.rm=T)
  presences.environment <- rasterLayers.i[,which(!colnames(rasterLayers.i) %in% c("x","y"))]
  dataCoords <- rasterLayers.i[,which(colnames(rasterLayers.i) %in% c("x","y"))]
  
  autocorrelationMaxDistance <- 500
  numSample <- 5000
  set.seed(42)
  numSample <- sample(1:nrow(dataCoords), min(numSample,nrow(dataCoords)), replace=FALSE)
  presences.environment <- presences.environment[numSample,]
  presences.environment <- presences.environment[,which(apply(presences.environment,2,var, na.rm=T) != 0)]
  presences.environment <- scale(presences.environment)
  dataCoords <- dataCoords[numSample,]
  
  # --------
  
  space <- spDists(as.matrix(dataCoords),as.matrix(dataCoords),longlat=TRUE)
  data <- as.matrix(ecodist::distance( presences.environment , method = "euclidean"))

  autocorrelationClassDistance <- autocorrelationClassDistance
  n.class <- round(autocorrelationMaxDistance/autocorrelationClassDistance)-1
  
  resultsMatrix <- data.frame(classdistanceFrom=seq(0,autocorrelationMaxDistance-autocorrelationClassDistance,autocorrelationClassDistance),
                              classdistanceTo=seq(autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationClassDistance),
                              R=NA,
                              pVal=NA)
  
  data.d <- as.vector(data)
  space.d <- as.vector(space)
  
  for( i in 1:nrow(resultsMatrix)) {
    
    d1 = resultsMatrix[i,1]
    d2 = resultsMatrix[i,2]
    
    remove <- which(space.d < d1 | space.d > d2)
    data.i <- data.d[-remove]
    space.i <- space.d[-remove]
    
    if( length(data.i) == 0 | length(space.i) == 0) { next }
    
    resultsMatrix[i,3] <- cor.test(data.i, space.i, method ="pearson")$estimate
    resultsMatrix[i,4] <- cor.test(data.i, space.i, method ="pearson")$p.value
    
    if(resultsMatrix[i,3] < 0) { resultsMatrix[i,3] <- 0; resultsMatrix[i,4] <- 1}
    
  }
  
  autocorrelationSignif <- 0.05
  minDistance = min(resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ])
  medianDistance = median(resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ])
  meanDistance = mean(resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ])
  sdDistance = sd(resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ])
  
  if( is.na(minDistance)) { minDistance <- autocorrelationMaxDistance }
  
  sac <- resultsMatrix[,c(2,3,4)]
  colnames(sac) <- c("Distance","R","signif")
  
  figure <- ggplot() +
    geom_line(data=sac, aes(x=Distance, y=R), linetype = "dashed", size=0.25)+
    geom_point(data=sac[sac$signif >= autocorrelationSignif,], aes(x=Distance, y=R), color="gray", size=2) +
    geom_point(data=sac[sac$signif < autocorrelationSignif,], aes(x=Distance, y=R), color="orange", size=3) + 
    labs(x = "Distance (km)") + 
    labs(y = "Pearson' correlation") +
    ggtitle("Multivariate spatial autocorrelation") +
    theme_light() +
    theme(
      panel.grid.major.y = element_blank(),
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_blank(),
      axis.title.y = element_text(vjust = +3),
      axis.title.x = element_text(vjust = -2),
      axis.ticks.y = element_blank()
    )
  
  return( list(figure=figure, minDistance = minDistance , medianDistance = medianDistance , meanDistance = meanDistance, sdDistance = sdDistance  ))
  
}

## -----------------------

spatialAutocorrelationBk <- function(rasterLayers) {
  
  occurrenceRecords.i <- occurrenceRecords
  coordinates(occurrenceRecords.i) <- ~Lon+Lat
  occurrenceRecords.i$pa <- rep(1,length(occurrenceRecords.i))
  crs(occurrenceRecords.i) <- crs(rasterLayers)
  
  rasterLayers.i <- crop(rasterLayers, extent(min(occurrenceRecords[,1])-5,max(occurrenceRecords[,1])+5,min(occurrenceRecords[,2])-5,max(occurrenceRecords[,2])+5))
  rasterLayers.i <- aggregate(rasterLayers.i,4,mean)
  spatialAutocorr <- cv_spatial_autocor(r=rasterLayers.i, x = occurrenceRecords.i, column="pa", rogress = FALSE, num_sample = 1000,  plot = TRUE) 
  
  range_table <- spatialAutocorr$range_table
  range_table$range <- range_table$range / 1000 
  range_table[which.max(range_table$range),"range"] <- range_table[which.max(range_table$range),"range"] / 1000
  
  figure <- ggplot(range_table, aes(x=reorder(layers,range), y=range)) +
    geom_segment( aes(x=reorder(layers,range), xend=layers, y=0, yend=range), size=0.5,  color="black") +
    geom_point( color="orange", size=5, alpha=1) +
    geom_hline(yintercept = median(range_table$range), linetype = "dashed", size=0.5,  color="black", alpha=0.5) +
    xlab("Predictor") +
    ylab("Spatial autocorrelation") +
    theme_light() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )

  # --------

  return( list(figure=figure, minDistance = min(range_table$range) , medianDistance = median(range_table$range)  ))
  
}

## -----------------------

spatialThinning <- function(occurrenceRecords,minDistance,verbose) {
  
  occurrenceRecords <- occurrenceRecords[which(!duplicated(occurrenceRecords)),]
  dataThinning <- data.frame(Name="Sp",occurrenceRecords)
  
  maxBlocks <- min(nrow(occurrenceRecords),10000)
  
  blocks <- data.frame(from=seq(from=1, to=nrow(dataThinning), by=maxBlocks)[-length(seq(from=1, to=nrow(dataThinning), by=maxBlocks))],
                       to=(seq(from=1, to=nrow(dataThinning), by=maxBlocks)-1)[-1])
  blocks[nrow(blocks),2] <- nrow(dataThinning)
  
  if( nrow(blocks) < 2 ) {
    
    blocks <- data.frame(from=1,
                         to=nrow(dataThinning))
    
    
  }
  
  coordinates.t <- data.frame()
  
  for(b in 1:nrow(blocks)) {
          
        thinError <- FALSE
        tryCatch( 
          coordinates.t.b <- thin(dataThinning[blocks[b,1]:blocks[b,2],],
                                lat.col = dataRecordsNames[2],
                                long.col = dataRecordsNames[1],
                                spec.col = "Name",
                                thin.par = minDistance,
                                reps = 1,
                                write.files = FALSE,
                                locs.thinned.list.return = TRUE,
                                verbose = FALSE)[[1]]
          , error=function(e) { thinError <<- TRUE })
        
        if( thinError ) { coordinates.t.b <- occurrenceRecords }
        
    coordinates.t <- rbind(coordinates.t,coordinates.t.b)
    
  }
  
  if(verbose) {
    
    cat( paste0("\n"))
    cat( paste0("\n"))    
    cat( paste0("Input Records: ",nrow(occurrenceRecords)))
    cat( paste0("\n"))
    cat( paste0("Final Records: ",nrow(coordinates.t)))
    
  }
  
  # Remove from main dataset of occurrences
  colnames( coordinates.t ) <- dataRecordsNames
  
  # Remove log file
  if(length(list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt"))){
    file.remove( list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt") ) 
  }
  
  return(coordinates.t)
  
}

## -----------------------

generateBackgroundInformation <- function(rasterLayers,occurrenceRecords,n,spatialAutocorrelationVal) {
  
  options(warn=-1)
  
  shape <- as.data.frame(subset(rasterLayers,1), xy=T,na.rm=T)[,1:2]
  
  parallelChunk <- data.frame(from=seq(1,nrow(shape),by=5000)[-length(seq(1,nrow(shape),by=5000))],
                              to=seq(1,nrow(shape),by=5000)[-1]-1)
  parallelChunk[nrow(parallelChunk),2] <- nrow(shape)
  
  if(nrow(parallelChunk) == 0) { parallelChunk <- data.frame(from=1,to=nrow(shape))}
  
  shape.i <- data.frame()
  for(par in 1:nrow(parallelChunk)) {
    
    shape.i.temp <- spatialThinning(shape[parallelChunk[par,1]:parallelChunk[par,2],],spatialAutocorrelationVal,verbose=FALSE)
    shape.i <- rbind(shape.i,shape.i.temp)
    
  }
  
  n <- min(n,nrow(shape.i))
  
  pa.environment <- raster::extract( rasterLayers,shape.i )
  pa.environment <- scale(pa.environment)
  pa.environment <- pa.environment[,which(!is.na(sapply(1:ncol(pa.environment) , function(x) { var(pa.environment[,x]) })))]
  pa.environment <- kmeans(pa.environment,centers=n,iter.max = 10, nstart = 1)$cluster
  
  pa.clusters <- unique( pa.environment )
  finalKMClusters <- sapply(1:n,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
  shape.i <- shape.i[finalKMClusters,]
  
  return(shape.i)
  options(warn=0)
  
}

## -----------------------

generatePseudoAbsences <- function(occurrenceRecords,rasterLayers,paType) {
  
  biasSurface <- NULL
  occurrenceRecords <- as.data.frame(occurrenceRecords)
  
  if( paType == "random" ) {

    pseudoAbsences <- xyFromCell( subset(rasterLayers,1) , Which( !is.na(subset(rasterLayers,1) ) , cells=T) )
  
    if(paRemoveOverOccurrence) {
      
      occurrenceRecords.sp <- occurrenceRecords
      coordinates(occurrenceRecords.sp) <- ~Lon+Lat
      crs(occurrenceRecords.sp) <- crs(rasterLayers)
      occurrenceRecords.sp <- raster::buffer(occurrenceRecords.sp, width=111000)
      
      pseudoAbsences.sp <- pseudoAbsences
      coordinates(pseudoAbsences.sp) <- ~x+y
      crs(pseudoAbsences.sp) <- crs(rasterLayers)
      pseudoAbsences <- pseudoAbsences[which(is.na(over(pseudoAbsences.sp,occurrenceRecords.sp))),]
    }

  }
  
  if( paType == "mess" ) {

    backgroundInformation <- xyFromCell( rasterLayers , Which(!is.na(subset(rasterLayers,1)), cell=TRUE) )

    proj <- data.frame(long=backgroundInformation[,1],lat=backgroundInformation[,2],raster::extract(rasterLayers, backgroundInformation))
    cal <- data.frame(long=occurrenceRecords[,1],lat=occurrenceRecords[,2],raster::extract(rasterLayers, occurrenceRecords))

    toKeep <- intersect(which(apply(var(proj),1,sum,na.rm=T) != 0) , which(apply(var(cal),1,sum,na.rm=T) != 0))
    
    if(length(toKeep)>0)  {
      proj <- proj[,toKeep]
      cal <- cal[,toKeep]
    }

    messSurface <- ecospat.mess( proj , cal , w="default")
    pseudoAbsencesClasses <- unique(messSurface[,"MESSneg"])[unique(messSurface[,"MESSneg"]) != 0]
    
    pseudoAbsences <- data.frame(Lon=backgroundInformation[which(messSurface[,"MESSneg"] != 0),1],
                                 Lat=backgroundInformation[which(messSurface[,"MESSneg"] != 0),2])
    
    occurrenceRecords.sp <- data.frame(occurrenceRecords)
    coordinates(occurrenceRecords.sp) <- ~Lon+Lat
    crs(occurrenceRecords.sp) <- crs(rasterLayers)
    occurrenceRecords.sp <- raster::buffer(occurrenceRecords.sp, width=111000)
    
    pseudoAbsences.sp <- pseudoAbsences
    coordinates(pseudoAbsences.sp) <- ~Lon+Lat
    crs(pseudoAbsences.sp) <- crs(rasterLayers)
    
    pseudoAbsences <- pseudoAbsences[which(is.na(over(pseudoAbsences.sp,occurrenceRecords.sp))),]
    
  }
  
  if( paType == "kernelDensity" ) {
    
    shape <- subset(rasterLayers,1)
    
    library(sm)
    bias <- cellFromXY(shape, as.matrix(occurrenceRecords))
    cells <- unique(sort(bias))
    kernelXY <- xyFromCell(shape, cells)
    samps <- as.numeric(table(bias))
    KDEsur <- sm.density(kernelXY, weights=samps, display="none", hmult=1, ylim=c(extent(shape)[3],extent(shape)[4]), xlim=c(extent(shape)[1],extent(shape)[2]), nbins=NA)
    KDErast <- SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
    KDErast <- SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, length(KDEsur$estimate))))
    KDErast <- raster(KDErast)
    biasSurface <- raster::resample(KDErast, shape, method="bilinear")

    Bias.surfaceThresh <- raster::extract(biasSurface,occurrenceRecords)
    Bias.surfaceThresh <- Bias.surfaceThresh[complete.cases(Bias.surfaceThresh)]
    Bias.surfaceThresh <- quantile(Bias.surfaceThresh, probs=0.025)
    biasSurface[biasSurface > Bias.surfaceThresh ] <- NA
    biasSurface <- (biasSurface - cellStats(biasSurface,min)) / cellStats(biasSurface - cellStats(biasSurface,min), max)
      
    pseudoAbsences <- as.data.frame( mask(biasSurface,shape),xy=T,na.rm=T)[,1:2]

    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    sink.points.poly <- raster::buffer(sink.points.poly, width=111000)

    if(paRemoveOverOccurrence) {
      pseudoAbsences.sp <- pseudoAbsences
      coordinates(pseudoAbsences.sp) <- ~x+y
      crs(pseudoAbsences.sp) <- crs(rasterLayers)
      pseudoAbsences <- pseudoAbsences[which(is.na(over(pseudoAbsences.sp,occurrenceRecords.sp))),]
    }

  }
  
  if( paType == "sre") {
    
    library(biomod2)
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    sink.points.poly$data <- rep(1,length(sink.points.poly))
    
    pseudoAbsences <- bm_PseudoAbsences(
      resp.var=terra::vect(sink.points.poly),
      expl.var= terra::rast(rasterLayers),
      nb.rep = 1,
      strategy = "sre",
      nb.absences = final.paRatio,
      sre.quant = 0 )
    
    pseudoAbsences <- pseudoAbsences$xy[which(is.na(pseudoAbsences$sp)),]
    
  }
  
  if( paType == "use") {
    
    library(USE)
    library(sf)
  
    prevalence <- nrow(occurrenceRecords) / max(paMinimum,nrow(occurrenceRecords))

    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- crs( rasterLayers )
    sink.points.poly$data <- rep(1,length(sink.points.poly))
    
    rasterLayers.i <- subset(rasterLayers,which(apply(raster::extract(rasterLayers,occurrenceRecords),2,function(x) { length(unique(range(x))) > 1  })))
    
    PCA_Eu <- RStoolbox::rasterPCA(img = rasterLayers, nSamples = NULL, nComp = 2, spca = T)
    PCstack <- raster::stack(PCA_Eu$map$PC1, PCA_Eu$map$PC2)
    PCstack.df <- as.data.frame(PCstack, xy = T, na.rm = T)
    PCstack.sp <- st_as_sf(PCstack.df, coords = c("PC1", "PC2"))
    Optres_Eu <- optimRes(sdf = PCstack.sp, grid.res = seq(1, 15), cr = 1, showOpt = T)
    Optres_Eu <- ifelse(!is.na(Optres_Eu$Opt_res),Optres_Eu$Opt_res,10)
    
    paThreshold <- 0.5
    
    pseudoAbsences <- paSampling(
      rasterLayers.i,
      pres = sink.points.poly,
      thres = paThreshold,
      grid.res = Optres_Eu, # Optres_Eu$Opt_res
      sub.ts = FALSE,
      prev = prevalence,
      plot_proc = FALSE,
      verbose = FALSE
    )

    pseudoAbsences <- data.frame(x=pseudoAbsences$x,y=pseudoAbsences$y)

  }
  
  if( paType == "mahalanobis" ) {
    
    solution <- function(v, len) {
      if (len<=1) {
        as.list(v)
      } else {
        append_each_to_list(solution(v, len-1), v)
      } 
    }
    append_each_to_list <- function(L, v) {
      purrr::flatten(lapply(v, 
                            function(n) lapply(L, function(l) c(l, n))
      ))
    }
    
    mahalanobis <- NULL
    library(adehabitatHS)
    rasterLayersSDF <- as.data.frame(rasterLayers, na.rm=T)
    rasterLayersSDFComb <- solution(c(1,0), ncol(rasterLayersSDF))
    comb.i <- 0
    
    while(is.null(mahalanobis))  {
      
      comb.i <- comb.i + 1
      
      rasterLayersSDF <- as.data.frame(rasterLayers, na.rm=T, xy=T)
      rasterLayersSDF <- rasterLayersSDF[complete.cases(rasterLayersSDF),]
      rasterLayersSDF.pts <- rasterLayersSDF[,c(1,2)]
      rasterLayersSDF <- rasterLayersSDF[,-c(1,2)]
      toKeep.col <- which(rasterLayersSDFComb[[comb.i]] == 1)
      
      if(length(toKeep.col) == 0) { break }
      
      rasterLayersSDF <- rasterLayersSDF[,toKeep.col]
      rasterLayersSDF <- SpatialPixelsDataFrame(rasterLayersSDF.pts, as.data.frame(rasterLayersSDF))
      
      tryCatch( 
        mahalanobis <- raster(x=mahasuhab(rasterLayersSDF,SpatialPointsDataFrame(occurrenceRecords,data=data.frame(data=rep(1,nrow(occurrenceRecords)))), type="probability"))
        , error=function(e) { Error <<- TRUE })
      
      if( ! is.null(mahalanobis) ) { break }
      if( length(rasterLayersSDFComb) == comb.i ) { break }
      
    }
  
    mahalanobisThresh <- raster::extract(aggregate(mahalanobis,8,mean),occurrenceRecords)
    mahalanobisThresh <- mahalanobisThresh[complete.cases(mahalanobisThresh)]
    mahalanobisThresh <- quantile(mahalanobisThresh, probs=0.025)
    mahalanobis[mahalanobis > mahalanobisThresh ] <- NA
    pseudoAbsences <- as.data.frame(mahalanobis,xy=T,na.rm=T)[,1:2]
    
    occurrenceRecords.sp <- as.data.frame(occurrenceRecords)
    coordinates(occurrenceRecords.sp) <- ~Lon+Lat
    crs(occurrenceRecords.sp) <- crs(rasterLayers)
    occurrenceRecords.sp <- raster::buffer(occurrenceRecords.sp, width=111000)
    
    if(paRemoveOverOccurrence) {
      pseudoAbsences.sp <- pseudoAbsences
      coordinates(pseudoAbsences.sp) <- ~x+y
      crs(pseudoAbsences.sp) <- crs(rasterLayers)
      pseudoAbsences <- pseudoAbsences[which(is.na(over(pseudoAbsences.sp,occurrenceRecords.sp))),]
    }
    
  }
  
  toRemove <- which(is.na(raster::extract(rasterLayers,pseudoAbsences)),arr.ind = TRUE)[,1]
  if(length(toRemove) > 0) { pseudoAbsences <- pseudoAbsences[-toRemove,] }
  
  plot(pseudoAbsences, col="gray")
  points(occurrenceRecords, col="black")
  colnames(pseudoAbsences) <- colnames(occurrenceRecords)
  return(list(records=pseudoAbsences,biasSurface=biasSurface))
  
  options(warn=0)
  
}

## -----------------------

resampleRecords <- function(records,rasterLayers,nRecords,EnvironmentStrat=TRUE,biasSurface=NULL) {
  
  if( !is.null(biasSurface) ) {
    
    sampleRec <- sample(1:nrow(records),size=min(nRecords,nrow(records)), replace = FALSE, prob =  raster::extract(biasSurface,records))
    records <- records[sampleRec,]
    
  }
  
  if( EnvironmentStrat ) {
    
    pa.environment <- raster::extract( rasterLayers, records)
    pa.environment <- pa.environment[,which(!is.na(pa.environment[1,]))]
    pa.environment <- as.data.frame(scale(pa.environment))
    
    pa.environment <- pa.environment[,which(!is.na(apply(pa.environment,2,var)))]
    
    nRecords <- min(nRecords,nrow(pa.environment)) -1
    
    pa.environment <- kmeans(pa.environment,centers=nRecords,iter.max = 10, nstart = 1)$cluster
    pa.clusters <- unique( pa.environment )
    finalKMClusters <- sapply(pa.clusters,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
    records <- records[finalKMClusters,]
    
  }

  return(records)

}

## -----------------------

predictDistribution <- function(model,rasterLayers) {
  
  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1
  cells <- Which(!is.na(shape), cell=TRUE)
  prediction <- matrix(NA,ncol=length(model$models),nrow=length(cells))
  
  for( m in 1:length(model$models) ) {

    matrixLayers <- as.data.frame(rasterLayers[cells])
    
    if( class(model$models[[m]]) == "xgb.Booster") {
      matrixLayers = xgb.DMatrix(data = data.matrix(matrixLayers[model$models[[m]]$feature_names]), label = rep(1,nrow(matrixLayers)))
    }
    
    prediction[,m] <- predict(model$models[[m]],matrixLayers, type="response")
  }

  shape[] <- NA
  rast <- shape
  rast[cells] <- apply(prediction,1,mean)
  rast[rast < 0] <- 0
  rast[rast > 1] <- 1

  rast.sd <- shape
  rast.sd[cells] <- apply(prediction,1,sd)
  rast.sd[rast.sd < 0] <- 0
  rast.sd[rast.sd > 1] <- 1
  
  res <- list(prediction=rast,prediciton.sd=rast.sd)
  return(res)
  
}

## -----------------------

SDMmodelPlot <- function(model,rasters,variable.to.plot,modelData,extrapolate,plotExtremeThreshold,plotRateChangeThreshold,yLimits) {
  
  if( ! is.null(yLimits) ) { yLimits <- round(yLimits + c(-0.055,0.055), digits = 1) }
  if( is.null(yLimits) ) { yLimits <- c(0,1) }
  
  yLimits[yLimits < 0] <- 0
  yLimits[yLimits > 1] <- 1
  
  options(warn=-1)
  
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  ## -----------------------
  
  data.to.plot.i <- raster::extract(subset(rasters,which(names(rasters) == variable.to.plot)),modelData@coords)
  data.to.plot <- raster::extract(rasters,modelData@coords)
  data.to.plot <- apply(data.to.plot,2,mean)
  data.to.plot <- do.call("rbind", replicate(100, data.to.plot, simplify = FALSE))
  data.to.plot[,variable.to.plot] <- seq(min(data.to.plot.i),max(data.to.plot.i),length.out=100)
  data.to.plot <- as.data.frame(data.to.plot)
  
  matrixEffect <- NULL
  
  for( m in 1:length(model$models)) {

    model.i <- model$models[[m]]
    
    if ( class(model.i) == "xgb.Booster") {

      matrixEffect.m <- data.to.plot
      matrixEffect.m <- matrixEffect.m[model.i$feature_names]
      matrixEffect.m <- xgb.DMatrix(data = data.matrix(matrixEffect.m), label = rep(0,nrow(matrixEffect.m)))
      matrixEffect.m <- data.frame(Effect= predict(model.i,matrixEffect.m, type = "response") ) 
      
    }
    
    if ( class(model.i) != "xgb.Booster") {
      
      matrixEffect.m <- data.frame(Effect= predict(model.i,data.to.plot)) 
      
    }
      
    matrixEffect.m <- logit2prob(matrixEffect.m)
      
    if(!is.null(matrixEffect)) { matrixEffect <- cbind(matrixEffect,data.frame(m=matrixEffect.m)) }
    if(is.null(matrixEffect)) { matrixEffect <- data.frame(m=matrixEffect.m) }
      
    }
    
  data.to.plot <- data.frame(Variable=data.to.plot[,variable.to.plot],value=apply(matrixEffect,1,mean,na.rm=T))
  data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
  
  if(data.to.plot[nrow(data.to.plot),2] > data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] < data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }
  if(data.to.plot[nrow(data.to.plot),2] < data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] > data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }
   
  ## -----------------------
  
  speciesDataUse <- raster::extract(subset(rasters,which(names(rasters) == variable.to.plot)),modelData@coords[modelData@pa==1,])
  speciesDataUse <- data.frame(x=speciesDataUse,y=min(data.to.plot$value))
  speciesDataUse <- speciesDataUse[sort(speciesDataUse[,1], index.return =T)$ix,]
  
  ## -----------------------
  
  if( ! extrapolate ) { data.to.plot <- data.to.plot[data.to.plot$Variable >= min(speciesDataUse$x, na.rm=T) & data.to.plot$Variable <= max(speciesDataUse$x, na.rm=T) , ]  }
  
  ## -----------------------
  
  partialPlot <- ggplot() +
    geom_point( data=speciesDataUse,aes(x=x, y=y), size=2,shape=15, fill="gray", color="gray") +
    geom_line( data=data.to.plot,aes(x=Variable, y=value),color="Black", size=0.5) +
    themePlot +
    xlab(variable.to.plot) + ylab("Effect on response (probability)")
  
  ## -----------------------
  
  span.vals <- data.to.plot$value
  span.vals.var <- numeric(length(span.vals)-1)
  for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- (abs(span.vals[i+1]) - abs(span.vals[i])) / (max(span.vals)-min(span.vals)) } 
  span.vals.var <- abs(span.vals.var)
  rateChangeThreshold <- which.max(span.vals.var)
  rateChangeThreshold <- data.to.plot[rateChangeThreshold,1]
  
  ## --
  
  if( data.to.plot[1,2] > data.to.plot[nrow(data.to.plot),2]  ) { tail <- "Max" } else { tail <- "Min" }
  
  if( tail == "Min") {
    
    usedThreshold <- min(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.05,na.rm=T)
    extremeThreshold <- max(data.to.plot[data.to.plot$value <= min(data.to.plot$value) + quantile(data.to.plot$value,0.05),"Variable"])
    if(extremeThreshold == Inf | extremeThreshold == -Inf) { extremeThreshold <- max(data.to.plot[data.to.plot$value <= min(data.to.plot$value),"Variable"]) }
    
  }
  
  if( tail == "Max") {
    
    usedThreshold <- max(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.95,na.rm=T)
    extremeThreshold <- min(data.to.plot[data.to.plot$value <= min(data.to.plot$value) - quantile(data.to.plot$value,0.05) ,"Variable"])
    if(extremeThreshold == Inf | extremeThreshold == -Inf) { extremeThreshold <- min(data.to.plot[data.to.plot$value <= min(data.to.plot$value) ,"Variable"]) }
    
  }
  
  if( data.to.plot[1,2] == data.to.plot[nrow(data.to.plot),2]  ) { rateChangeThreshold <- NA }
  
  if( plotExtremeThreshold) {
    partialPlot <- partialPlot + geom_vline(xintercept = extremeThreshold,linetype = "dashed")
  }
  if( plotRateChangeThreshold ) {
    partialPlot <- partialPlot + geom_vline(xintercept = rateChangeThreshold,linetype = "dotted")
  }
  
  if(length(rateChangeThreshold) == 0) { rateChangeThreshold <- NA }
  if(length(extremeThreshold) == 0) { extremeThreshold <- NA }
  if(length(usedThreshold) == 0) { usedThreshold <- NA }
  if(length(usedThresholdQ95) == 0) { usedThresholdQ95 <- NA }
  
  tippingPoints <- data.frame( maxRateChangeThreshold=rateChangeThreshold,extremeThreshold=extremeThreshold,usedThreshold=usedThreshold,usedThresholdQ95=usedThresholdQ95 )
  
  return(list(partialPlot=partialPlot,tippingPoints=tippingPoints, partialPlotData=data.to.plot))
  
}

## -----------------------

correctLayer <- function(rasters,name,type,value.is,value.will) {		
  
  name.i <- which(names(rasters) == name)
  temp.r <- rasters[[name.i]]		
  rasters <- dropLayer(rasters,name.i)
  
  if(type =="he") { temp.r[temp.r >= value.is] <- value.will }		
  if(type =="le") { temp.r[temp.r <= value.is] <- value.will }		
  if(type =="h") { temp.r[temp.r > value.is] <- value.will }		
  if(type =="l") { temp.r[temp.r < value.is] <- value.will }		
  
  names(temp.r) <- name
  rasters <- stack( rasters , temp.r )
  
  return(rasters)		
  
}		

## -----------------------

reclassifyPredicted <- function(prediction,reclassThreshold) {

    prediction[prediction > reclassThreshold] <- 1
    prediction[prediction <= reclassThreshold] <- 0
    
  return(prediction)
  
}

## -----------------------

getBackgroundExtent <- function(occurrenceRecords,rasterLayers,regionBuffer, bufferStep=1  ) {
  
  extentOccurrenceRecords <- c(min(occurrenceRecords[,1]),max(occurrenceRecords[,1]),min(occurrenceRecords[,2]),max(occurrenceRecords[,2]))
  
  backgroundInformation <- subset(rasterLayers,1)
  backgroundInformation[!is.na(backgroundInformation)] <- 1

  sink.points.poly <- as.data.frame(occurrenceRecords)
  coordinates( sink.points.poly ) <- dataRecordsNames
  proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
  sink.points.poly.list <- list()
  for( l in 1:length(sink.points.poly) ) {
    sink.points.poly.list <- c(sink.points.poly.list, list(gBuffer( sink.points.poly[l,], width=1 )))
  }
  sink.points.poly <- do.call("rbind",c(args = sink.points.poly.list, makeUniqueIDs = TRUE))
  row.names(sink.points.poly) <- as.character(1:length(sink.points.poly))
  sink.points.poly <- SpatialPolygonsDataFrame(sink.points.poly, data.frame(val=as.character(1:length(sink.points.poly))))
  sink.points.poly$val <- 1
  sink.points.poly.base <- gUnaryUnion(sink.points.poly, id = sink.points.poly$val)
  
  results <- data.frame()
  
  for(buff in seq(2,regionBuffer+1,by=bufferStep)) {
    
    sink.points.poly <- gBuffer( sink.points.poly.base, width=buff )
    sink.points.poly <- rasterize(sink.points.poly, backgroundInformation)
    sink.points.poly <- mask(sink.points.poly,backgroundInformation)
    backgroundInformation.i <- sample(Which(sink.points.poly == 1, cell=TRUE),min(length(Which(sink.points.poly == 1, cell=TRUE)),10000),replace=FALSE)
    backgroundInformation.i.coords <- xyFromCell(sink.points.poly,backgroundInformation.i)
      
    pa.environment <- rasterLayers[backgroundInformation.i]
    backgroundInformation.i.coords <- backgroundInformation.i.coords[which(complete.cases(pa.environment)),]
    
    pa.environment <- pa.environment[which(complete.cases(pa.environment)),]
    pa.environment <- pa.environment[,apply(pa.environment,2,var, na.rm=T) > 0]
    
    final.paRatio.i <- min(paMinimum,nrow(pa.environment)-1)
    
    pa.environment <- kmeans(pa.environment,centers=final.paRatio.i,iter.max = 10, nstart = 1)$cluster
    pa.clusters <- unique( pa.environment )
    finalKMClusters <- sapply(pa.clusters,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
    absences.i <- backgroundInformation.i.coords[finalKMClusters,]
    
    dataModel <- data.frame(occ=1,raster::extract(rasterLayers,occurrenceRecords[,1:2]))
    dataModel <- rbind(dataModel,data.frame(occ=0,raster::extract(rasterLayers,absences.i)))
    dataModel <- dataModel[complete.cases(dataModel),]
    dataModel <- dataModel[,c(1, which(apply(dataModel[,-1],2,var, na.rm=T) > 0 ))]
    
    predicted <- predict( glm(dataModel, family=binomial(link = "logit")),dataModel,type="response")
    predicted.accuracy <- modelPerformance(dataModel$occ,predicted,index="auc")
    
    results <- rbind(results,data.frame(buff=buff,performance=predicted.accuracy$auc))
    
  }
  
  ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2) }
  results[,2] <- as.numeric(ma(results[,2]))
  results <- results[complete.cases(results),]
  colnames(results) <- c("x","y")

  bestResult <- results[which.max(results$y),"x"]

  if( sum(results$y == max(results$y)) == length(results$y) ) { bestResult <- round(max(results[,"x"]) - 5)  }
  
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("Best buffer size distance: ",bestResult," degrees"))
  
  figure <- ggplot() +
    geom_vline(xintercept=bestResult, size=0.1) +
    geom_line(data = results, aes(x = x, y = y), size = 0.1, color="Black", linetype = "dashed") +
    geom_point(shape=19, data = results, aes(x = x, y = y), size = 1.75, color="Black") +
    ylab("Performance of preliminary models") + xlab("Buffer size (degree)") + themePlot
  
  return( list(figure=figure, bestResult = bestResult) )
  
}

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------












