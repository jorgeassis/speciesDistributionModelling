## ------------------------------------------------------------------------------

modelData <- prepareModelData(speciesData[speciesData$PA == 1, c("Lon","Lat")], speciesData[speciesData$PA == 0, c("Lon","Lat")], rasterLayers)
monotonicity <- monotonicityFull[names(monotonicityFull) %in% names(rasterLayers)]

if( modelType == "BRT" ) {
  

  model <- train("BRT", modelData, folds = folds)
  h <- list(n.trees=brtNTrees,interaction.depth = brtTreeDepth, shrinkage = brtLearning)
  exp1 <- gridSearch(model, hypers = h, metric = cvIndex)
  if(sum(is.na(exp1@results$test_TSS)) > 0) { exp1@results$test_TSS <- 1 }
  
  model <- train("BRT", modelData, 
                 n.trees = exp1@results$n.trees[which.max(exp1@results$test_TSS)],
                 shrinkage = exp1@results$shrinkage[which.max(exp1@results$test_TSS)],
                 bag.fraction =exp1@results$bag.fraction[which.max(exp1@results$test_TSS)],
                 interaction.depth = exp1@results$interaction.depth[which.max(exp1@results$test_TSS)], folds = folds)
  
}

if( modelType == "MBOOST" ) {
  
  model <- train("MBOOST", modelData, folds = folds)
  h <- list(df=mboostDF,shrinkage = mboostShrinkage , mstop=mboostMStop)
  exp1 <- gridSearch(model, hypers = h, metric = cvIndex)

  model <- train("MBOOST", modelData, 
                 df = exp1@results$df[which.max(exp1@results$test_TSS)],
                 shrinkage = exp1@results$shrinkage[which.max(exp1@results$test_TSS)],
                 mstop = exp1@results$mstop[which.max(exp1@results$test_TSS)], folds = folds)
  
}

if( modelType == "XGBOOST" ) {
  
  model <- train("XGBOOST", modelData, folds = folds) 
  h <- list(gamma=xgboostGamma, shrinkage = xgboostShrinkage  , depth = xgboostDepth, nrounds = xgboostRounds)
  exp1 <- gridSearch(model, hypers = h, metric = cvIndex)
  
  model <- train("XGBOOST", modelData, 
                 gamma = exp1@results$gamma[which.max(exp1@results$test_TSS)],
                 shrinkage = exp1@results$shrinkage[which.max(exp1@results$test_TSS)],
                 nrounds = exp1@results$nrounds[which.max(exp1@results$test_TSS)],
                 depth = exp1@results$depth[which.max(exp1@results$test_TSS)], folds = folds)
  
}

prediction <- predict(model, rasterLayers)
prediction[prediction > 1] <- 1
prediction[prediction < 0] <- 0
accur <- accuracyPredicted( prediction ,speciesData ,reclassificationIndex)
accur <- data.frame(accur,cvIndex=cvIndex, cvIndexValue=max(exp1@results[,which(grepl("test",names(exp1@results[which.max(exp1@results$test_TSS),])))]))

## -----------------------

modelFull <- list(type=modelType, model=model, performance=accur)

## -----------------------

gc(reset=TRUE)

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------