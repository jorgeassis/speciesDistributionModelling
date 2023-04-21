## ------------------------------------------------------------------------------

if( modelType == "BRT" ) {
  hyperParam <- list(n.trees=brtNTrees,interaction.depth = brtTreeDepth, shrinkage = brtLearning)
}

if( modelType == "MBOOST" ) {
  hyperParam <- list(df=mboostDF,shrinkage = mboostShrinkage , mstop=mboostMStop)
}

if( modelType == "XGBOOST" ) {
  hyperParam <-  list(gamma=xgboostGamma, shrinkage = xgboostShrinkage  , depth = xgboostDepth, nrounds = xgboostRounds)
}

gridMatrix <- modelGridSearch(modelType,hyperParam,modelData,cvFolds,monotonicity,cvIndex)
model <- modelTrainCV(modelType,modelData,cvFolds,monotonicity,gridMatrix,cvIndex)

## -----------------------

gc(reset=TRUE)

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------