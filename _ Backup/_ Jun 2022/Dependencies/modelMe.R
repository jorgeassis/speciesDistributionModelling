## -----------------------
closeAllConnections()
gc(reset=TRUE)
## -----------------------

sorensenIndex <- function(obs,pred) {
  
  tp <- sum( obs == 1 &  pred == 1 ) / sum( obs == 1 )
  fp <- sum( obs == 0 &  pred == 1 ) / sum( obs == 0 )
  fn <- sum( obs == 1 &  pred == 0 ) / sum( obs == 1 )
  
  sorensen.i <- (2 * tp) / ( fn + (2 * tp) + fp )
  return(sorensen.i)
  
}

## -----------------------

accuracy.area.c <- function(model,observed,predicted,predicted.map) {
  
  options(warn=-1)
  
  predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )
  predicted.accuracy$auc <- SDMTools::auc(observed ,predicted)
  
  predicted.distribution.area <- NA
  model.deviance <- Dsquared(glm(predicted~observed,family=binomial()))
  aicc <- NA
  
  if( class(model)[1] == "MaxEnt") { n.param <- get.params(model) ; aicc <- calc.aicc(n.param, presences.test, predicted.map)$AICc }
  
  if( class(model)[1] == "gbm") { model.deviance <- 1 - model$cv.statistics$deviance.mean }
  
  if(cvIndex == "tss" ) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
  
  if(cvIndex == "deviance") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$deviance),] }
  
  if(cvIndex == "aicc") { predicted.accuracy <- predicted.accuracy[which.min(predicted.accuracy$aicc),] }
  
  if(cvIndex == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
  
  if(cvIndex == "area") { 
    
    predicted.accuracy <- predicted.accuracy[predicted.accuracy$sensitivity >= 0.95,] 
    predicted.accuracy <- predicted.accuracy[nrow(predicted.accuracy),] 
    
    if( nrow(predicted.accuracy) == 0 ) { predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )[1,] }
    
    predicted.distribution <- crop(predicted.map,extent( min( occurrenceRecords[,"Lon"] ) , max( occurrenceRecords[,"Lon"] ) , min( occurrenceRecords[,"Lat"] ) , max( occurrenceRecords[,"Lat"] )))
    
    full.area <- predicted.distribution
    full.area[!is.na(full.area)] <- 1
    full.area <- raster::resample(r.area,full.area)
    full.area <- raster::mask(full.area,predicted.distribution)
    full.area <- sum(getValues(full.area),na.rm=T)
    
    predicted.distribution[ predicted.distribution > predicted.accuracy$threshold ] <- 1
    predicted.distribution[ predicted.distribution <= predicted.accuracy$threshold ] <- NA
    predicted.distribution.area <- raster::resample(r.area,predicted.distribution)
    predicted.distribution.area <- raster::mask(predicted.distribution.area,predicted.distribution)
    predicted.distribution.area <- 1 - (sum(getValues(predicted.distribution.area),na.rm=T) / full.area)
    
  }
  
  threshold <- predicted.accuracy$threshold 
  auc <- predicted.accuracy$AUC 
  specificity <- predicted.accuracy$specificity 
  sensitivity <- predicted.accuracy$sensitivity 
  tss <- predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 
  sorensen <- sorensenIndex( observed, ifelse(predicted >= threshold , 1 , 0) )
  
  if( length(threshold) == 0 ) { threshold <- NA }
  if( length(auc) == 0 ) { auc <- NA }
  if( length(specificity) == 0 ) { specificity <- NA }
  if( length(sensitivity) == 0 ) { sensitivity <- NA }
  if( length(tss) == 0 ) { tss <- NA }
  if( length(sorensen) == 0 ) { sorensen <- NA }
  
  predicted.accuracy <- data.frame( sorensen=sorensen,
                                    threshold = threshold ,
                                    auc = auc ,
                                    specificity = specificity ,
                                    sensitivity = sensitivity ,
                                    tss = tss ,
                                    area = predicted.distribution.area ,
                                    aicc = aicc,
                                    deviance = model.deviance )
  
  options(warn=0)
  
  return(predicted.accuracy)
  
}

## -----------------------

cv.k.span <- 1:cvKFolds
predictors <- names(rasterLayers)

## -----------------------

if( cvIndex == "area" ) { r.area <- raster::area(subset(rasterLayers,1)) }

## -----------------------------------------------------------------------------------

cross.validation.rounds <- list()

dataLayersMonotonocity.i <- sapply(predictors,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )

## ----------------

if( modelType == "MaxEnt" ) {
  
  stop("Error 02: Check Section")
  
  comb = expand.grid(cv.k = cv.k.span, feature.class=sapply(1:length(maxent.feature.class.span),function(x) { paste0(t(maxent.feature.class.span[x]),collapse=",")} ),betamultiplier=maxent.betamultiplier.span )
  
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------')
  cat('\n')
  cat('\n')
  cat('Model: MaxEnt \n')
  cat('Cross-validation rounds:', nrow(comb))
  cat('\n')
  cat('Number of parameters:', length(maxent.feature.class.span), '+' , length(maxent.betamultiplier.span))
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------')
  
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  cv.accuracy <- foreach(c = 1:nrow(comb), .packages = c("dismo", "rJava","SDMTools","ENMeval")) %dopar% {
    
    cv <- comb[c,1]
    feature <- as.character(unlist(strsplit(as.character(comb[c,2]),split=",")))
    beta <- comb[c,3]
    args <- c(paste0("betamultiplier=",beta),feature)
    
    # ----------------------
    
    if( TRUE %in% grepl("lat", colnames(crossValidation)) ) {
      
      presences.train <- presences[ presences$Lat <= crossValidation[cv,1] | presences$Lat >= crossValidation[cv,2] , ]
      absences.train <- absences[ absences$Lat <= crossValidation[cv,1] | absences$Lat >= crossValidation[cv,2] , ]
      presences.test <- presences[ presences$Lat >= crossValidation[cv,1] & presences$Lat <= crossValidation[cv,2] , ]
      absences.test <- absences[ absences$Lat >= crossValidation[cv,1] & absences$Lat <= crossValidation[cv,2] , ]
      
    }
    
    if( TRUE %in% grepl("lon", colnames(crossValidation)) ) {
      
      presences.train <- presences[ presences$Lon <= crossValidation[cv,1] | presences$Lon >= crossValidation[cv,2] , ]
      absences.train <- absences[ absences$Lon <= crossValidation[cv,1] | absences$Lon >= crossValidation[cv,2] , ]
      presences.test <- presences[ presences$Lon >= crossValidation[cv,1] & presences$Lon <= crossValidation[cv,2] , ]
      absences.test <- absences[ absences$Lon >= crossValidation[cv,1] & absences$Lon <= crossValidation[cv,2] , ]
      
    }  
    
    presences.train <- presences.train[!is.na(presences.train[,1]),]
    absences.train <- absences.train[!is.na(absences.train[,1]),]
    
    presences.test <- presences.test[!is.na(presences.test[,1]),]
    absences.test <- absences.test[!is.na(absences.test[,1]),]
    
    if( nrow(presences.train) > 0 &nrow(absences.train) > 0 & nrow(presences.test) > 0 & nrow(absences.test) > 0 ) { 
      
      train.dataset <- data.frame( PA = c( rep(1,nrow(presences.train)) , rep(0,nrow(absences.train)) ) , raster::extract( rasterLayers , rbind( presences.train, absences.train) ) )
      test.dataset <- data.frame( PA = c( rep(1,nrow(presences.test)) , rep(0,nrow(absences.test)) ) , raster::extract( rasterLayers , rbind( presences.test, absences.test) ) )
      
      train.dataset <- train.dataset[complete.cases(train.dataset),]
      test.dataset <- test.dataset[complete.cases(test.dataset),]
      
      model <- dismo::maxent(x = rasterLayers , p = presences.train, a=absences.train , args=c("-P","autorun=false",args) )
      
      observed <- test.dataset$PA
      predicted <- predict( model, test.dataset[,-1])
      
      if( cvIndex == "area" ) { 
        predicted.map <- predict( rasterLayers , model )
      } else { predicted.map <- NULL}
      
      predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
      
      predicted.accuracy <- data.frame( cv.round=cv,
                                        feature=paste(feature,collapse="+"),
                                        beta=beta,
                                        predicted.accuracy)
      
      return(predicted.accuracy)
      
    }
    
  }
  
  stopCluster(cl); rm(cl) ; gc(reset=TRUE)
  
  # ------------------
  
  cv.accuracy <- do.call(rbind,cv.accuracy)
  cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold),]
  
  # ------------------
  
  if(cvIndex == "aicc") { 
    best.model <- aggregate(list(aicc=cv.accuracy$aicc), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
    best.model.i <- which.min(best.model[,"aicc"])
  }
  if(cvIndex == "tss") { 
    best.model <- aggregate(list(tss=cv.accuracy$tss), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
    best.model.i <- which.max(best.model[,"tss"])
  }
  if(cvIndex == "auc") { 
    best.model <- aggregate(list(auc=cv.accuracy$auc), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
    best.model.i <- which.max(best.model[,"auc"])
  }
  if(cvIndex == "area") { 
    best.model <- aggregate(list(area=cv.accuracy$area), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
    best.model.i <- which.max(best.model[,"area"])
  }
  
  # ------------------
  
  best.model <- best.model[best.model.i,]
  
  best.model.feature <- as.character(best.model$feature)
  best.model.beta <- as.numeric(best.model$beta)
  
  best.model.metrics <- cv.accuracy[cv.accuracy$feature == best.model.feature & cv.accuracy$beta == best.model.beta ,  ]
  best.cross.validation <- best.model.metrics[,cvIndex]
  
  model <- dismo::maxent(x = rasterLayers , p = presences, a=absences , args=c("-P","autorun=false",c(paste0("betamultiplier=",best.model.beta),unlist(strsplit(best.model.feature, "[+]")))) )
  
  if( simplifyModels ) {
    
    model.importance <- summary.model(model,print.data=FALSE)
    model.simplift.vars <- as.character(model.importance$Variable)[model.importance$Permutation <= 1]
    
    model <- dismo::maxent(x = dropLayer(rasterLayers,which(names(rasterLayers) %in% model.simplift.vars)) , p = presences, a=absences , args=c("-P","autorun=false",c(paste0("betamultiplier=",best.model.beta),unlist(strsplit(best.model.feature, "[+]")))) )
    model.predicted <- predict( rasterLayers , model )
    
  }
  
}

## ------------------------------------------------------------------------------

if( modelType == "BRT" ) {
  
  if( length(names(rasterLayers)) == 1 ) {  
    
    dummy.raster <- rasterLayers
    dummy.raster[1:length(rep(1:2,length(dummy.raster) / 2))] <- rep(1:2,length(dummy.raster) / 2)
    names(dummy.raster) <- "layer"
    rasterLayers <- stack(rasterLayers,dummy.raster)
    dataLayersMonotonocity.i <- data.frame(data.frame(t(dataLayersMonotonocity.i)),layer=+1)
    predictors <- c(predictors,"layer")
    
  }
  
  comb = expand.grid(cv.k = cv.k.span, learning.complex=brtLearning,tree.depth=brtTreeDepth , bag=brtBagFraction )
  
  cl <- parallel::makeCluster(nCores)
  registerDoParallel(cl)
  
  cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind,  .export=c("brtMaxTrees","Dsquared"), .packages = c("dismo","SDMTools","ENMeval")) %dopar% {
    
    cv <- comb[c,1]
    l.rate <- comb[c,2]
    tree.c <- comb[c,3]
    bag <- comb[c,4]
    
    train.dataset <- speciesData[crossValidation[[cv]][[1]],]
    test.dataset <- speciesData[crossValidation[[cv]][[2]],]
    
    if( length(unique(train.dataset$PA)) == 2 & length(unique(test.dataset$PA)) == 2  ) { 
      
      if( (sum(train.dataset$PA == 1) / sum(train.dataset$PA == 0)) < 0.025 ) { 
        
        train.dataset <- rbind(train.dataset,do.call("rbind", replicate( 5 , train.dataset[train.dataset$PA == 1 ,], simplify = FALSE)))  
        
      }
      
      train.dataset <- data.frame( PA = train.dataset$PA , raster::extract( rasterLayers , train.dataset[,c("Lon","Lat")] ) )
      train.dataset[train.dataset == "NaN"] <- NA
      test.dataset <- data.frame( PA = test.dataset$PA , raster::extract( rasterLayers , test.dataset[,c("Lon","Lat")] ) )
      test.dataset[test.dataset == "NaN"] <- NA
      
      weights <- numeric(nrow(train.dataset))
      weights[!is.na(weights)] <- 1
      
      # Tune weights
      
      if( useWeights) {
        
        weights[which(train.dataset[,1] == 1)] <- 1
        weights[which(train.dataset[,1] == 0)] <- sum(train.dataset[,1] == 1) / sum(train.dataset[,1] == 0)
        
      }

      tryCatch( model <- gbm.step( data=train.dataset, 
                                   gbm.x = which(colnames(train.dataset) %in% predictors),
                                   gbm.y = 1, 
                                   family = "bernoulli", 
                                   plot.main = FALSE,
                                   tree.complexity = tree.c, 
                                   learning.rate = l.rate, 
                                   bag.fraction = bag, 
                                   n.folds=6,
                                   site.weights = weights,
                                   max.trees=brtMaxTrees,
                                   silent=TRUE,
                                   var.monotone = dataLayersMonotonocity.i ,
                                   verbose=FALSE) , error=function(e) { model <<- NULL })
      
      if( ! is.null(model) ) {
        
        num.tress <- model$gbm.call$best.trees
        observed <- test.dataset$PA
        predicted <- predict( model , test.dataset[,-1] , n.trees=num.tress,type="response")
        
        if( cvIndex == "area" ) { 
          predicted.map <- predict( rasterLayers , model , n.trees=num.tress,type="response")
        } else { predicted.map <- NULL}
        
        predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
        
      }
      
      if( is.null(model)) {
        
        predicted.accuracy <- data.frame( sorensen=NA, threshold=NA, auc=NA, specificity=NA, sensitivity=NA, tss=NA, area=NA, aicc=NA, deviance=NA )
        
      }
      
      predicted.accuracy <- data.frame( cv.round=cv,tree.c=tree.c,l.rate=l.rate,bag=bag,predicted.accuracy)
      return(predicted.accuracy)
      
    }
    
    if( length(unique(train.dataset$PA)) != 2 | length(unique(test.dataset$PA)) != 2 ) { return(NULL) }
    
  }
  
  stopCluster(cl); rm(cl)
  closeAllConnections()
  gc(reset=TRUE)
  
  # ------------------
  
  cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
  
  best.model <- cv.accuracy
  best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(bag=best.model$bag,tree.c=best.model$tree.c,l.rate=best.model$l.rate), mean)
  best.model.i <- best.model[sort(best.model$cvIndex, decreasing = TRUE, index.return = TRUE)$ix,]
  
  if( nrow(cv.accuracy) > 0 ) { 
    
    model <- NULL
    m.i <- 1
    
    best.model.tc <- best.model.i[m.i,"tree.c"]
    best.model.lr <- best.model.i[m.i,"l.rate"]
    best.model.bag <- best.model.i[m.i,"bag"]
    
    if( m.i == nrow(best.model.i) ) { best.model.lr <- 0.00001 ; best.model.tc <- 2 }
    
    best.model.metrics <- cv.accuracy[cv.accuracy$tree.c == best.model.tc & cv.accuracy$l.rate == best.model.lr & cv.accuracy$bag == best.model.bag ,  ]
    best.model.metrics <- best.model.metrics[best.model.metrics$threshold != 0 ,]
    best.cross.validation <- best.model.metrics[,cvIndex]
    
    train.dataset <- data.frame( PA = speciesData$PA , raster::extract( rasterLayers , speciesData[,c("Lon","Lat")] ) )
    train.dataset[train.dataset == "NaN"] <- NA
    
    if( sum(train.dataset$PA == 1) < 25 ) { train.dataset <- rbind(train.dataset,do.call("rbind", replicate( sum(train.dataset$PA == 1) * 10 , train.dataset[train.dataset$PA == 1 ,], simplify = FALSE)))  }
    
    weights <- numeric(nrow(train.dataset))
    weights[!is.na(weights)] <- 1
    
    if(useWeights) {
      weights[which(train.dataset[,1] == 1)] <- 1
      weights[which(train.dataset[,1] == 0)] <- sum(train.dataset[,1] == 1) / sum(train.dataset[,1] == 0)
    }

    tryCatch( model <- gbm.step( data=train.dataset, 
                                 gbm.x = which(colnames(train.dataset) %in% predictors),
                                 gbm.y = 1, 
                                 family = "bernoulli", 
                                 plot.main = FALSE,
                                 tree.complexity = best.model.tc, 
                                 learning.rate = best.model.lr, 
                                 bag.fraction = best.model.bag, 
                                 n.folds=10,
                                 step.size=50,
                                 site.weights = weights,
                                 max.trees=brtMaxTrees,
                                 silent=TRUE,
                                 var.monotone = dataLayersMonotonocity.i, 
                                 verbose=FALSE) , error=function(e) { Error <- TRUE })
    
    
  }
  
  if( nrow(cv.accuracy) == 0 ) {
    
    model <- NULL 
    best.cross.validation <- rep(0,10)
    best.model.metrics <- data.frame(threshold = 1 ,
                                     cv <- rep(0,10) ,
                                     auc= 0 ,
                                     specificity = 0 ,
                                     sensitivity = 0 ,
                                     tss = 0 ,
                                     aicc = 0 ,
                                     area = 1 ,
                                     deviance = 0 )
    
  }
  
}

## ------------------------------------------------------------------------------

if( modelType == "MBOOST" ) {
  
  comb = expand.grid(cv.k = cv.k.span, shrinkage=mboostShrinkage,df=mboostDF,mstop=mboostIterations )
  
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  cv.accuracy <- foreach(c = 1:nrow(comb), .combine = rbind, .export=c("brtMaxTrees","Dsquared"), .packages = c("dismo","SDMTools","ENMeval","mboost")) %dopar% {
    
    cv <- comb[c,1]
    shrinkage <- comb[c,2]
    df <- comb[c,3]
    mstop <- comb[c,4]
    
    train.dataset <- speciesData[crossValidation[[cv]][[1]],]
    test.dataset <- speciesData[crossValidation[[cv]][[2]],]
    
    if( length(unique(train.dataset$PA)) == 2 & length(unique(test.dataset$PA)) == 2  ) { 
      
      if( sum(train.dataset$PA == 1) < 10 ) { train.dataset <- rbind(train.dataset,do.call("rbind", replicate( sum(train.dataset$PA == 1) * 5 , train.dataset, simplify = FALSE)))  }
      
      train.dataset <- data.frame( PA = train.dataset$PA , raster::extract( rasterLayers , train.dataset[,c("Lon","Lat")] ) )
      train.dataset[train.dataset == "NaN"] <- NA
      test.dataset <- data.frame( PA = test.dataset$PA , raster::extract( rasterLayers , test.dataset[,c("Lon","Lat")] ) )
      test.dataset[test.dataset == "NaN"] <- NA
      
      colnames(train.dataset) <- c("PA",names(rasterLayers))
      colnames(test.dataset) <- c("PA",names(rasterLayers))
      
      train.dataset[,1] <- as.factor(train.dataset[,1])
      
      predictors <- colnames(train.dataset)[-1]
      response <- colnames(train.dataset)[1]
      
      constr_mono <- as.character(dataLayersMonotonocity.i)
      constr_mono[constr_mono == "-1"] <- "decreasing"
      constr_mono[constr_mono == "1"] <- "increasing"
      
      rhs <- paste(c(paste("bmono(", predictors, ", constraint = \"", constr_mono,"\", df = ",df,")", sep = "")),collapse = " + ")
      fm_mono <- as.formula(paste(response, " ~ ", rhs, collapse = ""))
      
      ctrl <- boost_control(mstop = mstop, trace = FALSE, nu=shrinkage,stopintern=TRUE)
      
      tryCatch( model <- mboost(fm_mono, data = train.dataset, control = ctrl, family = Binomial(type = "adaboost" , link = "logit" )) , error=function(e) { Error <<- TRUE })
      
      if( exists("model") ) { 
        
        observed <- test.dataset$PA
        predicted <- predict(model,test.dataset, type = "response")
        
        if( cvIndex == "area" ) { 
          predicted.map <- predict(rasterLayers,model)
          predicted.map <- predicted.map + ( min(getValues(predicted.map),na.rm=T) * (-1))
          predicted.map <- predicted.map / max(getValues(predicted.map),na.rm=T)
        } else { predicted.map <- NULL}
        
        predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
        
      }
      
      if( ! exists("model") ) { 
        
        predicted.accuracy <- data.frame( sorensen=NA, threshold=NA, auc=NA, specificity=NA, sensitivity=NA, tss=NA, area=NA, aicc=NA, deviance=NA )
        
      }
      
      predicted.accuracy <- data.frame( cv.round=cv,
                                        shrinkage=shrinkage,
                                        df=df,
                                        n.iterations=mstop,
                                        predicted.accuracy)
      
      return(predicted.accuracy)
      
    }
    
  }
  
  stopCluster(cl); rm(cl) ; gc(reset=TRUE)
  
  # ------------------
  
  cv.accuracy[cv.accuracy$threshold == 0 & !is.na(cv.accuracy$deviance) ,"threshold" ] <- 0.1
  best.model <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
  best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(shrinkage=best.model$shrinkage,df=best.model$df,n.iterations=best.model$n.iterations), mean)
  best.model.i <- best.model[sort(best.model$cvIndex, decreasing = TRUE, index.return = TRUE)$ix,]
  
  model <- NULL
  m.i <- 0
  
  while( is.null(model) ){
    
    m.i <- m.i + 1
    
    if(m.i > nrow(best.model.i) ) { break }
    
    best.model.shrinkage <- best.model.i[m.i,]$shrinkage
    best.model.df <- best.model.i[m.i,]$df
    best.model.iterations <- best.model.i[m.i,]$n.iterations
    
    best.model.metrics <- cv.accuracy[cv.accuracy$shrinkage == best.model.shrinkage & cv.accuracy$df == best.model.df & cv.accuracy$n.iterations == best.model.iterations ,  ]
    best.cross.validation <- best.model.metrics[,cvIndex]
    
    train.dataset <- data.frame( PA = speciesData$PA , raster::extract( rasterLayers , speciesData[,c("Lon","Lat")] ) )
    train.dataset[train.dataset == "NaN"] <- NA
    
    if( sum(train.dataset$PA == 1) < 8 ) { train.dataset <- rbind(train.dataset,do.call("rbind", replicate( sum(train.dataset$PA == 1) * 5 , train.dataset[train.dataset$PA == 1 ,], simplify = FALSE)))  }
    
    train.dataset <- train.dataset[complete.cases(train.dataset),]
    train.dataset[,1] <- as.factor(train.dataset[,1])
    
    colnames(train.dataset) <- c("PA",names(rasterLayers))
    
    predictors <- colnames(train.dataset)[-1]
    response <- colnames(train.dataset)[1]
    
    constr_mono <- as.character(dataLayersMonotonocity.i)
    constr_mono[constr_mono == "-1"] <- "decreasing"
    constr_mono[constr_mono == "1"] <- "increasing"
    
    rhs <- paste(c(paste("bmono(", predictors, ", constraint = \"", constr_mono,"\", df = ",best.model.df,")", sep = "")),collapse = " + ")
    fm_mono <- as.formula(paste(response, " ~ ", rhs, collapse = ""))
    
    ctrl <- boost_control(mstop = best.model.iterations, trace = FALSE, nu=best.model.shrinkage,stopintern=TRUE)
    
    tryCatch( model <- mboost(fm_mono, data = train.dataset, control = ctrl, family = Binomial(type = "adaboost" , link = "logit" )) , error=function(e) { Error <- TRUE })
    
  }
  
}

## ------------------------------------------------------------------------------

if( modelType == "MPNN" ) {
  
  comb = expand.grid(cv.k = cv.k.span, hidden=mpnnHidden,itereractions=mpnnItereractions )
  
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .packages = c("monmlp","dismo","SDMTools","ENMeval")) %dopar% {
    
    cv <- comb[c,1]
    hidden <- comb[c,2]
    itereractions <- comb[c,3]
    
    train.dataset <- speciesData[crossValidation[[cv]][[1]],]
    test.dataset <- speciesData[crossValidation[[cv]][[2]],]
    
    if( length(unique(train.dataset$PA)) == 2 & length(unique(test.dataset$PA)) == 2 ) { 
      
      train.dataset <- data.frame( PA = train.dataset$PA , raster::extract( rasterLayers , train.dataset[,c("Lon","Lat")] ) )
      train.dataset[train.dataset == "NaN"] <- NA
      test.dataset <- data.frame( PA = test.dataset$PA , raster::extract( rasterLayers , test.dataset[,c("Lon","Lat")] ) )
      test.dataset[test.dataset == "NaN"] <- NA
      
      correctLayers <- names(rasterLayers)[which(dataLayersMonotonocity.i == -1)]
      
      for(correctLayers.i in correctLayers) {
        train.dataset[,correctLayers.i] <- train.dataset[,correctLayers.i] * (-1)
        test.dataset[,correctLayers.i] <- test.dataset[,correctLayers.i] * (-1)
      }
      
      train.dataset <- train.dataset[complete.cases(train.dataset),]
      
      model <- monmlp.fit(x = as.matrix(train.dataset[,-1]), y = as.matrix(train.dataset[,1]),  monotone = 1:(ncol(train.dataset)-1),  hidden1 = hidden, bag = TRUE, iter.max = itereractions, iter.stopped = 10)
      
      if( ! is.null(model) ) {
        
        observed <- test.dataset$PA
        predicted <- as.numeric(monmlp.predict(as.matrix(test.dataset[,-1]), model))
        # odds <- exp(predicted)
        # predicted <- odds / (1 + odds)
        
        if( cvIndex == "area" ) { 
          predicted.map <- predictDistribution(rasterLayers,model, reclassToOne=FALSE)
        } else { predicted.map <- NULL}
        
        predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
        
      }
      
      if( is.null(model)) {
        
        predicted.accuracy <- data.frame( sorensen=0,threshold=0,auc=0,specificity=0,sensitivity=0,tss=0,area=0,aicc=+Inf,deviance=0)
        
      }
      
      predicted.accuracy <- data.frame( cv.round=cv,hidden=hidden,itereractions=itereractions,predicted.accuracy)
      return(predicted.accuracy)
      
    } else { return(NULL) }
    
  }
  
  stopCluster(cl); rm(cl) ; gc(reset=TRUE)
  
  # ------------------
  
  cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
  
  if( nrow(cv.accuracy) != 0 ) { 
    
    best.model <- cv.accuracy
    best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(hidden=best.model$hidden,itereractions=best.model$itereractions), mean)
    best.model.i <- which.max(best.model[,"cvIndex"])
    best.model <- best.model[best.model.i,]
    
    best.model.hidden <- best.model$hidden
    best.model.itereractions <- best.model$itereractions
    
    best.model.metrics <- cv.accuracy[cv.accuracy$hidden == best.model.hidden & cv.accuracy$itereractions == best.model.itereractions  ,  ]
    best.cross.validation <- best.model.metrics[,cvIndex]
    
    train.dataset <- data.frame( PA = speciesData$PA , raster::extract( rasterLayers , speciesData[,c("Lon","Lat")] ) )
    train.dataset[train.dataset == "NaN"] <- NA
    train.dataset <- train.dataset[complete.cases(train.dataset),]
    
    correctLayers <- names(rasterLayers)[which(dataLayersMonotonocity.i == -1)]
    for(correctLayers.i in correctLayers) {
      train.dataset[,correctLayers.i] <- train.dataset[,correctLayers.i] * (-1)
    }
    
    model <- monmlp.fit(x = as.matrix(train.dataset[,-1]), y = as.matrix(train.dataset[,1]),  monotone = 1:(ncol(train.dataset)-1),  hidden1 = best.model.hidden, bag = TRUE, iter.max = best.model.itereractions, iter.stopped = 10)
    
  }
  
  if( nrow(cv.accuracy) == 0 ) {
    
    model <- NULL 
    best.cross.validation <- rep(0,10)
    best.model.metrics <- data.frame(threshold = 1 ,
                                     cv <- rep(0,10) ,
                                     auc= 0 ,
                                     specificity = 0 ,
                                     sensitivity = 0 ,
                                     tss = 0 ,
                                     aicc = 0 ,
                                     area = 1 ,
                                     deviance = 0 )
    
  }
  
}

## -----------------------

modelFull <- list(model=model, 
             cv=best.cross.validation,
             auc=mean(best.model.metrics[,"auc"]),
             Specificity=mean(best.model.metrics[,"specificity"]),
             Sensitivity=mean(best.model.metrics[,"sensitivity"]),
             tss=mean(best.model.metrics[,"tss"]),
             area=mean(best.model.metrics[,"area"]),
             deviance=mean(best.model.metrics[,"deviance"]))

## -----------------------
closeAllConnections()
gc(reset=TRUE)
## -----------------------
