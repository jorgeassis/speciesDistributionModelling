

if( modelType == "XGBOOST" ) {
  relativeImportanceMatrix <- matrix(NA,ncol=length(model$models),nrow=nlayers(rasterLayers))
  for(m in 1:length(model$models)) {
    relativeImportance.m = as.data.frame(xgb.importance(model$models[[m]]$feature_names, model = model$models[[m]])[,1:2])
    relativeImportanceMatrix[,m] <- relativeImportance.m[match(names(rasterLayers),relativeImportance.m$Feature),"Gain"]
  }
  relativeImportance <- data.frame(Predictor=names(rasterLayers),relImportance=apply(relativeImportanceMatrix,1,mean,na.rm=T),relImportance.sd=apply(relativeImportanceMatrix,1,sd,na.rm=T))
  relativeImportance[is.na(relativeImportance)] <- 0
  relativeImportance$relImportance <- ( (relativeImportance$relImportance ) / sum(relativeImportance$relImportance) ) * 100
  relativeImportance$relImportance.sd <-  relativeImportance$relImportance.sd * 100
}

if( modelType == "BRT" ) {
  relativeImportanceMatrix <- matrix(NA,ncol=length(model$models),nrow=nlayers(rasterLayers))
  for(m in 1:length(model$models)) {
    relativeImportance.m <- summary(model$models[[m]])
    relativeImportanceMatrix[,m] <- relativeImportance.m[match(names(rasterLayers),relativeImportance.m$var),"rel.inf"]
  }
  relativeImportance <- data.frame(Predictor=names(rasterLayers),relImportance=apply(relativeImportanceMatrix,1,mean,na.rm=T),relImportance.sd=apply(relativeImportanceMatrix,1,sd,na.rm=T))
  relativeImportance[is.na(relativeImportance)] <- 0
}

if( modelType == "MBOOST" ) {
  relativeImportanceMatrix <- matrix(NA,ncol=length(model$models),nrow=nlayers(rasterLayers))
  for(m in 1:length(model$models)) {
    relativeImportanceMatrix[,m] <- sapply(1:nlayers(rasterLayers), function(x) { res <- varimp(model$models[[m]])[which(grepl(names(rasterLayers)[x],names(varimp(model$models[[m]]))))]; ifelse(length(res)>0,res,NA) } )
  }
  relativeImportanceMatrix[is.na(relativeImportanceMatrix)] <- 0
  relativeImportance <- data.frame(Predictor=names(rasterLayers),relImportance=apply(relativeImportanceMatrix,1,mean,na.rm=T),relImportance.sd=apply(relativeImportanceMatrix,1,sd,na.rm=T))
  relativeImportance$relImportance <- ((relativeImportance$relImportance - min(relativeImportance$relImportance) ) / sum(relativeImportance$relImportance) ) * 100
  relativeImportance$relImportance.sd <-  relativeImportance$relImportance.sd * 100
}

## ----------

relativeImportance$relImportance.sd <- relativeImportance$relImportance.sd / sqrt(cvKFolds)

relativeImportancePlot <- ggplot(data=relativeImportance[sort(relativeImportance[,2],decreasing = TRUE,index.return=TRUE)$ix,]) +
  geom_bar( aes(x= reorder(Predictor, relImportance) , y=relImportance), stat="identity", fill="#ffd45b", alpha=0.85) +
  geom_errorbar( aes(x= reorder(Predictor, relImportance), ymin=relImportance-ifelse(relImportance.sd == 0, NA, relImportance.sd), ymax=relImportance+ifelse(relImportance.sd == 0, NA, relImportance.sd)), width=.1, alpha=0.8, size=0.4, position=position_dodge(.9))  +
  coord_flip() + 
  theme_bw() +
  theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=13.5),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=13.5),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
  ) + labs(x = "Predictor") + 
  labs(y = "Relative importance (%)") + geom_hline(aes(yintercept=5 ),color="Black", linetype="dashed", size=0.3) +
  annotate("text", y = 5 + 1 , x = 1 , label = "5%" , hjust = 0) + geom_hline(aes(yintercept=0 ),color="Gray", size=0.3)

relativeContribution <- list(dataFrame=relativeImportance,plot=relativeImportancePlot)

## -----------------------