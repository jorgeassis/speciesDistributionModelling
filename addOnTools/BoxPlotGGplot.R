## --------------------------------------------------
## --------------------------------------------------

closeAllConnections()
gc(reset=TRUE)
library(gridExtra)

# -----------------------------
# Boxplot performance

data <- read.csv(paste0(stackResultsFolder,"/performance.csv"), sep=",")

data <- rbind(
  
  data.frame( Group= "Cross validation", Index = "Sensitivity" , Algorithm= "BRT" , Value = data$meanCV.sensitivity.BRT ),
  data.frame( Group= "Cross validation", Index = "AUC" , Algorithm= "BRT" , Value = data$meanCV.auc.BRT ),
  data.frame( Group= "Cross validation", Index = "Boyce" , Algorithm= "BRT" , Value = rescale(data$meanCV.boyce.BRT, to = c(0, 1), range = c(-1, 1)) ),
  
  data.frame( Group= "Cross validation", Index = "Sensitivity" , Algorithm= "AdaBoost" , Value = data$meanCV.sensitivity.MBOOST ),
  data.frame( Group= "Cross validation", Index = "AUC" , Algorithm= "AdaBoost" , Value = data$meanCV.auc.MBOOST ),
  data.frame( Group= "Cross validation", Index = "Boyce" , Algorithm= "AdaBoost" , Value = rescale(data$meanCV.boyce.MBOOST, to = c(0, 1), range = c(-1, 1)) ),
  
  data.frame( Group= "Cross validation", Index = "Sensitivity" , Algorithm= "XgBoost" , Value = data$meanCV.sensitivity.XGBOOST ),
  data.frame( Group= "Cross validation", Index = "AUC" , Algorithm= "XgBoost" , Value = data$meanCV.auc.XGBOOST ),
  data.frame( Group= "Cross validation", Index = "Boyce" , Algorithm= "XgBoost" , Value = rescale(data$meanCV.boyce.XGBOOST, to = c(0, 1), range = c(-1, 1)) ),
  
  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "zEnsemb." , Value = data$sensitivity.Ensemble ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "zEnsemb." , Value = data$auc.Ensemble ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "zEnsemb." , Value = rescale(data$boyce.Ensemble, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "zEnsemb. [Disp.]" , Value = data$sensitivity.EnsembleReachab ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "zEnsemb. [Disp.]" , Value = data$auc.EnsembleReachab ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "zEnsemb. [Disp.]" , Value = rescale(data$boyce.EnsembleReachab, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "BRT" , Value = data$sensitivity.BRT ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "BRT" , Value = data$auc.BRT ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "BRT" , Value = rescale(data$boyce.BRT, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "AdaBoost" , Value = data$sensitivity.MBOOST ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "AdaBoost" , Value = data$auc.MBOOST ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "AdaBoost" , Value = rescale(data$boyce.MBOOST, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "XgBoost" , Value = data$sensitivity.XGBOOST ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "XgBoost" , Value = data$auc.XGBOOST ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "XgBoost" , Value = rescale(data$boyce.XGBOOST, to = c(0, 1), range = c(-1, 1)) ))

data <- data[complete.cases(data),]

a <- ggplot(data[data$Group =="Cross validation",], aes(x=Value, y=Algorithm, fill=Index)) + 
      theme_bw() +

      geom_boxplot( lwd=0.25, outlier.size=0.4, outlier.colour="#414141") +
      coord_flip() +
      scale_fill_manual(name=NULL, values=c("#ffffe0","#e3a8aa","#f26885")) + # #FFD700","#FFB14E","#FA8775 #', ', '#', '#ffffe0
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.5","0.75","1"), limits=c(0, 1)) +
  xlab("Model performance") +
  ylab("Cross validation models") +
  # facet_grid(~Group, scales = "free_x", space = "free_x", switch = "x") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.placement = "outside") +
      theme(strip.background =element_rect(fill="white",colour = "white", size = 1)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  theme(axis.ticks = element_blank()) +
  theme(legend.position = c(0.2, 0.15))

a


b <- ggplot(data[data$Group =="Final prediction",], aes(x=Value, y=Algorithm, fill=Index)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(lwd=0.25, outlier.size=0.4, outlier.colour="#414141") +
  coord_flip() +
  scale_fill_manual(name=NULL, values=c("#ffffe0","#e3a8aa","#f26885")) + # #FFD700","#FFB14E","#FA8775 #', ', '#', '#ffffe0
  # facet_grid(~Group, scales = "free_x", space = "free_x", switch = "x") + 
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.5","0.75","1"), limits=c(0, 1)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.placement = "outside") +
  theme(strip.background =element_rect(fill="white",colour = "white", size = 1)) +
  theme(legend.position="none") + 
  ylab("Predictive models") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank() ) +
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  theme(axis.ticks = element_blank()) 

b

# Relative Contribution

data <- read.csv(paste0(stackResultsFolder,"/relativeContribution.csv"), sep=",")
data <- data[,which(grepl("ensemble",names(data)))]
names(data) <- gsub("\\.Ensemble","",names(data))
names(data) <- gsub("\\.ensemble","",names(data))

data <- tidyr::gather(data)
names(data) <- c("Predictor","Contribution")

data[data$Predictor == "TempMax","Predictor"] <- "Temp. Max."
data[data$Predictor == "TempMin","Predictor"] <- "Temp. Min."
data[data$Predictor == "seaIce","Predictor"] <- "Sea Ice"
data[data$Predictor == "CoastalExposure","Predictor"] <- "Exposure"

data$Predictor <- with(data, reorder(Predictor , -Contribution, median , na.rm=T))

c <- ggplot(data, aes(x=Predictor, y=Contribution, fill=Contribution)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_boxplot(width=0.5, fill="#93003a", lwd=0.25, outlier.size=0.4, outlier.colour="#414141") +
      scale_y_continuous(limits=c(0,75), breaks=c(0,5,25,50,75), labels=c("0","5","25","50","75")) +
      labs(y="Relative contribution (%)")+
  xlab("Predictor variable") +
  
      theme(legend.position="none") + geom_hline(yintercept=5,linetype=2, size=0.25)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  theme(axis.ticks = element_blank()) 

grid.arrange(a, b, c, nrow = 1, widths=c(0.4, 0.525, 0.5))

pdf(paste0(stackResultsFolder,"./performanceContributionPlot.pdf"), width=15, height=5)
grid.arrange(a, b, c, nrow = 1, widths=c(0.4, 0.525, 0.5))
dev.off()

