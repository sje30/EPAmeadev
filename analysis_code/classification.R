#Performs classification published in "Characterization of Early Cortical Neural Network Development in Multiwell Microelectrode Array Plates". Functions contained in classification_functions.R
library(pracma)
library(randomForest)
library(e1071)
library(doBy)

set.seed(123)

load("data_all.RData")
source("classification_functions.R")

ages<-c(5,7,9,12)
feature.mat<-data.all[,-c(13,14)]#Create feature matrix containing 12 features and DIV
feature.mat$DIV<-as.factor(feature.mat$DIV)
feature.mat.zeroed<-feature.mat
for (i in c(1:6, 9:11)){
  feature.mat.zeroed[is.na(feature.mat.zeroed[,i]),i]<-0
}
feature.mat.nona<-feature.mat.zeroed[-which(apply(feature.mat.zeroed, 1, function(x) any(is.na(x)))),]#Feature matrix with no NA values

classification<-NULL
#Classification using random forest. Output is average accuracy and importance of each factor
classification$random.forest<-boost.tree.age(feature.mat.nona) 

names(classification$random.forest$factors)<-c("Firing rate", "Within burst firing rate", "Bursting electrodes", "Burst rate", "Burst duration",  "Spikes in bursts", "CV of IBI" ,"CV of within burst ISI", "Network spike rate", "Network spike peak", "Network spike duration", "Correlation")

imp.order<-order(classification$random.forest$factors)  #Order of importance of variables
names.order<-names(classification$random.forest$factors[rev(imp.order)])

tune.out<-tune(svm, DIV~. , data=feature.mat.nona, kernel="radial", ranges = list(cost=c(0.001,0.01,0.1,1,5,10,100) , gamma=logseq(0.001, 10, 5))) #Find optimal cost and gamma values for SVM classification
best.cost<-tune.out$best.parameters$cost
best.gamma<-tune.out$best.parameters$gamma

#Classification using SVMs. Features are removed in order of their decreasing importance.
classification$svm.all<-data.frame(feature=names.order, accuracy=svm.class(feature.mat.nona, rem.vec=imp.order)) 

DIV.pair<-combn(ages, 2)
svm.pairs<-data.frame(feature=names.order)
#Classification using SVMs on wach pairwise combination of ages. 
for (i in 1:dim(DIV.pair)[2]) {
  DIVs<-DIV.pair[,i]
  svm.pairs<-cbind(svm.pairs, svm.class(feature.mat.nona, rem.vec=imp.order, DIVs=DIVs))
}
names(svm.pairs)<-c("feature", apply(DIV.pair, 2, function(x) paste("DIV", x[1], x[2], sep="")))
classification$svm.pairs<-svm.pairs

#Classification using only a fraction of the 48 wells
plate.feature.mat<-cbind(feature.mat.zeroed, Plate.id=data.all[,"Plate.id"])
plate.ids<-unique(plate.feature.mat$Plate.id)
SVM.wells<-lapply(c(1:5,seq(6,48,2)), function(x) frac.pred(plate.feature.mat, x)) #Classification using between 1 and 48 wells per plate
classification$SVM.frac.wells<-do.call("rbind", SVM.wells)



