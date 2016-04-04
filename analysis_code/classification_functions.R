#Functions for classification.R

#Function that performs random forest classification. Input is data frame of feature values, output is average prediction accuracy and importance of each feature.
boost.tree.age<-function(data.df) {
  fac.num<-dim(data.df)[2]-1
  ret<-NULL
  n.row<-nrow(data.df)
  #Set size of training set to 2/3 of total data
  Ntrain<-round(n.row*2/3)
  correct<-0
  imp<-rep(0,fac.num)
  #Perform 100 trials with different training sets 
  for (i in 1:100) {
    traindat<-sample(1:nrow(data.df), Ntrain)
    testdat<-data.df[-traindat ,]
    #Calculate boosted classification tree, using all 11 features and 500 trees
    tree.out<-randomForest(DIV~. ,data=data.df[traindat, ], mtry=fac.num,importance =TRUE, ntree=500, na.action=na.omit)
    tree.pred<-predict(tree.out,testdat,type="class")
    correct<- sum(tree.pred==testdat$DIV, na.rm=TRUE)/(n.row-Ntrain)+ correct
    imp<-importance(tree.out)[,"MeanDecreaseGini"]+imp
  }
  avg.imp<-imp/100
  ret$accuracy<-correct/100
  ret$factors<-avg.imp
  ret
}


#Function to perform classification using SVMs. Input is feature matrix, optional vector of features to not include in classification and optional vector of DIVs to include in classification. Output is average accuracy.
svm.class<-function(fmat, rem.vec=NULL, DIVs=NULL) {
  if (!is.null(DIVs)){
    fmat<-subset(fmat, DIV %in% DIVs)
  }
  acc<-svm.calc(fmat, "radial", "test.set")
  if (!is.null(rem.vec)){
    N<-length(rem.vec)-1
    acc<-c(sapply(N:1, function(x) svm.calc(fmat, "radial", "test.set", rem.vec[1:x])), acc)
  }
  acc
}


#Performs SVM classification
svm.calc<-function(data.df, kernel, test.type, rm.fac=NULL) {
  n.row<-nrow(data.df)
  accuracy<-NULL
  res.df<-NULL
  #Remove features included in the rm.fac vector
  if (!is.null(rm.fac)) {
    data.df<-data.df[,-1*rm.fac]
  }
  #Set cost and gamma to optimal cost from full model
  cst<-best.cost
  gamma<-best.gamma
  
  #Test model by using 2/3 of data as training set, and remaining as test set
  if (test.type=="test.set") {
    #Set training set to be 2/3 of all data
    Ntrain<-round(n.row*2/3)
    #Determing accuracy from 100 trials with random test sets
    for (i in 1:100) {
      trainnum<-sample(1:nrow(data.df), Ntrain)
      traindat<-data.df[trainnum, ]
      testdat<-data.df[-trainnum ,]
      accuracy[i]<-svm.test(traindat, testdat, kernel, cst, gamma)
    }
  } else if (test.type=="LOOCV")  {     #Test using leave one out cross validation
    for (i in 1:n.row) {
      testdat<-data.df[i,]
      traindat<-data.df[-i,]
      accuracy[i]<-svm.test(traindat, testdat, kernel, cst, gamma)
    }
  }
  mean(accuracy)
}


#Fits and tests the accuracy of each SVM classification
svm.test<-function(traindat, testdat, kernel, best.cost, best.gamma) {
  if (kernel=="radial") {
    svmfit<-svm(DIV~. , data=traindat, kernel="radial",gamma=best.gamma, cost=best.cost, na.action=na.omit)
  } else if (kernel=="linear") {
    svmfit<-svm(DIV~. , data=traindat, kernel="linear", cost=cst)
  } else if (kernel=="poly") {
    svmfit<-svm(DIV~. , data=traindat, kernel="polynomial", degree=2, cost=cst)
  }
  #Measure accuracy as proportion of test set that are correctly predicted
  test.pred<-predict(svmfit,testdat, type="class")
  na.row<-which(apply(testdat, 1, function(x) sum(is.na(x)))>0)
  if (length(na.row)>0){
    accuracy<- sum(test.pred==testdat$DIV[-na.row])/dim(testdat)[1]
  } else {
    accuracy<- sum(test.pred==testdat$DIV)/dim(testdat)[1]
  }
}

#SVM classification using a fraction of wells on each plate. Input is feature matrix and number of wells,
#output is average accuracy over 100 trials.
frac.pred<-function(fmat, num.well){
  all.acc<-NULL
  for (j in 1:100){
    samp.df<-NULL
    for (i in ages){
      samp.list<-lapply(plate.ids, function(x) find.data.sub(fmat,i, x, num.well))
      samp.df<-rbind(samp.df, do.call("rbind", samp.list))
    }
    na.rows<-union(which(is.na(samp.df$CV.of.IBI)), union(which(is.na(samp.df$CV.of.within.burst.ISI)), which(is.na(samp.df$Correlation))))
    acc<-svm.calc(samp.df[-na.rows,-14], "radial", "test.set", rm.fac=NULL)
    all.acc<-rbind(all.acc, data.frame(num.well=num.well, accuracy=acc))
  }
  all.acc
}

#Choose a subset of wells from particular DIV and plate. Input is feature mat, DIV, plate ID and number of wells. Output is data from chosen wells.
find.data.sub<-function(fmat, DIV.val, plate.id.val, num.well){
  sub.df<-subset(fmat,  Plate.id==plate.id.val & DIV==DIV.val)
  N<-dim(sub.df)[1]
  if(N>0){
    data.subset<-sub.df[sample(nrow(sub.df), num.well),]
    return(data.subset)
  } else {
    return(NULL)
  }
}