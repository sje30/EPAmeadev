#Contains functions for well_analysis.R script

#Function that creates data frame containing file names of HDF5 files and the DIV value for each file
make.age.df<-function(mea.data.files){
  age.df<-NULL
  for (i in ages){
    if ( i<10) {
      age.chars<-paste("DIV0", i, sep="")
    } else {
      age.chars<-paste("DIV", i, sep="")
    }
    find.files<- sapply(mea.data.files[,"files"], function(x) grepl(age.chars, x))
    file.names<- names(which(find.files))
    mea.names<-mea.data.files[which(find.files),"mea.key"]
    age.df<-rbind(age.df, data.frame(mea.name=mea.names, file.name = file.names, rep(i, length(file.names)) ))
  }
  age.df
}




#Function to create s objects with spike train length of 15 minutes or less. Input is original 
#s object, and output is shortened s object.
short.s.object<-function(s){
  rec.start<-s$rec.time[2]-900
  st<-sapply(s$spikes, function(x) (short.train(x, rec.start)))
  indxs<-as.vector(which(unlist(lapply(st, is.null))))
  if(length(indxs)>0) {
    s$spikes<-st[-indxs]
    chan<-s$channels[-indxs]
    s$channels<-chan
    s$NCells<-length(chan)
  } else {
    s$spikes<-st
  }
  s$nspikes<-sapply(s$spikes, length)
  s$rec.time[1]<-rec.start
  rec.length<-s$rec.time[2]-s$rec.time[1]
  s$meanfiringrate<-sapply(s$spikes, function(x) length(x)/rec.length)
  s
}

#Function to shorten spike train. Input is spike train and start time of the shortened recording.
#Output is shortened spike train
short.train<-function(st, rec.start) {
  st<-unlist(st)
  short<-st[which(st>rec.start)]
  if (length(short)==0){
    ret<-NULL
  } else{ 
    ret<- as.vector(short)
  }
  ret
}



#Performs specified analysis. Input is a list of arrays spike times, and type ("burst", "network.spike", 
#"fir.rate", "corr" or "theta.burst"). Output is data frame of the results of the specified analysis
well.analysis<- function(s.list, type) {
  if (type=="burst") { 
    ret.list<- mclapply(s.list, function(x) burst.analysis.well(x, type="burst"), mc.cores=num.cores)
    ret.trans<-lapply(ret.list, t)
    ret.mat<-do.call("rbind", ret.trans)
    well.names<-rownames(ret.mat)
    rownames(ret.mat)<-NULL
    ret<-data.frame(ret.mat, well.names,  stringsAsFactors=FALSE)
  } else if (type=="fir.rate") {
    rate.list<- lapply(s.list, well.fr)
    rate.df<-unlist(rate.list)
    rate.df[which(is.na(rate.df))]<-0
    burst.fr<-mclapply(s.list, function(x) burst.analysis.well(x, type="fir.rate"), mc.cores=num.cores)
    burst.fr.df<-unlist(burst.fr)
    well.names<-names(rate.df)
    ret<-data.frame(as.vector(rate.df), as.vector(burst.fr.df), well.names, stringsAsFactors=FALSE)
  } else if (type=="corr"){ 
    times<-c(0.05, 0.005, 0.001)
    ci.mean.df<-NULL
    for (i in 1:length(times)){
      ci.mean<-unlist(mclapply(s.list, function(x) well.corr(x, dt=times[i]), mc.cores=num.cores))
      ci.mean.df<-cbind(ci.mean.df, ci.mean)
    }
    well.names<-rownames(ci.mean.df)
    rownames(ci.mean.df)<-NULL
    ret<-data.frame(ci.mean.df, well.names, stringsAsFactors=FALSE)
  } else if (type=="theta.burst") {
    th.burst<-unlist(mclapply(s.list, tburst.well, mc.cores=num.cores))
    well.names<-names(th.burst)
    ret<-data.frame(th.burst, well.names)
  } else if (type=="network.spike") {
    ns.list<-mclapply(s.list, function(x) per.well.ns(x, ns.T, ns.N), mc.cores=num.cores)
    well.names<-as.vector(sapply(ns.list, rownames))
    ret<- data.frame(do.call("rbind", ns.list), well.names)
    rownames(ret)<-NULL
  }
  names(ret)<-c(df.names(type), "well")
  all.plate.nums<-data.frame(1:16, sapply(file.strings[1:16], plate.id))
  plate.num.reps<-as.vector(sapply(all.plate.nums[,1], function(x) rep(x, 48)))
  plate.num.list<-c(plate.num.reps, plate.num.reps[-((7*48+1):(8*48))], plate.num.reps, plate.num.reps[-((8*48+1):(9*48))])
  DIVs<-as.vector(sapply(age.list, function(x) rep(x, 48)))
  ret<-data.frame(ret, plate.num = plate.num.list, plate.id=all.plate.nums[plate.num.list,2], DIV=DIVs)
  ret
}

#Function to calculate burst related statistics. Input is s object 
#and type of analysis - either "fir.rate" or "burst".
burst.analysis.well <- function(s, type) {
  b<-spikes.to.bursts(s, "mi")
  s$allb <- b
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  names(wells) <- wells 
  if (type=="fir.rate") {
    sapply(wells, function(x) burst.fr.well(s, x))
  } else if (type=="burst"){
    sapply(wells, function(x) burst.allwells(s, x))
  }
}

#Function to calculate within burst firing rate on a well. Input is s object 
#and the well name. Output is mean within burst firing rate on that well.
burst.fr.well<-function(s, well) {
  indexes = names.to.indexes(names(s$spikes), well, allow.na=TRUE)
  bursts<-s$allb[indexes]
  not.nul<- which(sapply(bursts, function(y) dim(y)[1])>0)
  firing.rate<-logical(0) 
  for (j in not.nul) { 
    bl<-  bursts[[j]][,"len"]/ bursts[[j]][,"durn"]
    bl.mn<-mean(bl)
    firing.rate<-append(firing.rate, bl.mn)
  }
  mean(firing.rate)
}

#Function to calculate within burst statistics (rate, duration, % spikes in bursts and CV of IBI) on a well. 
#Input is s object and the well name. Output are mean burst statistics on this well.
burst.allwells<-function(s, well){
  bsum.list<-calc.burst.summary(s)
  indexes = names.to.indexes(names(s$spikes), well, allow.na=TRUE)
  well.bursts<-bsum.list[indexes,]
  bwell<-cbind(well.bursts[,c(6,8,12,19,7)], well.bursts[,16]/well.bursts[,15])
 apply(bwell, 2, function(x) mean(x,na.rm=TRUE))
}


#Calculates median bursts per minute, burst duration, percent of spikes in bursts and 
#CV of IBI from burst summary list
burst.stats<-function(bsum.list){
  bursts.pm<-sapply(sapply(bsum.list, "[[", 6), function(x) median(x[!x==0]))
  burst.dur <- sapply(sapply(bsum.list, "[[", 8), function(x) median(x[!x==0]))
  s.in.b<- sapply(sapply(bsum.list, "[[", 12), function(x) median(x[!x==0], na.rm=TRUE))
  cv.IBI <- sapply(sapply(bsum.list, "[[", 19), function(x) median(x[!x==0], na.rm=TRUE))
  stat.sum<-data.frame(bursts.pm=bursts.pm, burst.dur=burst.dur, s.in.b=s.in.b, cv.IBI=cv.IBI)
}

#Function to calculate mean firing rate on an array. Input is s object, output is 
#mean firing rate on all wells of the array
well.fr<- function(s) {
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  names(wells) <- wells
  mfr<-s$meanfiringrate
  sapply(wells, function(x) fr.per(x, mfr,s))
}

#Function to calculate mean firing rate on a well. Input is well name and vector of mean
#firing rate on all electrodes, output is mean firing rate on the well
fr.per<-function(well, mfr, s) {
  indexes = names.to.indexes(names(s$spikes), well, allow.na=TRUE)
  mean(mfr[indexes], na.rm=TRUE)
}

#Function to calculate correlation on each well of an array. Input is s object and dt value
#to be used. Output is vector of mean correlation on all wells of the array.
well.corr<-function(s, dt) {
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  names(wells) <- wells
  corr.mat<-tiling.allpairwise(s, dt=dt)
  mean.corrs<-sapply(wells, function(x) pick.corr.well(s, corr.mat, x))
  mean.corrs
}

#Function to calculate correlation on a well. Input is s object, correlation matrix on the array
#and well name. Output is mean correlation on the well.
pick.corr.well<-function(s, corr.mat, well) {
  indexes = names.to.indexes(names(s$spikes), well, allow.na=TRUE)
  corr.sub.mat<-corr.mat[indexes, indexes]
  m<- mean(corr.sub.mat[upper.tri(corr.sub.mat)])
  m
}

#Function to calculate fraction of electrodes theta bursting on all wells of an array. Input is
#s object, output is vector of fractions of theta bursting.
tburst.well <-function(s) {
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  tb<-sapply(wells, function(x) tburst.elec(s, x))
}

#Function to calculate the fraction of electrodes theta bursting on a well. Input is s object and well name,
#output is number representing fraction of theta bursting on that well
tburst.elec<-function(s, well){
  indexes = names.to.indexes(names(s$spikes), well, allow.na=TRUE)
  spks<-s$spikes[indexes]
  t.burst<-lapply(spks, smooth.isi.elec)
  theta.frac<- sum(unlist(t.burst))/length(t.burst)
}


#Function to calculate theta bursting. Input is vector of spike times from one electrode, output is TRUE
#if theta bursting is occuring at this electrode, otherwise false
smooth.isi.elec<-function(st) {
  allisi <- diff(unlist(st))
  x <- allisi
  if (length(x)>1) {
    den <- density(log(x))
    den$x <- exp(den$x)
    p <- peaks(den$y)
    pks<-which(p)
    theta.reg<-which(den$x<=(1/4) & den$x>=(1/10))
    ret<-sum(intersect(pks, theta.reg))>0 
  }else{
    ret<-FALSE
  }
  ret
}

#Function to calculate network spikes statistics. Input is s object, time threshold and electrode 
#threshold. Output is data frame containing number, median peak value and duration of network spikes
#on each well of the array
per.well.ns<-function(s, ns.T, ns.N) {
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  names(wells) <- wells
  ns.all<-mclapply(wells, function(well) {
    compute.ns(s, ns.T=ns.T, ns.N=ns.N, whichcells=well)
  }, mc.cores=num.cores)
  ns.sum<-sapply(ns.all, function(x) rbind(x$brief))
  ns.mat<-as.matrix(ns.sum)
  ns.mat[1,which(is.na(ns.mat[1,]))]<-0
  ns.mat.ret<-data.frame(ns.rate=ns.mat[1,]/15, ns.peak=ns.mat[2,], ns.dur=ns.mat[4,])
  ns.mat.ret
}

#Specifies the column names of data-frame for each type of analysis
df.names<-function(type) {
  if (type=="burst") {
    names<-c("Burst Rate (per min)", "Burst Duration (sec)", "% Spikes in Bursts", "CV of IBI", "Bursty", "CV of ISI")
  } else if (type=="network.spike") {
    names<-c("Network spike rate (per min)", "Network spike peak", "Network spike duration (sec)")
  }else if (type=="fir.rate"){
    names<-c("Firing rate (Hz)", "Within burst firing rate (Hz)")
  } else if (type=="corr") {
    names<-c("dt = 0.05", "dt= 0.005", "dt=0.001")
  } else if (type=="theta.burst") {
    names<-c("Fraction of electrodes theta bursting")
  } 
}


#Extracts plate ID from file name
plate.id<-function(file.name) {
  str_beg=max(gregexpr("/", file.name)[[1]])+1
  str_end=max(gregexpr("_D", file.name)[[1]])-1
  substr(file.name, start=str_beg, stop=str_end)
}



















