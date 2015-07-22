#Contains functions for Control_analysis.R script
# gets channel names of channels which fall within filter range
get.ids<-function( file, params  ){

  temp.all<-h5ls(file)
  
  temp<-h5read(file, 'summary.table')
  temp.mfr<-h5read(file,"sCount")/unlist(temp[5]-temp[4])
  temp.ind<-which(params$low<=temp.mfr &
                    temp.mfr<=params$high  )
  if (is.element("names",temp.all[,2])){
    ids=h5read(file, "names")[temp.ind]
  } else if (is.element("channels", temp.all[,2])){
    ids=h5read(file, "channels")[temp.ind]
  }
  ids
}




#Creates data frame containing file names of HDF5 files and the DIV value for each file
make.age.df<-function(mea.data.files){
  age.df<-NULL
  for (i in ages){
    if ( i<10) {
      age.chars<-paste("DIV0", i, sep="")
      } else {
        age.chars<-paste("DIV", i, sep="")
      }
    find.files<- sapply(mea.data.files[,"files"], function(x) grepl(age.chars, x))
    file.names<- mea.data.files[ which(find.files), ]
    age.df<-rbind(age.df, data.frame(file.name = file.names, rep(i, length(file.names)) ))
  }
  age.df
}


sum.by.well.fr.rate<-function(x){
  x$well<-substring(x$channels,1,2)
  t<-aggregate(list(x$meanfiringrate), by=list(x$well), median, na.rm=T)
  names(t)<-c("well","meanfiringrate")
  temp<-unlist(lapply( strsplit(basename(x$file),split="_"), function(x) x[2:4]))
  t$date<-rep(temp[1], length(t$well))
  t$Plate.SN<-rep(temp[2], length(t$well))
  t$DIV<-rep(temp[3], length(t$well))
  t
}

#input 2 lists and retrieve the file data
get.meta.data<-function(s.list, ret, type){
  file.name.t<-unlist( lapply(s.list, function(x){
    strsplit( basename(x$file), split=".h5")[1] }  ) )
  
  if ( type=="burst" ) {
    n.rows=unlist( lapply(ret, function(x) dim(x)[1] ) )
  } else if ( type=="network.spike"  ){
    n.rows=unlist(lapply(ret,function(x) dim(x)[1]))
  } else if (type=="fir.rate"){
    n.rows=unlist( lapply(ret, function(x) dim(x)[1] ) )
  } else if (type=="corr"){
    temp=sum(ret$well=="A1")
    n.rows=rep( dim(ret)[1]/temp, temp)
  } else if (type=="theta.burst"){
    n.rows=unlist(lapply(ret,function(x) dim(x)[1]))
  }
  
  file.name<-rep(file.name.t, times=n.rows)
  
  duration.t<-unlist( lapply(s.list, function(x){
    dur=(x$rec.time[2]-x$rec.time[1]); dur } ) )
  duration<-rep(duration.t, times=n.rows)
  
  date<-unlist(lapply(strsplit(file.name, split="_" ),function(x) x[2] ))
  plate<-unlist(lapply(strsplit(file.name, split="_" ),function(x) x[3] ))
  div<-unlist(lapply(strsplit(file.name, split="_" ),function(x) x[4] ))
  if ( type=="burst" ) {
    temp.df<-do.call(rbind.data.frame, ret )
    df<-cbind.data.frame(file.name ,
                         date, plate, div, duration, temp.df)
  } else if ( type=="network.spike"  ){
    
    temp<-do.call(rbind.data.frame, ret )
    df<-cbind.data.frame(file.name,
                         date, plate, div, duration, temp )
  } else if (type=="fir.rate") {
    temp.df<-do.call(rbind.data.frame, ret )
    df<-cbind.data.frame(file.name ,
                         date, plate, div, duration, well=temp.df[,"well"],
                         fir.rate=temp.df[,"meanfiringrate"] )
  } else if (type=="corr"){
    df<-cbind.data.frame(file.name ,
                         date, plate, div, duration, ret )
    
  } else if (type=="theta.burst"){
    temp.df<-do.call(rbind.data.frame, ret )
    df<-cbind.data.frame(file.name ,
                         date, plate, div, duration, temp.df )
  }
  
  df
}



#Performs specified analysis. Input is a list of arrays spike times, and type ("burst", "network.spike", 
#"fir.rate", "corr" or "theta.burst"). Output is data frame of the results of the specified analysis
analysis.fn<- function(s.list, type) {
  if (type=="burst") { 
    # list of plates by well burst summaries
    ret<- burst.analysis(s.list)
    df<-get.meta.data(s.list, ret, type )   
  } else if (type=="fir.rate" ) {
    #rate.list is list of by well firing rates
    rate.list<-lapply(s.list,  sum.by.well.fr.rate )
    df2<-get.meta.data( s.list, rate.list, type )
    df2[which(is.infinite(df2$meanfiringrate))]<-0
    
    df<-df2
   
  } else if (type=="corr"){ 
    
    times<-c(0.05, 0.005, 0.001)
    ci.mean.df1<-c()
    for (i in 1:length(times) ){
      ci.mean<-unlist(lapply(s.list, 
               function(x) well.corr(x, dt=times[i])) )
      ci.mean.df1<-cbind(ci.mean.df1, ci.mean)
    }
    suppressWarnings( ci.mean.df<-cbind.data.frame(names(ci.mean), ci.mean.df1) )
    colnames(ci.mean.df)<-c("well",times)
    df<-get.meta.data(s.list, ci.mean.df, type)

  } else if (type=="theta.burst") {
    
    th.burst<-mclapply(s.list, tburst.elec, mc.cores=num.cores)
    df<-get.meta.data(s.list, th.burst, type)
    
  } else if (type=="network.spike" ) {
    ns.list<-mclapply(s.list, function(x) well.ns(x, ns.T, ns.N), mc.cores=num.cores)
    df<-get.meta.data(s.list, ns.list, type)

  }
  
  df
}

#Function to calculate within burst firing rate. Input is list of bursts, 
#output is vector of firing rates
burst.fir.rate.by.well<-function(bursts, s) {
  not.nul<- which(sapply(bursts, function(y) dim(y)[1])>0)
  ch.want<-s$channels[not.nul]
  firing.rate<-logical(0) 
  for (j in not.nul) { 
    bl<-  bursts[[j]][,"len"]/ bursts[[j]][,"durn"]
    bl.mn<-mean(bl)
    firing.rate<-append(firing.rate, bl.mn)
  }
  temp<-cbind.data.frame(firing.rate, well=substring(ch.want,1,2))
  temp2<-aggregate(temp$firing.rate, by=list(temp$well), median, na.rm=T)
  names(temp2)<-c("well", "medfiringrate" )
  
  temp3<-unlist(lapply( strsplit(basename(s$file),split="_"), function(x) x[2:4]))
  temp2$date<-rep(temp3[1], length(temp2$well))
  temp2$Plate.SN<-rep(temp3[2], length(temp2$well))
  temp2$DIV<-rep(temp3[3], length(temp2$well))
  temp2
}


#Specifies the column names of data-frame for each type of analysis
df.names<-function(type) {
  if (type=="burst") {
    names<-c("Burst Rate (per min)", "Burst Duration (sec)", "% Spikes in Bursts", "CV of IBI")
  } else if (type=="network.spike") {
    names<-c("Network spike rate (per min)", "Network spike peak", "Network spike duration (sec)")
  }else if (type=="fir.rate"){
    names<-c("Firing rate (Hz)", "Within burst firing rate (Hz)", "Standard deviation of firing rate")
  } else if (type=="corr") {
    names<-c("dt = 0.05", "dt= 0.005", "dt=0.001")
  } else if (type=="theta.burst") {
    names<-c("Fraction of electrodes theta bursting")
  } 
}

#create wells from channels
sum.by.well<-function(b){
  
  b$well<-substring( as.character(b$channels),1,2)
  b2<-aggregate(list(b), by=list(b$well), FUN=mean, na.rm=TRUE )
  b2=b2[ ,c(1,3:21)]
  names(b2)[1]<-"well"
  b2
}

#Function to calculate burst statistics
#Input is list of spike trains, ouput is data frame of burst statistics
burst.analysis <- function(s.list) {
  #names a allb for each plate, a list with entry for each channel
  b.list<-mclapply(s.list , function(x) spikes.to.bursts(x, "mi"), mc.cores=num.cores)
  
  for (j in 1:length(s.list) ) {
    s.list[[j]]$allb <- b.list[[j]]
  }
  
  #bsum.list<-lapply(s.list, calc.burst.summary) #added in cv.ISIs
  #b.list<-mclapply(s.list , function(x) spikes.to.bursts(x, "mi"), mc.cores=num.cores)
  b.list<-lapply(s.list, calc.burst.summary.2 )
  #bsum.list<-lapply(s.list, calc.burst.summary.2 ) 
  # don't sum by well
  #bsum.list<-lapply(bsum.list, sum.by.well )
  bsum.list<-lapply(bsum.list, sum.by.well )

  bsum.list
}


# calc.burst.summary but with extra CV.isis
calc.burst.summary.2<-function(s, bursty.threshold = 1){
  allb <- s$allb
  channels <- s$channels
  spikes <- as.vector(s$nspikes)
  duration <- s$rec.time[2] - s$rec.time[1]
  mean.freq <- round(spikes/duration, 3)
  nbursts <- sapply(allb, num.bursts)
  bursts.per.sec <- round(nbursts/duration, 3)
  bursts.per.min <- bursts.per.sec * 60
  bursty = ifelse(bursts.per.min >= bursty.threshold, 1, 0)
  durations <- burstinfo(allb, "durn")
  mean.dur <- round(sapply(durations, mean), 3)
  sd.dur <- round(sapply(durations, sd), 3)
  ISIs = calc.all.isi(s, allb)
  mean.ISIs = sapply(ISIs, mean)
  sd.ISIs = unlist(sapply(ISIs, sd, na.rm = TRUE))
  cv.ISIs = sd.ISIs/mean.ISIs
  ns <- burstinfo(allb, "len")
  mean.spikes <- round(sapply(ns, mean), 3)
  sd.spikes <- round(sapply(ns, sd), 3)
  total.spikes.in.burst <- sapply(ns, sum)
  per.spikes.in.burst <- round(100 * (total.spikes.in.burst/spikes), 
                               3)
  si <- burstinfo(allb, "SI")
  mean.si <- round(sapply(si, mean), 3)
  IBIs <- calc.all.ibi(s, allb)
  mean.IBIs <- sapply(IBIs, mean)
  sd.IBIs <- sapply(IBIs, sd, na.rm = TRUE)
  cv.IBIs <- round(sd.IBIs/mean.IBIs, 3)
  mean.IBIs <- round(mean.IBIs, 3)
  sd.IBIs <- round(sd.IBIs, 3)
  df <- data.frame(channels = channels, spikes = spikes, mean.freq = mean.freq, 
                   nbursts = nbursts, bursts.per.sec = bursts.per.sec, bursts.per.min = bursts.per.min, 
                   bursty = bursty, mean.dur = mean.dur, sd.dur = sd.dur, 
                   mean.spikes = mean.spikes, sd.spikes = sd.spikes, per.spikes.in.burst = per.spikes.in.burst, 
                   per.spikes.out.burst = round(100 - per.spikes.in.burst, 
                                                3), 
                   mean.si = mean.si, mean.isis = mean.ISIs, sd.mean.isis = sd.ISIs, 
                   cv.ISIs = cv.ISIs,
                   mean.IBIs = mean.IBIs, sd.IBIs = sd.IBIs, cv.IBIs = cv.IBIs)
  df
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

#Function to calculate within burst firing rate. Input is list of bursts, 
#output is vector of firing rates
burst.fir.rate<-function(bursts) {
  not.nul<- which(sapply(bursts, function(y) dim(y)[1])>0)
  firing.rate<-logical(0) 
  for (j in not.nul) { 
    bl<-  bursts[[j]][,"len"]/ bursts[[j]][,"durn"]
    bl.mn<-mean(bl)
    firing.rate<-append(firing.rate, bl.mn)
  }
  firing.rate
}
  

#Function to calculate theta bursting. Input is vector of spike times from one electrode, output is TRUE if
#theta bursting is occuring at this electrode, otherwise false
smooth.isi.elec<-function(s) {
  allisi <- diff(unlist(s))
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

#Calculates fraction of electrodes theta bursting on any array. Input is list of spike trains from one array,
#output is number representing fraction of theta bursting on that array
tburst.elec<-function(s){
  spks<-s$spikes
  t.burst<-lapply(spks, smooth.isi.elec)
  temp.well<-substring(names(t.burst),1,2)
  theta.frac<-aggregate(x= unlist(t.burst), by=list(temp.well), FUN=function(x) sum(x)/length(x)  )
  names(theta.frac)<-c("well","theta.frac")
  theta.frac
}
  



#Well-wise calculation of tiling correlation, adapted from Catherine's code
tiling.allpairwise.wells <- function(s, dt=0.05, whichcells=NULL) {
  indexes = names.to.indexes(names(s$spikes), whichcells, allow.na=TRUE)
  N <- length(indexes)
  indices=array(0,dim=c(N,N))
  
  if(N>0){
    for(i in 1:N) {
      for(j in 1:N) {
        
        z <- .C("run_TM",as.integer(length(s$spikes[indexes][[i]])),
                as.integer(length(s$spikes[indexes][[j]])),
                as.double(dt),
                as.double(s$rec.time[[2]]),
                index=as.double(1),
                as.double(as.vector(s$spikes[indexes][[i]])),
                as.double(as.vector(s$spikes[indexes][[j]])),
                PACKAGE="sjemea")
        indices[i,j] <- z[[5]]
      }
    }
  }
  indices
}

#Calculates the mean correlation per well, 
#Input is spike train and dt value, output is mean correlation value.
well.corr<-function(s, dt) {
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  names(wells) <- wells
  corr.all<-mclapply( wells, function(well) {
    tiling.allpairwise.wells(s, dt=dt, whichcells=well)
      }, mc.cores=num.cores)

  mean.corrs<-sapply(corr.all, function(m) mean(m[upper.tri(m)], na.rm=T) )
  #total.corr<-mean(mean.corrs, na.rm=TRUE)
  mean.corrs
}

#Computes network spikes on a well-wise basis. Input is spike train, time bin and
#electrode number threshold. Output is network spike rate, peak and duration.
well.ns<-function(s, ns.T, ns.N) {
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  names(wells) <- wells
  ns.all<-mclapply(wells, function(well) {
    compute.ns(s, ns.T=ns.T, ns.N=ns.N, whichcells=well)
      }, mc.cores=num.cores)
  ns.sum<-sapply(ns.all, function(x) rbind(x$brief))
  ns.mat<-as.matrix(ns.sum)
  ns.mat[1,which(is.na(ns.mat[1,]))]<-0
  # remove this function that averages across wells
  #ns.stats<-as.matrix(apply(ns.mat, 1, function(x) mean(x, na.rm=TRUE)))
  #ns.stat.ret<-c(ns.stats[1]/15, ns.stats[2], ns.stats[4])
  ns.mat2<-cbind.data.frame(colnames(ns.mat), t(ns.mat) )
  colnames(ns.mat2)<-c("well", "ns.n","peak.m", "peak.sd", "durn.m","durn.sd" )
  ns.mat2
}

#Computes the standard deviation of well-wise mean firing rate
sd.fr<-function(s){
  plateinfo <- plateinfo(s$layout$array)
  wells <- plateinfo$wells
  names(wells) <- wells
  well.nms<-substr(names(s$meanfiringrate), 1, 2)
  mn.fr<-NULL
  for (j in 1:length(wells)){
    indexes<-which(well.nms==wells[[j]])
    mn.fr<-c(mn.fr, mean(s$meanfiringrate[indexes], na.rm=TRUE))
  }
  sd.ret<-sd(mn.fr, na.rm=TRUE)
}

#Creates pdfs containing strip charts for each of the quantities. Input is data frame of quantities to be plotted
#and type of analysis. 
plot.fn<-function(data.df, type) {
  fname<-paste(type, "pdf", sep=".")
  pdf(file = fname, width = 15, height = 12)
  par(mfrow = c(2, 2))
  main<-NULL
  nms<-names(data.df)
  lnth<-length(nms)
  if (type=="corr") {
    main<-nms[1:lnth]
    nms<-rep("Mean correlation", lnth)
  }
  for (i in 1:lnth) {
    data.df2<-data.frame(n=data.df[,i], age=age.list)
    stripchart(n~age, data= data.df2 ,  vertical=TRUE, method="jitter", col="deepskyblue", pch=16,  xaxt="n", xlab = "DIV", ylab=nms[i])
    axis(1, at = seq(0.5, 5.5, 1), labels = FALSE, tick = TRUE)
    axis(1, at=seq(1, 5, 1), labels = ages, tick=FALSE)
    boxplot(n~age, data=data.df2,  outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0, add=TRUE, medcol="blue4", axes=FALSE)
  }
  dev.off()
}

