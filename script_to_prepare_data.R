#Collates all EPA data into a RData file
# author: Ellese
# modificaitons: DIana 
#

#Performs 5 types of analysis on hdf5 files, and outputs pdf plots:
# 1. "burst" - burst analysis
# 2. "network.spike" - network spikes analysis
# 3. "fir.rate" - firing rate and within burst firing rate
# 4. "corr"  -  correlation analysis (using tiling correlation)
# 5. "theta.burst" - theta burst analysis


library(rhdf5)
library(sjemea)
library(parallel)

#set switches
mac=T

##The available number of cores
num.cores<-1 ##INPUT number of cores here
  
##The directory in which the HDF5 files are stored


  ##INPUT data directory here, in the form "/Users/Ellese/PhD/EPA/Data/"
if (mac){
  functions.file<-"/Users/dh2744/Dropbox/HistoricOntogenyAnalysis/PreparedData/functions_prepared_data.R" 
  mea.data.dir<-"/Users/dh2744/Dropbox/HistoricOntogenyAnalysis/Data/forStephenEglen/"
} else {
  drive.letter<-'C:/Users/diana_000/Dropbox'
  mea.data.dir <- paste(drive.letter,
                        '/HistoricOntogenyAnalysis/Data/forStephenEglen/',
                        sep="")
  functions.file<-paste(drive.letter,
                        "/HistoricOntogenyAnalysis/PreparedData/",
                        "functions_prepared_data.R",
                        sep="")
  
}
source( functions.file )

#Sets global variables
mea.data.files <- as.data.frame( list.files(mea.data.dir) )
colnames(mea.data.files)<-"files"; 
rownames(mea.data.files)<-mea.data.files$files

ages<-c(2,5,7,9,12)
types<-c("burst", "network.spike", "fir.rate", "corr", "theta.burst")
mi.par <- list(beg.isi=0.100, end.isi=0.250, min.ibi=0.800, min.durn=0.05, min.spikes=6)
ns.T<-0.003 #Network spike time bin
ns.N<-5 #Network spike threshold

#Reads in data files 
age.df<-make.age.df(mea.data.files)
age.list<-age.df[,2]
ind.not2<-which(age.df[,2]!=2 )
file.strings<-sapply(age.df[ind.not2 ,1] , toString)
age.list<-age.list[ind.not2]

# provide min and max filter for channels
f<-list(); f[[4]]=0
f[[1]]$low=0.001575; f[[1]]$high[1]=2.4;
f[[2]]$low=0.0081; f[[2]]$high[1]=3.3;
f[[3]]$low=0.02; f[[3]]$high[1]=3.6;
f[[4]]$low=0.0628; f[[4]]$high=4.2
names(f)<-as.character(unique(sort(age.list)))
#remove DIV 2

want.create.data=F
if (want.create.data){
  s.list<-list();
  for( i in 1:length(file.strings) ){
    age=age.list[i]
    ids<-get.ids(file.strings[i], params=f[[as.character(age)]] )
    #new 3-1-2015 data <- h5read(path.expand(h5file), name = "/")
    temp.time<-h5read(file.strings[i], name = "/" )
    rec.duration=(temp.time$summary.table$time.max-temp.time$summary.table$time.min)
    # check if recording duration greater than 15 minutes
    dur.want=60*15 # 10 minutes
    if ( rec.duration>dur.want ){
      #reduce time
      end=as.numeric( temp.time$summary.table$time.max )
      beg=as.numeric(end-dur.want)   
    } else {
      end=as.numeric( temp.time$summary.table$time.max )
      beg=as.numeric( temp.time$summary.table$time.min )
    }
    temp<-h5.read.spikes(file.strings[i], ids=ids, beg=beg, end=end )
    s.list[[i]]<-temp
  }
} else{
  if (mac){
    data.slist<-"/Users/dh2744/Dropbox/HistoricOntogenyAnalysis/PreparedData/s_list_3-2-15.RData"
    load(data.slist)
  } else {
    load(paste(drive.letter,
               "/HistoricOntogenyAnalysis/PreparedData/s_list_3-2-15.RData",
               sep="") )
  }
  
}



#Perform analysis for all 5 types of analysis. data.all stores data returned from all
#5 types of analysis
data.all<-NULL
for (i in 1:length(types)) {
  data.all[[i]]<-analysis.fn(s.list, types[i])
  #plot.fn(data.all[[i]], types[i])
}
names(data.all)<-types

total<-merge(data.all[[1]],data.all[[2]],by=c("file.name","date", "plate","div","duration","well") )
total<-merge(total, data.all[[3]], by=c("file.name","date", "plate","div","duration","well") )
total<-merge(total, data.all[[4]], by=c("file.name","date", "plate","div","duration","well") )
total<-merge(total, data.all[[5]], by=c("file.name","date", "plate","div","duration","well") )
# sort all data by total
total<-total[with(total, order(file.name, div, well)),]
#Saves data in data_all.RData
if (mac){
  save(data.all, file="/Users/dh2744/Dropbox/HistoricOntogenyAnalysis/PreparedData/data_all_wellLevel.RData")  
  save(total, file="/Users/dh2744/Dropbox/HistoricOntogenyAnalysis/PreparedData/date_all_merged.RData")  
} else {
  save(data.all, file="data_all_diana.RData")
}


