#Performs well-wise data analysis
#R CMD BATCH well_analysis.R

#Performs 5 types of analysis on hdf5 files, and outputs pdf plots:
# 1. "burst" - burst analysis
# 2. "network.spike" - network spikes analysis
# 3. "fir.rate" - firing rate and within burst firing rate
# 4. "corr"  -  correlation analysis (using tiling correlation)
# 5. "theta.burst" - theta burst analysis

library(rhdf5)
library(sjemea)
library(parallel)
library(ggplot2)

##PARAMETERS TO CHANGE
##Change this to the available number of cores
num.cores<-4

##Change this to the directory in which the HDF5 files are stored
mea.data.dir <- "../allH5Files"

source("well_analysis_functions.R")

mea.data.files <- make.meafile.cache(mea.data.dir)
ages<-c(5,7,9,12)
types<-c("burst", "network.spike", "fir.rate", "corr")
mi.par <- list(beg.isi=0.100, end.isi=0.250, min.ibi=0.800, min.durn=0.05, min.spikes=6)
ns.T<-0.003 #Network spike time bin
ns.N<-5 #Network spike threshold

#Reads in data files 
age.df<-make.age.df(mea.data.files)
#Remove outlying plate, which is the 24th plate
age.df<-age.df[-24,]
age.list<-age.df[,3]
file.strings<-sapply(age.df[,2], toString)

#Creates list of s objects
s.list<-mclapply(file.strings, h5.read.spikes, mc.cores=num.cores)
#Shortens all spike trains to 15 minutes of less
s.list.short<-mclapply(s.list, short.s.object, mc.cores=num.cores)

#data.all.wells is a list of data frames. Each element is one type of feature. Each data frame
#contains the data values for each well, the well name, plate number and DIV
data.all.wells<-NULL
for (i in 1:length(types)) {
  data.all.wells[[i]]<-well.analysis(s.list.short, types[i])
}
names(data.all.wells)<-types

data.all<-cbind(data.all.wells[[3]][,c(1:2)], data.all.wells[[1]][,c(5,1:4,6)], data.all.wells[[2]][,1:3], data.all.wells[[4]][, c(1, 4,6,7)])
names(data.all)<-c("Firing.rate", "Within.burst.firing.rate", "Bursting.electrodes", "Burst.rate", "Burst.duration", "Spikes.in.bursts", "CV.of.IBI", "CV.of.within.burst.ISI", "Network.spike.rate", "Network.spike.peak", "Network.spike.duration", "Correlation", "Well", "Plate.id", "DIV" )
save(data.all, file="data_wellwise.RData")
