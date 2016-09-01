#Contains functions for plotting data from dataframe data.all. Latest data available at https://github.com/sje30/EPAmeadev/blob/master/data_wellwise_14Aug.RData
library(doBy)
library(grid)
library(gridExtra)
library(ggplot2)

## Load in the data file.  This datafile comes from ./well_analysis.R
load("data_wellwise.RData")
##load("../data_wellwise_14Aug.RData")

plate.vals<-summaryBy(.~Plate.id+DIV, data = data.all, FUN = function(x) { m = median(x, na.rm=TRUE)}) 
names(plate.vals)<-c("plate.num", names(data.all)[c(15, 1:12)])
plot.names<-c("Mean firing rate (Hz)", "Within burst firing rate (Hz)", "Fraction of bursting electrodes", "Burst rate (/min)", "Burst duration (s)", "% spikes in bursts", "CV of IBI", "CV of within burst ISI", "Network spike rate (/min)", "Network spike peak","Network spike duration (s)", "Mean correlation" )
ages<-unique(data.all$DIV)

#Function to define plot colours
colour.hue.gg<-function(n, lum=70, chrom=50) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=lum, c=chrom)[1:n]
}

plot.col<-colour.hue.gg(5)[-c(2)] #Plot colours

#Calculates statistical significance between the different ages pairwise
stat.test<-function(data.df) {
  p.val<-sapply(ages[1:3], function(x) wilcox.fn(data.df, x))
  q.val<-p.adjust(p.val, method="fdr")
  star<-rep("",3)
  star[which(q.val<=0.05)]<-"*"
  star
}

#Performs Wilcoxon signed rank test
wilcox.fn <- function(data.df, age1) {
  age1.indx<-which(ages==age1)
  p.value<-wilcox.test(val~DIV, data=data.df, subset = (DIV==age1 | DIV==ages[age1.indx+1]))$p.value
}

#Function to create individual box plots for each feature with significance stars. Input is the plate.vals dataframe, the column number of the feature to the plotted and optional y limits
epa.plots<-function(plate.vals, col.num, ylims=NULL){
  plot.df<-data.frame(DIV=as.factor(plate.vals[,2]), val=plate.vals[,col.num])
  star<-stat.test(plot.df)
  plt<-qplot(DIV, val, data = plot.df, geom="boxplot", fill=DIV) + ylab(plot.names[col.num-2])+xlab("") +theme_bw()+theme(axis.line = element_line(colour = "black"),axis.title.y=element_text(size=30, vjust=0.8)) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=30)) +theme(legend.position="none")+scale_fill_manual(values = plot.col)+ylim(ifelse(rep(is.null(ylims),2), c(min(plot.df[,2]), max(plot.df[,2])), ylims)) +annotate("text",x=1.5,y=ylims[2]*0.99,label=star[1], size=14)+ annotate("text",x=2.5,y=ylims[2]*0.99,label=star[2], size=14)+ annotate("text",x=3.5,y=ylims[2]*0.99,label=star[3], size=14)
  for (i in which(star!="")){
    plt<-plt+annotate(x=c(i+0.01,i+0.01,i+0.99,i+0.99),y=c(ylims[2]-0.04*(ylims[2]-ylims[1]),ylims[2]-0.02*(ylims[2]-ylims[1]),ylims[2]-0.02*(ylims[2]-ylims[1]),ylims[2]-0.04*(ylims[2]-ylims[1])),"path")
  }
  plt
}

#Twelve panel figure containing box plots for each feature
tiff("feat_all.tiff", height=2200, width=2750, res = 70)
grid.arrange(epa.plots(plate.vals, 3, ylims=c(0,2))+theme(plot.margin=unit(c(0.3,1.5,1,1.5), "cm")), epa.plots(plate.vals, 6, c(0,5.5))+theme(plot.margin=unit(c(0.3,1.5,1,1.5), "cm")), epa.plots(plate.vals, 7, c(0, 1.7))+theme(plot.margin=unit(c(0.3,1.5,1,1.5), "cm")), epa.plots(plate.vals, 5, c(0,1.1))+theme(plot.margin=unit(c(0.5,1.5,1,1.5), "cm")), epa.plots(plate.vals, 4, c(0, 135))+theme(plot.margin=unit(c(0.5,1.5,1,1.5), "cm")), epa.plots(plate.vals, 8, c(0,100))+theme(plot.margin=unit(c(0.5,1.5,1,1.5), "cm")), epa.plots(plate.vals, 9, c(0.3, 1.6))+theme(plot.margin=unit(c(0.5,1.5,1,1.5), "cm")), epa.plots(plate.vals, 10, c(0.8, 3))+theme(plot.margin=unit(c(0.5,1.5,1,1.5), "cm")), epa.plots(plate.vals, 11, c(0,10))+theme(plot.margin=unit(c(0.5,1.5,1,1.5), "cm")), epa.plots(plate.vals, 13, c(0, 0.038))+theme(plot.margin=unit(c(0.5,1.5,1.5,1.5), "cm")), epa.plots(plate.vals, 12, c(4, 12))+theme(plot.margin=unit(c(0.5,1.5,1.5,1.5), "cm")), epa.plots(plate.vals, 14, c(0,1))+theme(plot.margin=unit(c(0.5,1.5,1.5,1.5), "cm")), heights= rep(unit(0.24, "npc"),4), widths = rep(unit(0.32, "npc"),3) , nrow=4, ncol=3)
grid.text(LETTERS[1:12], gp=gpar(fontsize=50),
          x=unit(rep(c(0.03, 0.355, 0.675),4), "npc"),
          y=unit(c(rep(0.985,3), rep(0.74, 3), rep(0.5, 3), rep(0.26,3)), "npc"))
grid.text("Days in Vitro", x = unit(0.5, "npc"), y = unit(0.03, "npc"), just = "centre", gp=gpar(fontsize=40))
dev.off()
