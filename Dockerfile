FROM rocker/hadleyverse
MAINTAINER Stephen Eglen <sje30@cam.ac.uk>
RUN apt-get update


ENV PROJ /home/rstudio/epameadev1
RUN mkdir $PROJ
#RUN git clone https://github.com/sje30/EPAmeadev.git $PROJ

WORKDIR $PROJ
RUN wget https://github.com/sje30/EPAmeadev/archive/master.zip
#RUN unzip master.zip
#RUN rm master.zip

RUN ls
# RUN Rscript installs.R


WORKDIR $PROJ/analysis_code
# RUN Rscript well_analysis.R
# RUN Rscript EPAmeadev_feature_plots.R

## building locally:
## docker build -t sje30/epameadev .


## Docker seems to have problems with large files:
## downloading the .zip is much smaller than the git repo, (which
## includes .git):
## https://github.com/docker/docker/issues/2292


