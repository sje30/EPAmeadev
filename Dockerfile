FROM rocker/hadleyverse
MAINTAINER Stephen Eglen <sje30@cam.ac.uk>
RUN apt-get update


ENV PROJ /home/rstudio/epameadev
RUN mkdir $PROJ
RUN git clone https://github.com/sje30/EPAmeadev.git $PROJ

WORKDIR $PROJ
RUN ls -R
##RUN Rscript installs.R


WORKDIR $PROJ/analysis_code
## RUN R CMD BATCH well_analysis.R

