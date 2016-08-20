FROM rocker/hadleyverse
MAINTAINER Stephen Eglen <sje30@cam.ac.uk>
RUN apt-get update
RUN Rscript -e 'install.packages(c("doBy", "ggplot2", "", "xtable"))'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("rhdf5")'
