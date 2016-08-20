FROM rocker/hadleyverse
MAINTAINER Stephen Eglen <sje30@cam.ac.uk>
RUN apt-get update
RUN Rscript -e 'install.packages(c("doBy", "ggplot2", "grid", "gridExtra", "xtable"))'
RUN Rscript -e 'install.packages(c("pracma", "randomForest", "e1071"))'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("rhdf5")'
RUN Rscript -e 'devtools::install_github("sje30/sjemea",build_vignettes=TRUE)'


ENV PROJ /home/rstudio/EPAmeadev
RUN mkdir $PROJ
RUN git clone https://github.com/sje30/EPAmeadev.git $PROJ

WORKDIR $PROJ/analysis_code
RUN R CMD BATCH well_analysis.R

