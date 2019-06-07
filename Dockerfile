## We choosed the docker image "Ubuntu rolling" as base
FROM ubuntu:rolling

ARG DEBIAN_FRONTEND=noninteractive

## install linux packages (necessary for R packages used by SMAC)
RUN apt-get update
RUN apt-get -y install apt-utils
# packages important for GEOquery
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install libssl-dev
RUN apt-get -y install curl
RUN apt-get -y install r-cran-xml

RUN apt-get -y install libfftw3-dev
RUN apt-get -y install libv8-3.14-dev
RUN apt-get -y install libxml2-dev

# packages important for plotly
RUN apt-get -y install pandoc
RUN apt-get -y install imagemagick

## installing python and related libraries (useful for SMAC)
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN pip3 install argparse
RUN pip3 install python-dateutil
RUN pip3 install biopython
## adding autocomplete for python script
RUN pip3 install argcomplete
RUN activate-global-python-argcomplete --dest=- > ~/.bash_completion.d
RUN echo 'source ~/.bash_completion.d' >> ~/.bashrc
RUN /bin/bash -c 'source /etc/profile'

## Loading R environment
RUN apt-get -y install r-base

## setup R configs
# initiliasing the cran source
# RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
# installing required packages from CRAN
RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('randomcoloR')"
RUN Rscript -e "install.packages('psych')"

RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('plotly')"
RUN Rscript -e "install.packages('htmlwidgets')"
RUN Rscript -e "install.packages('heatmaply')"
RUN Rscript -e "install.packages('visNetwork')"

RUN Rscript -e "install.packages('mclust')"
RUN Rscript -e "install.packages('RankAggreg')"

RUN Rscript -e "install.packages('heatmaply')"
RUN Rscript -e "install.packages('randomcoloR')"
RUN Rscript -e "install.packages('BiocManager')"

# installing required packages from Bioconductor
RUN Rscript -e "BiocManager::install('GEOquery')"
RUN Rscript -e "BiocManager::install('genefu')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"
RUN Rscript -e "BiocManager::install('circlize')"


# installing ESTIMATE (R) package
RUN Rscript -e "install.packages('utils')" -e "library(utils)" -e "rforge <- 'http://r-forge.r-project.org'" -e "install.packages('estimate', repos=rforge, dependencies=TRUE)"

# adding SMAC scripts to the container
ADD smac smac
# changing directory to SMAC
WORKDIR "/smac"
# defining entrypoint for SMAC
ENTRYPOINT ["python3","smac.py"]
