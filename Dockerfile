## We choosed the docker image "Ubuntu rolling" as base
FROM ubuntu:rolling
ARG DEBIAN_FRONTEND=noninteractive
## Install required packages
RUN echo 'updating package list'
RUN apt-get -y update
RUN apt-get -y upgrade
RUN echo 'Installing linux packages >> ESSENTIAL'
RUN apt-get -y install apt-utils
RUN apt-get -y install wget
## Packages for plotly
RUN apt-get -y install pandoc
RUN apt-get -y install imagemagick
## Packages for GEOparse
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install libssl-dev
RUN apt-get -y install libxml2-dev
RUN apt-get -y install curl

RUN echo 'Installing python >> BASE'
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN echo 'Installing python >> PACKAGES'
RUN pip3 install argparse
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install GEOparse
RUN pip3 install tenacity
RUN pip3 install biopython
RUN pip3 install python-dateutil
RUN pip3 install func-timeout
## adding autocomplete for python script
RUN pip3 install argcomplete
RUN activate-global-python-argcomplete --dest=- > ~/.bash_completion.d
RUN echo 'source ~/.bash_completion.d' >> ~/.bashrc
RUN /bin/bash -c 'source /etc/profile'
RUN echo 'Installing R >> BASE'
RUN apt-get -y install r-base
RUN echo 'Installing R >> PACKAGES'
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
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

## installing required packages from Bioconductor"
RUN Rscript -e "BiocManager::install('genefu')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"
RUN Rscript -e "BiocManager::install('circlize')"

## Packages for functional enrichment analysis
RUN Rscript -e "BiocManager::install('clusterProfiler')"
RUN Rscript -e "BiocManager::install('ReactomePA')"
RUN Rscript -e "BiocManager::install('rWikiPathways')"
RUN Rscript -e "BiocManager::install('DOSE')"
RUN Rscript -e "BiocManager::install('GSEABase')"
RUN Rscript -e "install.packages('magrittr')"
RUN Rscript -e "install.packages('msigdbr')"
RUN Rscript -e "BiocManager::install('enrichplot')"
RUN Rscript -e "BiocManager::install('org.Hs.eg.db')"

## installing ESTIMATE (R) package
RUN Rscript -e "install.packages('utils')" -e "library(utils)" -e "rforge <- 'http://r-forge.r-project.org'" -e "install.packages('estimate', repos=rforge, dependencies=TRUE)"

# adding SMAC scripts to the container
ADD SMAC2 SMAC2
# changing directory to SMAC
WORKDIR "/SMAC2"
# defining entrypoint for SMAC
ENTRYPOINT ["python3","main.py"]
