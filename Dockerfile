FROM rocker/r-base:4.1.3
MAINTAINER Marcin Kierczak <marcin.kierczak@scilifelab.se>

RUN apt update && \
    apt install -y \
    libcurl4-openssl-dev \
    libxml2-dev
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz', type='source', repos=NULL)"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz', type='source', repos=NULL)"
RUN R -e "install.packages(c('BiocManager' ,'remotes'))"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('rtracklayer')"
RUN R -e "remotes::install_github('cgmisc-team/cgmisc')"

ENTRYPOINT ["R"]
CMD ["--help"]
