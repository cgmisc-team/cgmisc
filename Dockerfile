FROM rocker/r-base:4.1.3
MAINTAINER Marcin Kierczak <marcin.kierczak@scilifelab.se>

RUN apt update && apt install -y r-base r-cran-devtools libcurl4-openssl-dev
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz', type='source', repos=NULL)"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz', type='source', repos=NULL)"
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('rtracklayer')"
RUN R -e "devtools::install_github('cgmisc-team/cgmisc')"

ENTRYPOINT ["R"]
CMD ["--help"]
