<!-- badges: start -->
[![R build status](https://github.com/cgmisc-team/cgmisc/workflows/R-CMD-check/badge.svg)](https://github.com/cgmisc-team/cgmisc/actions)
<!-- badges: end -->

## Introduction

**cgmisc** is a R package that enables enhanced data analysis and visualisation of results from GWAS. The package contains several utilities and modules that complement and enhance the functionality of existing softwares. It also provides several tools for advanced visualisation of genomic data and utilises the power of the R language to aid in preparation of publication-quality figures. Some of the package functions are specific for the domestic dog (*Canis familiaris*) data.

## Release philosophy
Beginning from version 2.9.11, we are no longer using releases system. Instead, we maintain cgmisc in the CD/CI manner. From time to time, major versions will be frozen and available as source packages. Otherwise, track commit messages to know what has changed.

## Pre-requisites
`cgmisc` enchances functionalities of `GenABEL` package which is, unfortunately, no longer supported. Thus you will need to install it manualy from source available on CRAN Package Archives:

*  `install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz", type='source', repos=NULL)`  

*  `install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz", type='source', repos=NULL)`  

## Installation 
We recommend installing *cgmisc* by using:
`devtools::install_github('cgmisc-team/cgmisc')`

To install using the tarball, open a terminal and type: 
`R CMD INSTALL cgmisc_[version].tar.gz`

## How to cite `cgmisc`
Kierczak M, Jablonska J, Forsberg SKG, Bianchi M, Tengvall K, Pettersson M, Scholz V, Meadows JRS, Jern P, Carlborg O Lindblad-Toh K. cgmisc: enhanced genome-wide association analyses and visualization. Bioinformatics. Oxford University Press; 2015;31: 3830-3831. 
doi:10.1093/bioinformatics/btv426

## Selected publications using `cgmisc`

Here we list some publications where cgmisc has been helpful:

* [Bianchi M. et al. Whole-genome genotyping and resequencing reveal the association of a deletion in the complex interferon alpha gene cluster with hypothyroidism in dogs. BMC Genomics 2020.](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6700-3)

* [Weich K. et al. Pigment Intensity in Dogs is Associated with a Copy Number Variant Upstream of KITLG. Genes 2020.](https://doi.org/10.3390/genes11010075)

* [Mikkola LI et al. Novel protective and risk loci in hip dysplasia in German Shepherds. PLoS Genetics 2019.](https://doi.org/10.1371/journal.pgen.1008197)

* [Raymond B et al. Genome-wide association study for bone strength in laying hens. J Anim Sci 2018. 96 (7) 2525-2535.](https://doi.org/10.1093/jas/sky157)

* [Brown EA et al. FGF4 retrogene on CFA12 is responsible for chondrodystrophy and intervertebral disc disease in dogs. PNAS 2017. 114 (43) 11476-11481.](https://doi.org/10.1073/pnas.1709082114)

