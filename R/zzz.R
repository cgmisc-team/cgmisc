#Displays a msg 
options(warn=-1)

.onAttach <- function(lib,pkg = "cgmisc") {
  packageStartupMessage("\n Package cgmisc contains miscellaneous functions, useful for extending
genome-wide association study (GWAS) analyses. \n")
  
  Desc <- packageDescription(pkg)
  Version <- Desc$Version
  Date <- Desc$Date
  Name <- Desc$Package
  Author <-Desc$Author
  Title <- Desc$Description
  License <- Desc$License
  
  packageStartupMessage(paste("Package Name:",Name,
                            "\n","Version:", Version,
                            "\n", "Date:",Date,
                            "\n","Author:",Author,
                            "\n","License",License,"\n\n", Title,"\n"))
  
  #install Bioconductor packages if necessary
  
  source("http://bioconductor.org/biocLite.R", echo = FALSE, verbose = FALSE)
  
  if(!require(GenomicRanges)){
    biocLite("GenomicRanges")
  }
  
  if(!require(rtracklayer)){
    biocLite("rtracklayer")
  }
  
  #chceck versions of dependencies
  
  packages <- c(GenABEL='1.7.4', ggplot2='1.0.1', grid='3.1.1', genetics='1.3.8.1')
  inst.pkgs <- installed.packages()[,3]
  installed <- sort(inst.pkgs[which(names(inst.pkgs) %in% names(packages))])
  required <- sort(packages[names(installed)])
  a <- installed[which(sort(installed)!=sort(required))]
  b <- required[which(sort(installed)!=sort(required))]

  packageStartupMessage(
  if(length(a)!=0){
    paste("\nWARNING! The following package version differs from the recommendend one:", names(a),a, "->", names(b),b,"\n")
  })
}