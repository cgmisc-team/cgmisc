#Displays a msg 

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
                            "\n","License",License,"\n\n", Title))
  
  #install Bioconductor packages if necessary
  
  source("http://bioconductor.org/biocLite.R", echo = FALSE, verbose = FALSE)
  
  if(!require(GenomicRanges)){
    biocLite("GenomicRanges")
  }
  
  if(!require(rtracklayer)){
    biocLite("rtracklayer")
  }

}