##'@title Open UCSC GenomeBrowser for selected region
##' 
##' @description \code{open.region.UCSC} redirects you to the UCSC GenomeBrowser and shows the specified region
##' a \code{\link[GenABEL]{gwaa.data-class}} object.  
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@slu.se}> , Jagoda Jablonska <\email{jagoda100jablonska@@gmail.com}>
##' @param chr number of chromosome you want to see
##' @param coords a vector of two values, if the argument is missing the entire chromosome will be taken into consideration
##' @param assembly type of assembly (camFam2 or camFam3)
##' @return null
##' @export open.region.ucsc

open.region.ucsc <- function(chr, coords, assembly = "canFam3"){
  options("scipen" = 10)  
  if(!missing(coords)){
    browseURL(url = paste('https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=',assembly,'&position=chr',chr,'%3A',coords[1],'-',coords[2],'&hgsid=199102533_HTfCx4awReCJ9z6stFytIxPVgGa2', sep = ''))
  }else{
    browseURL(url =paste("https://genome-euro.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Dog&db=",assembly,"&position=chr",chr,"&hgt.positionInput=chr7&hgt.suggestTrack=refGene&Submit=submit&hgsid=198319470_ttWacdiriYPZdD6ssreWAGjkASKc",sep=""))
  }
}
open.region.UCSC <- open.region.ucsc

