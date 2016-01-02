##' @title Fix sex chromosome naming
##' 
##' @description Function is renaming sex chromosome to "X" and pseudoautosomal region (PAR) to "parX".
##' One can provide the PAR of interest by setting chr, start.coord and stop.coors parameters. The canine genome 
##' (canFam2 and canFam3) is already implemented and can be accessed by setting the "assembly" parameter to either
##' "canFam2" or "canFam3" depending on the preference.
##' Canine seudoautosomal region is 340,578 bp to 6,642,728 bp in canFam3 and 314,872 bp to 6,598,016 bp in canFam2.
##'   
##' @param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}} 
##' @details Currently supported canine assemblies are: "canFam3","canFam2"
##' @return a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}}
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>, Matteo Bianchi <\email{Matteo.Bianchi@@imbim.uu.se}>
##' @keywords canine X chromosome sex pseudoautosomal PAR
##' @return an object of \code{\link[GenABEL]{gwaa.data-class}} class
##' @examples 
##' \dontrun{
##' fixed.data <- Xfix.canfam(data, assembly='canFam3')
##' }
##' @seealso \code{\link[GenABEL]{gwaa.data-class}}
##' @export chr.x.fix.canfam
chr.x.fix.canfam <- function(data, assembly="", chr=39, start.coord=0, stop.coord=0) { 
  tmp <- data
  if (assembly == "canFam3") {
    region1 <- which(tmp@gtdata@chromosome==chr & tmp@gtdata@map > 340578  & tmp@gtdata@map < 6642728)
  } else if (assembly == "canFam2") {
    region1 <- which(tmp@gtdata@chromosome==chr & tmp@gtdata@map > 314872  & tmp@gtdata@map < 6598016)
  } else {
    #print("ERROR! assembly not supported!")
    region1 <- which(tmp@gtdata@chromosome==chr & tmp@gtdata@map > start.coord  & tmp@gtdata@map < stop.coord)
  }
  chrom <- try(chromosome(tmp))
  chrom[region1] <- "parX"
  chrom[chrom==as.character(chr)] <- "X"
  tmp@gtdata@chromosome <- as.factor(chrom)
  tmp
}
