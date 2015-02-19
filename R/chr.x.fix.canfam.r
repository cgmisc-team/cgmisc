##' @title Fix canine X chromosome naming
##' 
##' @description Function is renaming canine sex chromosome from "39" to "X". 
##' Pseudoautosomal (PAR) region (340,578 bp to 6,642,728 bp canFam3 and 314,872 bp to 6,598,016 bp canFam2) on the 
##' canine X chromosome is renamed to "parX"   
##' @param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}} 
##' @details Currently supported assemblies are: "canFam3","canFam2"
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
chr.x.fix.canfam <- function(data, assembly="canFam3") { 
  tmp <- data
  if (assembly == "canFam3") {
    region1 <- which(tmp@gtdata@chromosome==39 & tmp@gtdata@map > 340578  & tmp@gtdata@map < 6642728)
  } else if (assembly == "canFam2") {
    region1 <- which(tmp@gtdata@chromosome==39 & tmp@gtdata@map > 314872  & tmp@gtdata@map < 6598016)
  } else {
    print("ERROR! assembly not supported!")
  }
  chrom <- try(chromosome(tmp))
  chrom[region1] <- "parX"
  chrom[chrom=="39"] <- "X"
  tmp@gtdata@chromosome <- as.factor(chrom)
  tmp
}
