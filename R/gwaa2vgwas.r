##' @title Convert GenABEL data structure to vGWAS data structure
##' 
##' @description Converts  GenABEL data (\code{\link[GenABEL]{gwaa.data-class}} object) to \code{\link[vGWAS]{vGWAS}} input format.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@slu.se}>
##' @param data an object of the \code{\link[GenABEL]{gwaa.data-class}} 
##' @param trait.name name of the trait of interest 
##' @return a vGWAS readable list
##' @examples
##'  \dontrun{vGWAS.input.data <- gwaa2vgwas(data=cgmiscdat1, "trait1")}
##' @keywords GenABEL, vGWAS, convert, gwaa
##' @seealso \code{\link[GenABEL]{GenABEL}}, \code{\link[vGWAS]{vGWAS}}
##' @export
gwaa2vgwas <- function(data, trait.name) {
  if (!require(vGWAS)) {
    stop("ERROR! Library vGWAS couldn't be loaded.")
  }
  chr <- as.numeric(chromosome(gtdata(data)))
  map <- as.numeric(map(gtdata(data)))
  geno <- sub("/","", as.matrix(as.character((data))))
  pheno <- phdata(data)[,trait.name]
  result <- list(chr=chr, map=map, geno=geno, pheno=pheno)
  result
}