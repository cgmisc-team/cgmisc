##' @title Prepare a PHASE input file from gwaa data
##' 
##' @description Function for preparing PHASE and fastPHASE input files from \code{\link[GenABEL]{gwaa.data-class}} object. 
##' It can be used on a single entire chromosome or on a specified chromosomal region.
##' @param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}},
##' @param chr chromosome number (name) to be used,
##' @param coords a vector of two coordinates that define the region. If not set, the entire chromosome will be used,
##' @param outFile a complete path and name of output file.
##' @return NULL
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @keywords fastPHASE PHASE haplotype genetics
##' @examples 
##' \dontrun{
##' create.Haploview.info(data=data.qc2, chr=2, coords=c(3030587,5030587), outFile="~user/test.info")
##' }
##' @seealso \code{\link[GenABEL]{gwaa.data-class}}
##' @export create.haploview.info create.Haploview.info
##' @aliases create.Haploview.info

create.haploview.info <- function(data, chr, coords, outFile) {
  #data  <- data@gtdata???
  if(!missing(coords)){
    region <- which(data@gtdata@chromosome == chr & map(data) > coords[1] & map(data) < coords[2])
  }
  else{
    region <- which(data@gtdata@chromosome == chr)
  }
  data.region <- data[, region]
  result <- data.frame(marker=data.region@gtdata@snpnames, coord=data.region@gtdata@map)
  write.table(result, outFile, append=T, row.names=F, col.names=F, quote=F, na="?", sep="\t")
}

create.Haploview.info <- create.haploview.info
