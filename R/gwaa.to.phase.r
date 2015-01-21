##' @title Prepare a PHASE input file from gwaa data.
##' 
##' @description Function for preparing PHASE and fastPHASE input files from \code{\link[GenABEL]{gwaa.data-class}} object. 
##' It can be used on a single entire chromosome or on a specified chromosomal region.
##' @param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}},
##' @param chr chromosome nubmber (name) to be used,
##' @param coords a vector of coordinates that define the region. If not set, the entire chromosome will be used,
##' @param outFile a complete path and name of the desired output file.
##' @return NULL
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @keywords fastPHASE PHASE haplotype genetics
##' @examples 
##' \dontrun{
##' gwaa2PHASE(data=data.qc2, chr=2, coords=c(3030587,5030587), outFile="~user/test.phase")
##' }
##' @seealso \code{\link[GenABEL]{gwaa.data-class}}
##' @export gwaa.to.phase
gwaa.to.phase <- function(data, chr, coords = NULL, outFile){
  if(!is.null(coords) & length(coords) == 2){
    region <- which(data@gtdata@chromosome == chr & data@gtdata@map > coords[1] & data@gtdata@map < coords[2])
  }else if(is.null(coords)){
    region <- which(data@gtdata@chromosome == chr)
  }else{
    region <- which(data@gtdata@chromosome == chr & data@gtdata@map %in% coords)
  }
  data.region <- data[, region]
  nSNPs <- data.region@gtdata@nsnps
  nIDs <- data.region@gtdata@nids
  map <- data.region@gtdata@map
  geno1 <- gsub(".*/", "", as.character.gwaa.data(data.region))
  geno2 <- gsub("/.*", "", as.character.gwaa.data(data.region))
  
  write(c(nIDs, nSNPs), outFile, sep="\n")
  write.table(t(c("P", map)), outFile, append=T, quote=F, row.names=F, col.names=F)
  write.table(t(rep("S", length(map))), outFile, append=T, quote=F, row.names=F, col.names=F,sep="")
  for(id in data.region@gtdata@idnames){
    write(id, outFile, append=T)
    write.table(t(geno1[id,]), outFile, append=T, row.names=F, col.names=F, quote=F, na="?")
    write.table(t(geno2[id,]), outFile, append=T, row.names=F, col.names=F, quote=F, na="?")
  }
}
gwaa.to.PHASE <- gwaa.to.phase
gwaa2phase <- gwaa.to.phase
gwaa2PHASE <- gwaa.to.phase


