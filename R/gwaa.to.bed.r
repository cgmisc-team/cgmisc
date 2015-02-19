##'@title Create a BED file with GWAS p-values in a region
##' 
##' @description Returns a file containing GWAS p-values for a region in bed format. It can be further used as
##' UCSC Genome Browser custom track.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param gwas a \code{\link[GenABEL]{scan.gwaa-class}} object produced by GenABEL
##' @param chr chromosome to be considered
##' @param range a vector of genomic ranges
##' @param fname name of the output file. By default: output.bed 
##' @return a bed file with -log10 of corrected p-value ("Pc1df") 
##' @examples
##'  \dontrun{
##'  chr <- 36
##'  range <- c(14e6,24e6)
##'  gwas <- data.qc1
##'  gwas.region.bed(chr, range, gwas, fname="my_output.bed")
##'  }
##' @keywords bed, gwas, p-values
##' @export gwaa.to.bed
gwaa.to.bed  <- function(chr, range, gwas, fname="output.bed") {
  range <- format(range, scientific=F)
  cat(paste("browser position chr",chr,":",range[1],"-",range[2], sep=""), file=fname, append=F, sep="\n")
  cat("track type=bedGraph color=0,128,0 visiblity=full name='-log10(pval)' description='GWAS result' visibility=2", file=fname, append=T, sep="\n")
  tmp <- which(gwas@annotation$Chromosome == chr & gwas@annotation$Position >= range[1] & gwas@annotation$Position <= range[2])
  for (i in tmp) {
    pos <- gwas@annotation$Position[i]
    pval <- -log10(gwas@results$Pc1df[i])
    t <- paste("chr",chr," ", pos, " ", pos, " ", pval, sep="")
    cat(t, file=fname, append=T, sep="\n")
  }
}
gwaa2bed <- gwaa.to.bed